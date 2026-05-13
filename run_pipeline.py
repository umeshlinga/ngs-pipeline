"""
run_pipeline.py
---------------
NGS Pipeline Orchestrator

End-to-end pipeline: FASTQ QC → Alignment → Variant Calling → Annotation.
Reads sample sheet CSV and processes each sample sequentially or in parallel.

Usage
-----
    python run_pipeline.py --config config/config.yaml --samples sample_sheet.csv
    python run_pipeline.py --config config/config.yaml --samples sample_sheet.csv --threads 8
"""

import argparse
import csv
import json
import logging
import sys
import time
import yaml
from datetime import datetime
from pathlib import Path
from typing import List, Dict

from pipeline.qc.fastq_qc import FastqQC, QCMetrics
from pipeline.alignment.aligner import SequenceAligner
from pipeline.variant_calling.variant_caller import VariantCaller
from pipeline.annotation.annotator import VCFAnnotator

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(f"pipeline_run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
    ],
)
logger = logging.getLogger(__name__)


def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_sample_sheet(sample_sheet_path: str) -> List[Dict]:
    """
    Load samples from a CSV file.

    Expected columns: sample_id, fastq_r1, fastq_r2 (optional)
    """
    samples = []
    with open(sample_sheet_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            samples.append(row)
    logger.info(f"Loaded {len(samples)} samples from {sample_sheet_path}")
    return samples


def run_sample(sample: Dict, config: dict) -> Dict:
    """
    Run the full NGS pipeline for a single sample.

    Parameters
    ----------
    sample : dict
        Sample metadata with keys: sample_id, fastq_r1, fastq_r2.
    config : dict
        Pipeline configuration.

    Returns
    -------
    dict
        Per-sample results summary.
    """
    sample_id = sample["sample_id"]
    fastq_r1 = sample["fastq_r1"]
    fastq_r2 = sample.get("fastq_r2")

    logger.info(f"{'='*60}")
    logger.info(f"Processing sample: {sample_id}")
    logger.info(f"{'='*60}")

    start_time = time.time()
    results = {"sample_id": sample_id, "status": "failed", "steps": {}}

    # ── Step 1: FASTQ QC ──────────────────────────────────────────
    logger.info(f"[{sample_id}] Step 1/4: FASTQ Quality Control")
    try:
        qc = FastqQC(
            min_quality=config["qc"]["min_quality"],
            min_length=config["qc"]["min_length"],
        )
        qc_result = qc.run(fastq_r1, sample_id=sample_id)
        results["steps"]["qc"] = qc_result.to_dict()

        if qc_result.flagged:
            logger.warning(f"[{sample_id}] QC flagged: {qc_result.warnings}")
            if config["qc"].get("fail_on_flag", False):
                results["status"] = "qc_failed"
                return results
    except Exception as e:
        logger.error(f"[{sample_id}] QC failed: {e}")
        results["steps"]["qc"] = {"error": str(e)}
        return results

    # ── Step 2: Alignment ─────────────────────────────────────────
    logger.info(f"[{sample_id}] Step 2/4: Sequence Alignment")
    try:
        aligner = SequenceAligner(
            reference=config["reference"]["genome"],
            aligner=config["alignment"]["tool"],
            threads=config["alignment"]["threads"],
            output_dir=config["output"]["aligned_dir"],
        )
        align_result = aligner.align(fastq_r1, sample_id, fastq_r2=fastq_r2)
        results["steps"]["alignment"] = align_result.summary

        if not align_result.success:
            logger.error(f"[{sample_id}] Alignment failed")
            return results

        if align_result.mapping_rate < config["alignment"]["min_mapping_rate"]:
            logger.warning(
                f"[{sample_id}] Low mapping rate: {align_result.mapping_rate*100:.1f}%"
            )
    except Exception as e:
        logger.error(f"[{sample_id}] Alignment error: {e}")
        results["steps"]["alignment"] = {"error": str(e)}
        return results

    # ── Step 3: Variant Calling ───────────────────────────────────
    logger.info(f"[{sample_id}] Step 3/4: Variant Calling")
    try:
        caller = VariantCaller(
            reference=config["reference"]["genome"],
            output_dir=config["output"]["variants_dir"],
            min_base_quality=config["variant_calling"]["min_base_quality"],
        )
        vc_result = caller.call_variants(
            bam_path=align_result.bam_path,
            sample_id=sample_id,
            intervals=config["variant_calling"].get("intervals"),
        )
        results["steps"]["variant_calling"] = vc_result.summary

        if not vc_result.success:
            logger.error(f"[{sample_id}] Variant calling failed")
            return results
    except Exception as e:
        logger.error(f"[{sample_id}] Variant calling error: {e}")
        results["steps"]["variant_calling"] = {"error": str(e)}
        return results

    # ── Step 4: Annotation ────────────────────────────────────────
    logger.info(f"[{sample_id}] Step 4/4: Variant Annotation")
    try:
        annotator = VCFAnnotator(
            output_dir=config["output"]["variants_dir"],
            genome_build=config["reference"]["build"],
        )
        annotated_vcf = annotator.annotate_vcf_local(
            vcf_path=vc_result.vcf_path,
            sample_id=sample_id,
        )
        results["steps"]["annotation"] = {
            "annotated_vcf": annotated_vcf,
            "success": True,
        }
    except Exception as e:
        logger.error(f"[{sample_id}] Annotation error: {e}")
        results["steps"]["annotation"] = {"error": str(e)}

    elapsed = round(time.time() - start_time, 2)
    results["status"] = "completed"
    results["elapsed_seconds"] = elapsed
    logger.info(f"[{sample_id}] Completed in {elapsed}s")
    return results


def main():
    parser = argparse.ArgumentParser(description="NGS End-to-End Pipeline")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--samples", required=True, help="Path to sample sheet CSV")
    parser.add_argument("--threads", type=int, default=None, help="Override thread count")
    args = parser.parse_args()

    config = load_config(args.config)
    if args.threads:
        config["alignment"]["threads"] = args.threads

    samples = load_sample_sheet(args.samples)
    all_results = []

    logger.info(f"Starting NGS pipeline — {len(samples)} sample(s)")
    pipeline_start = time.time()

    for sample in samples:
        result = run_sample(sample, config)
        all_results.append(result)

    # Write run summary
    summary_path = Path(config["output"].get("summary_dir", ".")) / "pipeline_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w") as f:
        json.dump({
            "run_date": datetime.now().isoformat(),
            "total_samples": len(samples),
            "completed": sum(1 for r in all_results if r["status"] == "completed"),
            "failed": sum(1 for r in all_results if r["status"] == "failed"),
            "total_elapsed_seconds": round(time.time() - pipeline_start, 2),
            "samples": all_results,
        }, f, indent=2)

    logger.info(f"Pipeline complete. Summary: {summary_path}")
    completed = sum(1 for r in all_results if r["status"] == "completed")
    logger.info(f"Results: {completed}/{len(samples)} samples completed successfully")


if __name__ == "__main__":
    main()
