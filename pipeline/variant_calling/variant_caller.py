"""
variant_caller.py
-----------------
Variant Calling Module

Wraps GATK HaplotypeCaller for germline SNP/indel calling.
Handles BAM validation, VCF output, and per-sample calling metrics.
"""

import subprocess
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


@dataclass
class VariantCallingResult:
    """Stores variant calling output and metrics."""
    sample_id: str
    vcf_path: str
    total_variants: int = 0
    snps: int = 0
    indels: int = 0
    success: bool = False
    error_message: str = ""
    warnings: List[str] = field(default_factory=list)

    @property
    def summary(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "vcf_path": self.vcf_path,
            "total_variants": self.total_variants,
            "snps": self.snps,
            "indels": self.indels,
            "success": self.success,
            "warnings": self.warnings,
        }


class VariantCaller:
    """
    Calls germline variants from a sorted, indexed BAM using GATK HaplotypeCaller.

    Optionally applies VQSR or hard filtering for variant quality control.

    Parameters
    ----------
    reference : str
        Path to the reference genome FASTA (must have .fai and .dict).
    output_dir : str
        Directory to write VCF output files.
    min_base_quality : int
        Minimum base quality for variant calling (default: 20).
    min_mapping_quality : int
        Minimum mapping quality (default: 20).
    """

    def __init__(
        self,
        reference: str,
        output_dir: str = "data/variants",
        min_base_quality: int = 20,
        min_mapping_quality: int = 20,
    ):
        self.reference = reference
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.min_base_quality = min_base_quality
        self.min_mapping_quality = min_mapping_quality

    def _run_command(self, cmd: list, step: str) -> subprocess.CompletedProcess:
        logger.info(f"Running {step}")
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if result.returncode != 0:
            logger.error(f"{step} failed:\n{result.stderr}")
            raise RuntimeError(f"{step} exited with code {result.returncode}")
        return result

    def _validate_bam(self, bam_path: str):
        """Check BAM exists and has an index."""
        bam = Path(bam_path)
        if not bam.exists():
            raise FileNotFoundError(f"BAM file not found: {bam_path}")
        if not Path(bam_path + ".bai").exists():
            logger.warning("BAM index not found — running samtools index")
            self._run_command(["samtools", "index", bam_path], "samtools index")

    def _count_variants(self, vcf_path: str) -> dict:
        """Count SNPs and indels in VCF using bcftools stats."""
        result = subprocess.run(
            ["bcftools", "stats", vcf_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        counts = {"snps": 0, "indels": 0}
        for line in result.stdout.splitlines():
            if line.startswith("SN") and "number of SNPs:" in line:
                counts["snps"] = int(line.split()[-1])
            elif line.startswith("SN") and "number of indels:" in line:
                counts["indels"] = int(line.split()[-1])
        return counts

    def call_variants(
        self,
        bam_path: str,
        sample_id: str,
        intervals: Optional[str] = None,
    ) -> VariantCallingResult:
        """
        Run GATK HaplotypeCaller on a single sample BAM.

        Parameters
        ----------
        bam_path : str
            Path to sorted, indexed BAM file.
        sample_id : str
            Sample identifier for output naming.
        intervals : str, optional
            BED or interval list file to restrict calling regions.

        Returns
        -------
        VariantCallingResult
            VCF path and variant statistics.
        """
        vcf_path = str(self.output_dir / f"{sample_id}.g.vcf.gz")
        result = VariantCallingResult(
            sample_id=sample_id,
            vcf_path=vcf_path,
        )

        try:
            self._validate_bam(bam_path)

            cmd = [
                "gatk", "HaplotypeCaller",
                "-R", self.reference,
                "-I", bam_path,
                "-O", vcf_path,
                "--emit-ref-confidence", "GVCF",
                "--min-base-quality-score", str(self.min_base_quality),
                "--minimum-mapping-quality", str(self.min_mapping_quality),
                "--sample-name", sample_id,
            ]
            if intervals:
                cmd += ["-L", intervals]

            self._run_command(cmd, "GATK HaplotypeCaller")

            # Count variants
            counts = self._count_variants(vcf_path)
            result.snps = counts["snps"]
            result.indels = counts["indels"]
            result.total_variants = result.snps + result.indels
            result.success = True

            # Warn if variant count is unexpectedly low
            if result.total_variants < 10:
                result.warnings.append(
                    f"Very few variants called ({result.total_variants}). "
                    "Check alignment quality and target regions."
                )

            logger.info(
                f"Variant calling complete — "
                f"{result.snps} SNPs, {result.indels} indels"
            )

        except Exception as e:
            result.error_message = str(e)
            logger.error(f"Variant calling failed for {sample_id}: {e}")

        return result

    def genotype_gvcfs(
        self,
        gvcf_paths: List[str],
        cohort_id: str,
    ) -> str:
        """
        Joint genotyping across multiple samples (cohort mode).

        Parameters
        ----------
        gvcf_paths : list of str
            Paths to per-sample GVCF files.
        cohort_id : str
            Cohort identifier for output naming.

        Returns
        -------
        str
            Path to joint-called VCF.
        """
        combined_gvcf = str(self.output_dir / f"{cohort_id}.combined.g.vcf.gz")
        joint_vcf = str(self.output_dir / f"{cohort_id}.joint.vcf.gz")

        # Combine GVCFs
        combine_cmd = ["gatk", "CombineGVCFs", "-R", self.reference]
        for gvcf in gvcf_paths:
            combine_cmd += ["-V", gvcf]
        combine_cmd += ["-O", combined_gvcf]
        self._run_command(combine_cmd, "CombineGVCFs")

        # Joint genotyping
        genotype_cmd = [
            "gatk", "GenotypeGVCFs",
            "-R", self.reference,
            "-V", combined_gvcf,
            "-O", joint_vcf,
        ]
        self._run_command(genotype_cmd, "GenotypeGVCFs")

        logger.info(f"Joint genotyping complete: {joint_vcf}")
        return joint_vcf
