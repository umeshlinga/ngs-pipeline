"""
annotator.py
------------
VCF Annotation Module

Annotates variants using Ensembl VEP via REST API and local VEP tool.
Extracts gene names, consequences, allele frequencies, and clinical significance.
"""

import json
import logging
import requests
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Dict

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ENSEMBL_REST_URL = "https://rest.ensembl.org/vep/human/hgvs"


@dataclass
class VariantAnnotation:
    """Stores annotation for a single variant."""
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str = ""
    consequence: str = ""
    impact: str = ""
    allele_frequency: float = 0.0
    clinvar_significance: str = ""
    transcript_id: str = ""
    protein_change: str = ""
    raw: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "gene": self.gene,
            "consequence": self.consequence,
            "impact": self.impact,
            "allele_frequency": self.allele_frequency,
            "clinvar_significance": self.clinvar_significance,
            "transcript_id": self.transcript_id,
            "protein_change": self.protein_change,
        }


class VCFAnnotator:
    """
    Annotates VCF files using Ensembl VEP.

    Supports both REST API annotation (small variant sets)
    and local VEP tool annotation (large cohort VCFs).

    Parameters
    ----------
    output_dir : str
        Directory to write annotated output files.
    genome_build : str
        Reference genome build (default: 'GRCh38').
    use_rest_api : bool
        Use Ensembl REST API for annotation (default: False).
    """

    def __init__(
        self,
        output_dir: str = "data/variants",
        genome_build: str = "GRCh38",
        use_rest_api: bool = False,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.genome_build = genome_build
        self.use_rest_api = use_rest_api

    def _parse_vcf_variants(self, vcf_path: str) -> List[Dict]:
        """Extract variant records from a VCF file."""
        variants = []
        with open(vcf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                variants.append({
                    "chrom": parts[0],
                    "pos": int(parts[1]),
                    "id": parts[2],
                    "ref": parts[3],
                    "alt": parts[4],
                    "qual": parts[5] if len(parts) > 5 else ".",
                    "filter": parts[6] if len(parts) > 6 else ".",
                    "info": parts[7] if len(parts) > 7 else ".",
                })
        return variants

    def _annotate_via_rest(self, hgvs_notation: str) -> dict:
        """Query Ensembl REST API for a single variant."""
        url = f"{ENSEMBL_REST_URL}/{hgvs_notation}"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        try:
            response = requests.get(url, headers=headers, timeout=10)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            logger.warning(f"REST API annotation failed for {hgvs_notation}: {e}")
            return {}

    def _parse_vep_annotation(self, vep_result: dict, variant: dict) -> VariantAnnotation:
        """Parse VEP REST API response into VariantAnnotation."""
        annotation = VariantAnnotation(
            chrom=variant["chrom"],
            pos=variant["pos"],
            ref=variant["ref"],
            alt=variant["alt"],
            raw=vep_result,
        )
        if not vep_result:
            return annotation

        # Extract most severe consequence
        if isinstance(vep_result, list) and vep_result:
            top = vep_result[0]
            annotation.consequence = top.get("most_severe_consequence", "")

            # Extract transcript consequences
            tc = top.get("transcript_consequences", [])
            if tc:
                best = tc[0]
                annotation.gene = best.get("gene_symbol", "")
                annotation.transcript_id = best.get("transcript_id", "")
                annotation.impact = best.get("impact", "")
                annotation.protein_change = best.get("amino_acids", "")

            # Extract allele frequency from colocated variants
            for coloc in top.get("colocated_variants", []):
                freq = coloc.get("frequencies", {})
                if freq:
                    first_allele = list(freq.values())
                    if first_allele and isinstance(first_allele[0], dict):
                        annotation.allele_frequency = list(
                            first_allele[0].values()
                        )[0] or 0.0

        return annotation

    def annotate_vcf_local(self, vcf_path: str, sample_id: str) -> str:
        """
        Annotate a VCF using local VEP installation.

        Parameters
        ----------
        vcf_path : str
            Input VCF path.
        sample_id : str
            Sample identifier for output naming.

        Returns
        -------
        str
            Path to annotated VCF.
        """
        output_vcf = str(self.output_dir / f"{sample_id}.annotated.vcf.gz")
        cmd = [
            "vep",
            "--input_file", vcf_path,
            "--output_file", output_vcf,
            "--format", "vcf",
            "--vcf",
            "--assembly", self.genome_build,
            "--everything",
            "--fork", "4",
            "--compress_output", "bgzip",
            "--force_overwrite",
        ]
        logger.info(f"Running VEP annotation for {sample_id}")
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if result.returncode != 0:
            raise RuntimeError(f"VEP annotation failed:\n{result.stderr}")
        logger.info(f"Annotated VCF written: {output_vcf}")
        return output_vcf

    def annotate_variants_api(
        self, vcf_path: str, max_variants: int = 100
    ) -> List[VariantAnnotation]:
        """
        Annotate variants from a VCF using Ensembl REST API.

        Best suited for small variant sets (< 200 variants).

        Parameters
        ----------
        vcf_path : str
            Path to input VCF file.
        max_variants : int
            Maximum number of variants to annotate (default: 100).

        Returns
        -------
        list of VariantAnnotation
        """
        variants = self._parse_vcf_variants(vcf_path)[:max_variants]
        annotations = []

        for var in variants:
            hgvs = f"{var['chrom']}:g.{var['pos']}{var['ref']}>{var['alt']}"
            raw = self._annotate_via_rest(hgvs)
            annotation = self._parse_vep_annotation(raw, var)
            annotations.append(annotation)

        logger.info(f"Annotated {len(annotations)} variants via REST API")
        return annotations

    def write_annotation_report(
        self, annotations: List[VariantAnnotation], sample_id: str
    ) -> str:
        """
        Write annotations to a structured JSON report.

        Parameters
        ----------
        annotations : list of VariantAnnotation
        sample_id : str

        Returns
        -------
        str
            Path to the JSON report file.
        """
        report_path = str(self.output_dir / f"{sample_id}.annotation_report.json")
        report = {
            "sample_id": sample_id,
            "total_variants_annotated": len(annotations),
            "variants": [a.to_dict() for a in annotations],
        }
        with open(report_path, "w") as f:
            json.dump(report, f, indent=2)
        logger.info(f"Annotation report saved: {report_path}")
        return report_path
