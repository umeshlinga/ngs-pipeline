"""
test_pipeline.py
----------------
Unit tests for the NGS pipeline modules.
"""

import os
import gzip
import json
import pytest
import tempfile
from unittest.mock import patch, MagicMock

from pipeline.qc.fastq_qc import FastqQC, QCMetrics, run_batch_qc
from pipeline.annotation.annotator import VCFAnnotator, VariantAnnotation


# ── Fixtures ──────────────────────────────────────────────────────

def make_fastq(path: str, num_reads: int = 10, quality_char: str = "I", read_len: int = 100):
    """Write a synthetic FASTQ file."""
    seq = "ATCGATCGATCG" * (read_len // 12) + "ATCG"
    seq = seq[:read_len]
    qual = quality_char * read_len
    with open(path, "w") as f:
        for i in range(num_reads):
            f.write(f"@read_{i}\n{seq}\n+\n{qual}\n")


def make_fastq_gz(path: str, num_reads: int = 10, quality_char: str = "I"):
    """Write a synthetic gzipped FASTQ file."""
    seq = "ATCGATCGATCGATCG" * 6
    qual = quality_char * len(seq)
    with gzip.open(path, "wt") as f:
        for i in range(num_reads):
            f.write(f"@read_{i}\n{seq}\n+\n{qual}\n")


def make_vcf(path: str):
    """Write a minimal VCF file."""
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("chr1\t100\t.\tA\tT\t50\tPASS\t.\n")
        f.write("chr1\t200\t.\tGG\tG\t40\tPASS\t.\n")


# ── QC Tests ──────────────────────────────────────────────────────

class TestFastqQC:

    def test_basic_qc_pass(self, tmp_path):
        """High quality reads should pass QC."""
        fq = str(tmp_path / "sample.fastq")
        make_fastq(fq, num_reads=20, quality_char="I", read_len=100)
        qc = FastqQC(min_quality=20, min_length=50)
        metrics = qc.run(fq, sample_id="test_sample")

        assert metrics.total_reads == 20
        assert metrics.reads_passing_filter == 20
        assert metrics.pass_rate == 100.0
        assert not metrics.flagged
        assert metrics.mean_quality_score > 20

    def test_low_quality_reads_flagged(self, tmp_path):
        """Low quality reads should be filtered and sample flagged."""
        fq = str(tmp_path / "low_qual.fastq")
        make_fastq(fq, num_reads=20, quality_char="!", read_len=100)  # Phred=0
        qc = FastqQC(min_quality=20, min_length=50)
        metrics = qc.run(fq, sample_id="low_qual")

        assert metrics.low_quality_reads == 20
        assert metrics.reads_passing_filter == 0
        assert metrics.flagged
        assert len(metrics.warnings) > 0

    def test_short_reads_filtered(self, tmp_path):
        """Reads below min_length should not pass."""
        fq = str(tmp_path / "short.fastq")
        make_fastq(fq, num_reads=10, quality_char="I", read_len=30)
        qc = FastqQC(min_quality=20, min_length=50)
        metrics = qc.run(fq, sample_id="short")

        assert metrics.reads_passing_filter == 0

    def test_gzipped_fastq(self, tmp_path):
        """Gzipped FASTQ files should be processed correctly."""
        fq_gz = str(tmp_path / "sample.fastq.gz")
        make_fastq_gz(fq_gz, num_reads=15, quality_char="I")
        qc = FastqQC(min_quality=20, min_length=50)
        metrics = qc.run(fq_gz, sample_id="gz_sample")

        assert metrics.total_reads == 15
        assert not metrics.flagged

    def test_missing_file_raises(self):
        """Non-existent file should raise FileNotFoundError."""
        qc = FastqQC()
        with pytest.raises(FileNotFoundError):
            qc.run("/nonexistent/path/sample.fastq")

    def test_qc_metrics_to_dict(self):
        """QCMetrics.to_dict() should return all expected keys."""
        m = QCMetrics(sample_id="s1", total_reads=100, reads_passing_filter=95)
        d = m.to_dict()
        assert "sample_id" in d
        assert "pass_rate_pct" in d
        assert d["pass_rate_pct"] == 95.0

    def test_batch_qc(self, tmp_path):
        """Batch QC should return one result per file."""
        paths = []
        for i in range(3):
            fq = str(tmp_path / f"sample_{i}.fastq")
            make_fastq(fq, num_reads=5)
            paths.append(fq)

        results = run_batch_qc(paths, min_quality=20, min_length=50)
        assert len(results) == 3
        assert all(isinstance(r, QCMetrics) for r in results)


# ── Annotator Tests ───────────────────────────────────────────────

class TestVCFAnnotator:

    def test_parse_vcf_variants(self, tmp_path):
        """VCF parser should extract correct number of variants."""
        vcf = str(tmp_path / "test.vcf")
        make_vcf(vcf)
        annotator = VCFAnnotator(output_dir=str(tmp_path))
        variants = annotator._parse_vcf_variants(vcf)
        assert len(variants) == 2
        assert variants[0]["chrom"] == "chr1"
        assert variants[0]["pos"] == 100
        assert variants[0]["ref"] == "A"
        assert variants[0]["alt"] == "T"

    def test_write_annotation_report(self, tmp_path):
        """Annotation report should be valid JSON with expected structure."""
        annotator = VCFAnnotator(output_dir=str(tmp_path))
        annotations = [
            VariantAnnotation(
                chrom="chr1", pos=100, ref="A", alt="T",
                gene="BRCA1", consequence="missense_variant",
            )
        ]
        report_path = annotator.write_annotation_report(annotations, "sample1")
        assert os.path.exists(report_path)

        with open(report_path) as f:
            report = json.load(f)

        assert report["sample_id"] == "sample1"
        assert report["total_variants_annotated"] == 1
        assert report["variants"][0]["gene"] == "BRCA1"

    def test_annotation_to_dict(self):
        """VariantAnnotation.to_dict() should include all fields."""
        ann = VariantAnnotation(
            chrom="chr2", pos=500, ref="G", alt="A",
            gene="TP53", consequence="stop_gained", impact="HIGH",
        )
        d = ann.to_dict()
        assert d["gene"] == "TP53"
        assert d["consequence"] == "stop_gained"
        assert d["impact"] == "HIGH"

    @patch("pipeline.annotation.annotator.requests.get")
    def test_rest_api_annotation(self, mock_get, tmp_path):
        """REST API annotation should parse VEP response correctly."""
        mock_response = MagicMock()
        mock_response.json.return_value = [{
            "most_severe_consequence": "missense_variant",
            "transcript_consequences": [{
                "gene_symbol": "EGFR",
                "transcript_id": "ENST00000275493",
                "impact": "MODERATE",
                "amino_acids": "V/L",
            }],
            "colocated_variants": [],
        }]
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        annotator = VCFAnnotator(output_dir=str(tmp_path), use_rest_api=True)
        result = annotator._annotate_via_rest("chr7:g.55174772A>T")
        parsed = annotator._parse_vep_annotation(result, {
            "chrom": "chr7", "pos": 55174772, "ref": "A", "alt": "T"
        })

        assert parsed.gene == "EGFR"
        assert parsed.consequence == "missense_variant"
        assert parsed.impact == "MODERATE"
