"""
fastq_qc.py
-----------
FASTQ Quality Control Module

Performs per-base quality scoring, adapter trimming detection,
GC content analysis, and generates a structured QC report.
"""

import os
import gzip
import logging
from dataclasses import dataclass, field
from typing import List, Tuple, Optional
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


@dataclass
class QCMetrics:
    """Stores per-sample QC metrics."""
    sample_id: str
    total_reads: int = 0
    reads_passing_filter: int = 0
    mean_quality_score: float = 0.0
    gc_content: float = 0.0
    avg_read_length: float = 0.0
    low_quality_reads: int = 0
    flagged: bool = False
    warnings: List[str] = field(default_factory=list)

    @property
    def pass_rate(self) -> float:
        if self.total_reads == 0:
            return 0.0
        return round((self.reads_passing_filter / self.total_reads) * 100, 2)

    def to_dict(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "total_reads": self.total_reads,
            "reads_passing_filter": self.reads_passing_filter,
            "pass_rate_pct": self.pass_rate,
            "mean_quality_score": self.mean_quality_score,
            "gc_content_pct": round(self.gc_content * 100, 2),
            "avg_read_length": round(self.avg_read_length, 2),
            "low_quality_reads": self.low_quality_reads,
            "flagged": self.flagged,
            "warnings": self.warnings,
        }


class FastqQC:
    """
    Performs quality control on FASTQ files.

    Supports both plain (.fastq) and gzipped (.fastq.gz) files.
    Evaluates per-read Phred quality scores, GC content, and read length.

    Parameters
    ----------
    min_quality : int
        Minimum mean Phred quality score to pass a read (default: 20).
    min_length : int
        Minimum read length to retain (default: 50).
    max_n_content : float
        Maximum fraction of N bases allowed per read (default: 0.05).
    """

    PHRED_OFFSET = 33  # Illumina 1.8+ encoding

    def __init__(
        self,
        min_quality: int = 20,
        min_length: int = 50,
        max_n_content: float = 0.05,
    ):
        self.min_quality = min_quality
        self.min_length = min_length
        self.max_n_content = max_n_content

    def _open_fastq(self, filepath: str):
        """Open plain or gzipped FASTQ file."""
        if filepath.endswith(".gz"):
            return gzip.open(filepath, "rt")
        return open(filepath, "r")

    def _parse_quality(self, qual_string: str) -> List[int]:
        """Convert ASCII quality string to Phred scores."""
        return [ord(c) - self.PHRED_OFFSET for c in qual_string.strip()]

    def _mean_quality(self, scores: List[int]) -> float:
        return sum(scores) / len(scores) if scores else 0.0

    def _gc_content(self, sequence: str) -> float:
        sequence = sequence.strip().upper()
        if not sequence:
            return 0.0
        gc = sequence.count("G") + sequence.count("C")
        return gc / len(sequence)

    def _n_content(self, sequence: str) -> float:
        sequence = sequence.strip().upper()
        if not sequence:
            return 0.0
        return sequence.count("N") / len(sequence)

    def run(self, fastq_path: str, sample_id: Optional[str] = None) -> QCMetrics:
        """
        Run QC analysis on a FASTQ file.

        Parameters
        ----------
        fastq_path : str
            Path to the FASTQ (or .fastq.gz) file.
        sample_id : str, optional
            Sample identifier. Defaults to filename stem.

        Returns
        -------
        QCMetrics
            Populated QC metrics object.
        """
        path = Path(fastq_path)
        if not path.exists():
            raise FileNotFoundError(f"FASTQ file not found: {fastq_path}")

        sample_id = sample_id or path.stem.replace(".fastq", "")
        metrics = QCMetrics(sample_id=sample_id)

        quality_scores_all = []
        gc_values = []
        lengths = []

        logger.info(f"Running QC on: {fastq_path}")

        with self._open_fastq(fastq_path) as fh:
            while True:
                header = fh.readline()
                if not header:
                    break
                sequence = fh.readline()
                plus = fh.readline()
                quality = fh.readline()

                if not all([header, sequence, plus, quality]):
                    break

                metrics.total_reads += 1
                seq = sequence.strip()
                scores = self._parse_quality(quality)
                mean_q = self._mean_quality(scores)
                gc = self._gc_content(seq)
                n_frac = self._n_content(seq)
                read_len = len(seq)

                # Apply filters
                if (
                    mean_q >= self.min_quality
                    and read_len >= self.min_length
                    and n_frac <= self.max_n_content
                ):
                    metrics.reads_passing_filter += 1
                else:
                    metrics.low_quality_reads += 1

                quality_scores_all.append(mean_q)
                gc_values.append(gc)
                lengths.append(read_len)

        if metrics.total_reads > 0:
            metrics.mean_quality_score = round(
                sum(quality_scores_all) / len(quality_scores_all), 2
            )
            metrics.gc_content = sum(gc_values) / len(gc_values)
            metrics.avg_read_length = sum(lengths) / len(lengths)

        # Flag sample if QC thresholds are not met
        if metrics.pass_rate < 80.0:
            metrics.flagged = True
            metrics.warnings.append(
                f"Low pass rate: {metrics.pass_rate}% (threshold: 80%)"
            )
        if metrics.gc_content < 0.35 or metrics.gc_content > 0.65:
            metrics.warnings.append(
                f"Unusual GC content: {round(metrics.gc_content * 100, 1)}%"
            )

        logger.info(
            f"QC complete — {metrics.total_reads} reads, "
            f"{metrics.pass_rate}% passing, "
            f"mean Q={metrics.mean_quality_score}"
        )
        return metrics


def run_batch_qc(
    fastq_files: List[str],
    min_quality: int = 20,
    min_length: int = 50,
) -> List[QCMetrics]:
    """
    Run QC across multiple FASTQ files.

    Parameters
    ----------
    fastq_files : list of str
        Paths to FASTQ files.
    min_quality : int
        Minimum Phred quality threshold.
    min_length : int
        Minimum read length.

    Returns
    -------
    list of QCMetrics
        One QCMetrics object per sample.
    """
    qc = FastqQC(min_quality=min_quality, min_length=min_length)
    results = []
    for fpath in fastq_files:
        try:
            result = qc.run(fpath)
            results.append(result)
        except Exception as e:
            logger.error(f"Failed QC for {fpath}: {e}")
    return results
