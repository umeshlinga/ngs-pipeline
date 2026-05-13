"""
aligner.py
----------
Sequence Alignment Module

Wraps BWA-MEM and STAR aligners with subprocess calls,
handles SAM → BAM conversion, sorting, and indexing via samtools.
"""

import os
import subprocess
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Literal

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

AlignerType = Literal["bwa", "star"]


@dataclass
class AlignmentResult:
    """Stores alignment output metadata."""
    sample_id: str
    bam_path: str
    aligner: str
    total_reads: int = 0
    mapped_reads: int = 0
    mapping_rate: float = 0.0
    success: bool = False
    error_message: str = ""

    @property
    def summary(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "bam_path": self.bam_path,
            "aligner": self.aligner,
            "total_reads": self.total_reads,
            "mapped_reads": self.mapped_reads,
            "mapping_rate_pct": round(self.mapping_rate * 100, 2),
            "success": self.success,
        }


class SequenceAligner:
    """
    Aligns FASTQ reads to a reference genome using BWA-MEM or STAR.

    After alignment, converts SAM → BAM, sorts, and indexes using samtools.

    Parameters
    ----------
    reference : str
        Path to the reference genome FASTA (BWA) or STAR genome directory.
    aligner : str
        Aligner to use: 'bwa' or 'star'.
    threads : int
        Number of CPU threads (default: 4).
    output_dir : str
        Directory to write BAM output files.
    """

    def __init__(
        self,
        reference: str,
        aligner: AlignerType = "bwa",
        threads: int = 4,
        output_dir: str = "data/aligned",
    ):
        self.reference = reference
        self.aligner = aligner
        self.threads = threads
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def _run_command(self, cmd: list, step: str) -> subprocess.CompletedProcess:
        """Execute a shell command and log output."""
        logger.info(f"Running {step}: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result.returncode != 0:
            logger.error(f"{step} failed:\n{result.stderr}")
            raise RuntimeError(f"{step} failed with code {result.returncode}")
        return result

    def _bwa_align(self, fastq_r1: str, fastq_r2: Optional[str], sam_path: str):
        """Run BWA-MEM alignment."""
        cmd = [
            "bwa", "mem",
            "-t", str(self.threads),
            self.reference,
            fastq_r1,
        ]
        if fastq_r2:
            cmd.append(fastq_r2)

        with open(sam_path, "w") as sam_out:
            logger.info(f"BWA-MEM aligning to {self.reference}")
            result = subprocess.run(cmd, stdout=sam_out, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"BWA-MEM failed:\n{result.stderr}")

    def _star_align(self, fastq_r1: str, fastq_r2: Optional[str], sample_id: str):
        """Run STAR alignment for RNA-seq."""
        cmd = [
            "STAR",
            "--runThreadN", str(self.threads),
            "--genomeDir", self.reference,
            "--readFilesIn", fastq_r1,
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--outFileNamePrefix", str(self.output_dir / sample_id) + "/",
            "--outSAMattributes", "NH", "HI", "AS", "NM",
        ]
        if fastq_r2:
            cmd.append(fastq_r2)
        if fastq_r1.endswith(".gz"):
            cmd += ["--readFilesCommand", "zcat"]

        self._run_command(cmd, "STAR alignment")

    def _sam_to_sorted_bam(self, sam_path: str, bam_path: str):
        """Convert SAM → sorted BAM → index using samtools."""
        # SAM → BAM
        self._run_command(
            ["samtools", "view", "-@", str(self.threads), "-bS", sam_path, "-o", bam_path + ".unsorted.bam"],
            "samtools view"
        )
        # Sort BAM
        self._run_command(
            ["samtools", "sort", "-@", str(self.threads), bam_path + ".unsorted.bam", "-o", bam_path],
            "samtools sort"
        )
        # Index BAM
        self._run_command(
            ["samtools", "index", bam_path],
            "samtools index"
        )
        # Cleanup unsorted
        os.remove(bam_path + ".unsorted.bam")
        logger.info(f"Sorted BAM written: {bam_path}")

    def _flagstat(self, bam_path: str) -> dict:
        """Parse samtools flagstat output for mapping metrics."""
        result = subprocess.run(
            ["samtools", "flagstat", bam_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        metrics = {"total": 0, "mapped": 0}
        for line in result.stdout.splitlines():
            if "in total" in line:
                metrics["total"] = int(line.split()[0])
            elif "mapped (" in line:
                metrics["mapped"] = int(line.split()[0])
        return metrics

    def align(
        self,
        fastq_r1: str,
        sample_id: str,
        fastq_r2: Optional[str] = None,
    ) -> AlignmentResult:
        """
        Align reads for a single sample.

        Parameters
        ----------
        fastq_r1 : str
            Path to R1 FASTQ file.
        sample_id : str
            Sample identifier used for output naming.
        fastq_r2 : str, optional
            Path to R2 FASTQ file (paired-end).

        Returns
        -------
        AlignmentResult
            Alignment output metadata and QC metrics.
        """
        result = AlignmentResult(
            sample_id=sample_id,
            bam_path=str(self.output_dir / f"{sample_id}.sorted.bam"),
            aligner=self.aligner,
        )

        try:
            if self.aligner == "bwa":
                sam_path = str(self.output_dir / f"{sample_id}.sam")
                self._bwa_align(fastq_r1, fastq_r2, sam_path)
                self._sam_to_sorted_bam(sam_path, result.bam_path)
                os.remove(sam_path)
            elif self.aligner == "star":
                self._star_align(fastq_r1, fastq_r2, sample_id)
                result.bam_path = str(
                    self.output_dir / sample_id / "Aligned.sortedByCoord.out.bam"
                )

            # Collect mapping statistics
            stats = self._flagstat(result.bam_path)
            result.total_reads = stats["total"]
            result.mapped_reads = stats["mapped"]
            if result.total_reads > 0:
                result.mapping_rate = result.mapped_reads / result.total_reads
            result.success = True

            logger.info(
                f"Alignment complete — {result.mapped_reads}/{result.total_reads} "
                f"reads mapped ({result.mapping_rate*100:.1f}%)"
            )

        except Exception as e:
            result.error_message = str(e)
            logger.error(f"Alignment failed for {sample_id}: {e}")

        return result
