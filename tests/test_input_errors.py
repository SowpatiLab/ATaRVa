"""Fatal input errors must exit non-zero and report on stderr; help screens still exit 0."""
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

CORE = Path(__file__).resolve().parent.parent / "ATARVA" / "core.py"
DATA = Path(__file__).resolve().parent / "test_data"
FASTA = DATA / "test_ref.fa"
BAM = DATA / "test_reads.bam"
BED = DATA / "test_region.bed.gz"


def run(args, cwd):
    return subprocess.run([sys.executable, str(CORE), *args], cwd=cwd, capture_output=True, text=True)


def genotype(tmp_path, fasta=FASTA, bam=BAM, bed=BED):
    return run(["genotype", "-f", str(fasta), "-b", str(bam), "-r", str(bed),
                "--vcf", str(tmp_path / "out.vcf")], tmp_path)


def write_bam(path, sort_order):
    hd = {"VN": "1.6"}
    if sort_order:
        hd["SO"] = sort_order
    header = {"HD": hd, "SQ": [{"SN": "chr1", "LN": 5000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header):
        pass
    return path


def assert_fatal(result, tmp_path, message):
    assert result.returncode != 0, f"exit code 0 for fatal error; stdout: {result.stdout}"
    assert message in result.stderr, f"{message!r} not on stderr: {result.stderr}"
    assert not (tmp_path / "out.vcf").exists()


def test_missing_fasta(tmp_path):
    assert_fatal(genotype(tmp_path, fasta=tmp_path / "nofile.fa"), tmp_path, "not a valid FASTA file")


def test_missing_bam(tmp_path):
    assert_fatal(genotype(tmp_path, bam=tmp_path / "nofile.bam"), tmp_path, "not a valid alignment file")


def test_missing_regions(tmp_path):
    assert_fatal(genotype(tmp_path, bed=tmp_path / "nofile.bed.gz"), tmp_path, "not a valid tabix file")


@pytest.mark.parametrize("sort_order, message", [
    ("unsorted", "It should be sorted by 'coordinate'"),
    (None, "No sort order specified in the header"),
])
def test_unsorted_bam(tmp_path, sort_order, message):
    bam = write_bam(tmp_path / "unsorted.bam", sort_order)
    assert_fatal(genotype(tmp_path, bam=bam), tmp_path, message)


@pytest.mark.parametrize("args", [[], ["genotype"], ["merge"]])
def test_help_exits_zero(tmp_path, args):
    result = run(args, tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip()
