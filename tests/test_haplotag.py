"""Regression tests for --haplotag phasing.

Builds a tiny synthetic reference (one CAG repeat locus) and a BAM with HP/PS-tagged reads,
runs the pipeline into tmp_path and checks the freshly written VCF record.
"""
import random
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

CORE = Path(__file__).resolve().parent.parent / "ATARVA" / "core.py"

FLANK = 1000
MOTIF = "CAG"
REF_UNITS = 10                      # reference allele: (CAG)10 = 30 bp
EXP_UNITS = 14                      # expanded allele:  (CAG)14 = 42 bp
REPEAT_START = FLANK
REPEAT_END = FLANK + REF_UNITS * len(MOTIF)

rng = random.Random(7)
# flanks without CAG-like repeats: draw from A/T/C only
LEFT = "".join(rng.choice("ATC") for _ in range(FLANK))
RIGHT = "".join(rng.choice("ATC") for _ in range(FLANK))
REF_SEQ = LEFT + MOTIF * REF_UNITS + RIGHT


def make_read(name, start, expanded, hp=None, ps=None, qual=30):
    """Read from `start` to the end of the contig; expanded reads carry an insertion at the repeat end."""
    ins = MOTIF * (EXP_UNITS - REF_UNITS) if expanded else ""
    seq = REF_SEQ[start:REPEAT_END] + ins + REF_SEQ[REPEAT_END:]
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = seq
    a.flag = 0
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = 60
    left_m = REPEAT_END - start
    right_m = len(REF_SEQ) - REPEAT_END
    if ins:
        a.cigartuples = [(0, left_m), (1, len(ins)), (0, right_m)]
    else:
        a.cigartuples = [(0, left_m + right_m)]
    a.query_qualities = pysam.qualitystring_to_array(chr(qual + 33) * len(seq))
    tags = [("NM", len(ins)), ("MD", str(left_m + right_m))]
    if hp is not None:
        tags.append(("HP", hp))
    if ps is not None:
        tags.append(("PS", ps))
    a.set_tags(tags)
    return a


def build_inputs(tmp_path, reads):
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chrT\n{REF_SEQ}\n")
    pysam.faidx(str(fa))

    bed = tmp_path / "loci.bed"
    bed.write_text(f"chrT\t{REPEAT_START}\t{REPEAT_END}\t{MOTIF}\t{len(MOTIF)}\n")
    bed_gz = pysam.tabix_index(str(bed), preset="bed", force=True)

    unsorted = tmp_path / "unsorted.bam"
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chrT", "LN": len(REF_SEQ)}]}
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for r in reads:
            out.write(r)
    bam = tmp_path / "reads.bam"
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    return fa, Path(bed_gz), bam


def run_atarva(tmp_path, reads, extra_args):
    fa, bed, bam = build_inputs(tmp_path, reads)
    vcf = tmp_path / "out.vcf"
    result = subprocess.run(
        [sys.executable, str(CORE), "genotype", "-f", str(fa), "-b", str(bam), "-r", str(bed),
         "-o", str(vcf), "-t", "1", *extra_args],
        capture_output=True, text=True,
    )
    # ATaRVa can exit 0 after an error, so check the output as well as the exit code
    assert result.returncode == 0, f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    assert vcf.exists(), f"no VCF written\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    records = [l.rstrip("\n").split("\t") for l in vcf.read_text().splitlines() if not l.startswith("#")]
    assert len(records) == 1, f"expected 1 record, got {len(records)}\nSTDERR:\n{result.stderr}"
    rec = records[0]
    return dict(zip(rec[8].split(":"), rec[9].split(":")))


def phased_reads(n_hp1, n_hp2, ps=500, start=0):
    """HP1 reads carry the expanded allele, HP2 reads the reference allele."""
    reads = []
    for i in range(n_hp1):
        reads.append(make_read(f"hp1_{i}", start + 2 * i, expanded=True, hp=1, ps=ps))
    for i in range(n_hp2):
        reads.append(make_read(f"hp2_{i}", start + 2 * i + 1, expanded=False, hp=2, ps=ps))
    return reads


def test_haplotag_phases_locus(tmp_path):
    """All reads tagged, one phase set: the call comes from the haplotags, HP1 allele first."""
    s = run_atarva(tmp_path, phased_reads(10, 10), ["--haplotag", "HP"])
    assert s["PS"] == "500"
    assert s["GT"] == "1|0"
    assert s["AL"] == "42,30"
    assert s["SD"] == "10,10"


def test_haplotag_falls_back_when_too_many_untagged(tmp_path):
    """>15% untagged reads: the locus is not phased by haplotag (no PS) but is still genotyped."""
    reads = phased_reads(8, 8)
    reads += [make_read(f"untagged_{i}", 40 + i, expanded=(i % 2 == 0)) for i in range(4)]  # 4/20 = 20%
    s = run_atarva(tmp_path, reads, ["--haplotag", "HP"])
    assert s["PS"] == "."
    assert sorted(int(x) for x in s["AL"].split(",")) == [30, 42]


def test_haplotag_ignores_phase_set_of_evicted_read(tmp_path):
    """At --max-reads, the evicted low-quality read's PS must not count (it would make the locus multi-PS)."""
    low = make_read("lowq_ps999", 0, expanded=True, hp=1, ps=999, qual=10)
    reads = [low] + phased_reads(10, 10, start=1)       # 21 reads, max 20: the 21st evicts the low-quality read
    s = run_atarva(tmp_path, reads, ["--haplotag", "HP", "--max-reads", "20"])
    assert s["PS"] == "500"
    assert s["SD"] == "10,10"


def test_no_haplotag_flag_ignores_tags(tmp_path):
    """Without --haplotag the tags are not used: no PS."""
    s = run_atarva(tmp_path, phased_reads(10, 10), [])
    assert s["PS"] == "."
    assert sorted(int(x) for x in s["AL"].split(",")) == [30, 42]
