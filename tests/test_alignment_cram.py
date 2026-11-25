"""Tests for BAM/CRAM helpers (scripts/lib/alignment.py) and scripts/bam_to_cram.sh, on synthetic reads only."""

from __future__ import annotations

import hashlib
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

SCRIPT_DIR = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.alignment import (
    caller_cram_problem,
    cram_version,
    index_stats,
    is_cram,
    open_alignment,
    parse_idxstats,
    t1k_regions,
)

HAVE_SAMTOOLS = shutil.which("samtools") is not None
needs_samtools = pytest.mark.skipif(not HAVE_SAMTOOLS, reason="samtools not installed")

# Synthetic genome: contig "1" long enough for the converter's 1:1000000-1100000 check, contig "X", and a decoy
CONTIGS = [("1", 1_200_000), ("X", 300_000), ("Un_SYN0001v1_decoy", 5_000)]


def _write_bam(path: Path, n_pairs: int = 400, seed: int = 7) -> Path:
    rng = random.Random(seed)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": n, "LN": ln} for n, ln in CONTIGS],
              "RG": [{"ID": "rg1", "SM": "SYNTHETIC"}]}
    reads = []
    for i in range(n_pairs):
        tid = 0 if i % 4 else (1 if i % 8 else 2)
        span = CONTIGS[tid][1]
        lo = 1_000_000 if tid == 0 and i % 2 else 0
        pos = rng.randrange(lo, min(span - 400, lo + 150_000))
        for mate in (0, 1):
            a = pysam.AlignedSegment()
            a.query_name = f"pair{i}"
            a.query_sequence = "".join(rng.choice("ACGT") for _ in range(100))
            a.flag = 1 | 2 | (64 if mate == 0 else 128) | (16 if mate else 32)
            a.reference_id = tid
            a.reference_start = pos + mate * 200
            a.mapping_quality = 60
            a.cigar = [(0, 100)]
            a.next_reference_id = tid
            a.next_reference_start = pos + (1 - mate) * 200
            a.template_length = 300 if mate == 0 else -300
            a.query_qualities = pysam.qualitystring_to_array("".join(rng.choice("#+5?I") for _ in range(100)))
            a.set_tag("RG", "rg1")
            reads.append(a)
    unmapped = pysam.AlignedSegment()
    unmapped.query_name = "lonely"
    unmapped.query_sequence = "ACGT" * 25
    unmapped.flag = 4
    unmapped.query_qualities = pysam.qualitystring_to_array("I" * 100)
    reads.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for r in reads + [unmapped]:
            out.write(r)
    pysam.index(str(path))
    return path


def _md5(p: Path) -> str:
    return hashlib.md5(p.read_bytes()).hexdigest()


def test_cram_version_and_caller_check(tmp_path):
    c30, c31, fake, bam = tmp_path / "a.cram", tmp_path / "b.cram", tmp_path / "c.cram", tmp_path / "d.bam"
    c30.write_bytes(b"CRAM\x03\x00" + b"\0" * 20)
    c31.write_bytes(b"CRAM\x03\x01" + b"\0" * 20)
    fake.write_bytes(b"BAM\x01")
    bam.write_bytes(b"\x1f\x8b")
    assert cram_version(c30) == (3, 0) and cram_version(c31) == (3, 1)
    assert cram_version(fake) is None and cram_version(tmp_path / "missing.cram") is None
    assert caller_cram_problem(c30) is None
    assert "CRAM 3.1" in caller_cram_problem(c31) and "bam_to_cram.sh" in caller_cram_problem(c31)
    assert "does not start with a CRAM header" in caller_cram_problem(fake)
    assert caller_cram_problem(bam) is None
    assert is_cram(c30) and not is_cram(bam)


def test_parse_idxstats_drops_unplaced_row():
    text = "1\t1200000\t300\t2\nX\t300000\t50\t0\n*\t0\t0\t1\n"
    assert parse_idxstats(text) == [("1", 1200000, 300), ("X", 300000, 50)]


def test_t1k_regions_merges_genes_and_adds_non_main_contigs(tmp_path):
    coord = tmp_path / "coord.fa"
    coord.write_text(
        ">GENEA*01 6 10000 12000 +\nACGT\n>GENEA*02 6 11500 13000 +\nACGT\n>GENEB*01 6 50000 51000 +\nACGT\n"
        ">GENEC*01 19 2000 3000 -\nACGT\n>GENED*01 7 100 200 +\nACGT\n", encoding="utf-8")
    contigs = [("6", 60_000), ("19", 2_500), ("1", 1000), ("chrX", 500), ("MT", 16_569),
               ("Un_SYN1v1_decoy", 900), ("EBV", 171_823), ("1_SYN2v1_random", 700)]
    got = t1k_regions(coord, contigs, pad=1000)
    assert got[:3] == [("19", 1000, 2500), ("6", 9000, 14000), ("6", 49000, 52000)]  # gene on '7' (no contig) dropped
    assert got[3:] == [("Un_SYN1v1_decoy", 0, 900), ("EBV", 0, 171_823), ("1_SYN2v1_random", 0, 700)]


def test_index_stats_and_open_alignment_on_bam(tmp_path):
    bam = _write_bam(tmp_path / "s.bam")
    stats = dict((c, m) for c, _, m in index_stats(bam))
    assert stats["1"] > 0 and stats["X"] > 0 and stats["Un_SYN0001v1_decoy"] > 0
    with open_alignment(bam, tmp_path / "no_such.fa") as af:
        assert af.references == tuple(n for n, _ in CONTIGS)


@needs_samtools
def test_index_stats_reads_cram_counts_that_pysam_reports_as_zero(tmp_path):
    bam = _write_bam(tmp_path / "s.bam")
    cram = tmp_path / "s.cram"
    subprocess.run(["samtools", "view", "-C", "--output-fmt-option", "no_ref=1", "--output-fmt-option", "version=3.0",
                    "-o", str(cram), str(bam)], check=True)
    subprocess.run(["samtools", "index", str(cram)], check=True)
    assert cram_version(cram) == (3, 0)
    with pysam.AlignmentFile(str(cram)) as af:
        assert all(s.mapped == 0 for s in af.get_index_statistics())  # the trap index_stats avoids
    assert index_stats(cram) == index_stats(bam)
    with open_alignment(cram) as af:
        assert af.count("1", 1_000_000, 1_100_000) == pysam.AlignmentFile(str(bam)).count("1", 1_000_000, 1_100_000)


CONVERTER = SCRIPT_DIR / "bam_to_cram.sh"


@needs_samtools
def test_converter_writes_verified_cram_and_leaves_bam_untouched(tmp_path):
    (tmp_path / "in").mkdir()
    bam = _write_bam(tmp_path / "in" / "sample.bam")
    before = (_md5(bam), bam.stat().st_mtime_ns, _md5(Path(str(bam) + ".bai")))
    out = tmp_path / "DNA backup" / "out"  # the real BAMs live under a folder name with a space
    proc = subprocess.run(["bash", str(CONVERTER), "--out-dir", str(out), "--threads", "2", str(bam)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr + proc.stdout
    cram = out / "sample.cram"
    assert cram.exists() and (out / "sample.cram.crai").exists()
    assert not list(out.glob("*.part*")) and not list(out.glob(".sample.verify.*"))
    assert cram_version(cram) == (3, 0)
    report = (out / "sample.cram.verify.txt").read_text(encoding="utf-8")
    for line in ("quickcheck:  PASS", "flagstat:    identical", "idxstats:    identical", "checksum -a: identical"):
        assert line in report
    assert "Un_SYN0001v1_decoy" in report  # the non-main contig was region-checked too
    md5_line = (out / "sample.cram.md5").read_text(encoding="utf-8").splitlines()[0]
    assert md5_line.split()[0] == _md5(cram)
    assert (_md5(bam), bam.stat().st_mtime_ns, _md5(Path(str(bam) + ".bai"))) == before
    # a second run refuses to overwrite without --force
    again = subprocess.run(["bash", str(CONVERTER), "--out-dir", str(out), str(bam)], capture_output=True, text=True)
    assert again.returncode != 0 and "exists" in again.stderr


@needs_samtools
def test_converter_refuses_reference_missing_a_contig(tmp_path):
    (tmp_path / "in").mkdir()
    bam = _write_bam(tmp_path / "in" / "sample.bam")
    fa = tmp_path / "partial.fa"
    rng = random.Random(1)
    with open(fa, "w", encoding="utf-8", newline="\n") as fh:  # no decoy contig
        for name, ln in CONTIGS[:2]:
            fh.write(f">{name}\n")
            seq = "".join(rng.choice("ACGT") for _ in range(ln))
            fh.writelines(seq[i:i + 80] + "\n" for i in range(0, ln, 80))
    pysam.faidx(str(fa))
    out = tmp_path / "out"
    proc = subprocess.run(["bash", str(CONVERTER), "--out-dir", str(out), "--reference", str(fa), str(bam)],
                          capture_output=True, text=True)
    assert proc.returncode != 0
    assert "missing from" in proc.stderr and "embed_ref=2" in proc.stderr
    assert not out.exists() or not any(out.iterdir())
