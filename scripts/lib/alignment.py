"""BAM / CRAM helpers shared by Step 3 and Step 10.

CRAM differences that matter here:
  - Decoding a reference-based CRAM needs its FASTA; a no_ref CRAM (what scripts/bam_to_cram.sh writes by default)
    needs none. pysam and the callers are given the reference whenever the input is a CRAM.
  - A CRAM index (.crai) has no per-contig read counts: pysam's get_index_statistics() returns zeros for a CRAM
    instead of failing, so the counts come from `samtools idxstats`, which reads them from the slice headers.
  - ExpansionHunter v5 bundles an htslib that cannot decode CRAM 3.1 codecs, so CRAMs must be version 3.0.
  - T1K's bam-extractor (samtools 0.1.19 API) cannot read CRAM at all and opens its input twice, so it cannot read a
    pipe either; Step 10 hands it a temporary BAM instead (see t1k_input_bam).
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any

# The newest CRAM version every Step 10 caller can decode (ExpansionHunter v5: 3.0)
MAX_CALLER_CRAM_VERSION = (3, 0)


def is_cram(path: Path) -> bool:
    return Path(path).suffix.lower() == ".cram"


def cram_version(path: Path) -> tuple[int, int] | None:
    """(major, minor) from a CRAM file definition (bytes 'CRAM', major, minor), or None if not a CRAM."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(6)
    except OSError:
        return None
    if len(head) < 6 or head[:4] != b"CRAM":
        return None
    return head[4], head[5]


def open_alignment(path: Path, reference: Path | None = None):
    """pysam.AlignmentFile for a BAM or CRAM; a CRAM gets the reference FASTA when one is given and exists."""
    import pysam

    kwargs: dict[str, Any] = {}
    if is_cram(path) and reference is not None and Path(reference).exists():
        kwargs["reference_filename"] = str(reference)
    return pysam.AlignmentFile(str(path), "rc" if is_cram(path) else "rb", **kwargs)


def parse_idxstats(text: str) -> list[tuple[str, int, int]]:
    """(contig, length, mapped) rows from `samtools idxstats` output, without the '*' unplaced row."""
    rows = []
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) >= 3 and f[0] != "*":
            rows.append((f[0], int(f[1]), int(f[2])))
    return rows


def index_stats(path: Path, samtools: str | None = None) -> list[tuple[str, int, int]]:
    """Per-contig (contig, length, mapped reads) for a BAM (from the .bai) or a CRAM (via samtools idxstats)."""
    if is_cram(path):
        exe = samtools or shutil.which("samtools") or "samtools"
        proc = subprocess.run([exe, "idxstats", str(path)], capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"samtools idxstats failed on {Path(path).name}: {proc.stderr.strip()[-400:]}")
        return parse_idxstats(proc.stdout)
    with open_alignment(path) as af:
        return [(s.contig, af.get_reference_length(s.contig), s.mapped) for s in af.get_index_statistics()]


def caller_cram_problem(path: Path) -> str | None:
    """Why a CRAM cannot go to ExpansionHunter (too new a CRAM version), or None when it can."""
    if not is_cram(path):
        return None
    ver = cram_version(path)
    if ver is None:
        return f"{Path(path).name} does not start with a CRAM header"
    if ver > MAX_CALLER_CRAM_VERSION:
        return (f"{Path(path).name} is CRAM {ver[0]}.{ver[1]}, which ExpansionHunter v5 cannot decode; "
                "re-encode it as CRAM 3.0 with scripts/bam_to_cram.sh")
    return None


def t1k_regions(coord_fa: Path, header_contigs: list[tuple[str, int]], pad: int = 1000) -> list[tuple[str, int, int]]:
    """Regions whose reads reproduce T1K's own whole-file extraction from a BAM.

    T1K picks candidate reads by k-mer from the whole file, not only from its gene coordinates; most of those that lie
    outside the HLA genes are on the decoy / unplaced / EBV contigs. So: the coordinate-file genes (+/- pad, merged)
    plus every contig that is not a main chromosome. Unmapped reads are added separately.
    """
    main = {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}
    lengths = dict(header_contigs)
    genes: dict[str, list[tuple[int, int]]] = {}
    with open(coord_fa, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            f = line[1:].split()
            if len(f) >= 4 and f[2].isdigit() and f[3].isdigit():
                genes.setdefault(f[1], []).append((max(0, int(f[2]) - pad), int(f[3]) + pad))
    out: list[tuple[str, int, int]] = []
    for chrom in sorted(genes):
        if chrom not in lengths:
            continue
        cur = None
        for s, e in sorted(genes[chrom]):
            if cur and s <= cur[1]:
                cur = (cur[0], max(cur[1], e))
            else:
                if cur:
                    out.append((chrom, *cur))
                cur = (s, e)
        if cur:
            out.append((chrom, cur[0], min(cur[1], lengths[chrom])))
    out += [(c, 0, n) for c, n in header_contigs if c.removeprefix("chr") not in main and c not in genes]
    return out


def t1k_input_bam(cram: Path, coord_fa: Path, out_dir: Path, reference: Path | None = None, threads: int = 4,
                  samtools: str | None = None) -> Path:
    """Write the reads T1K needs from a CRAM to out_dir/t1k_input.bam (indexed) and return its path."""
    exe = samtools or shutil.which("samtools") or "samtools"
    out_dir.mkdir(parents=True, exist_ok=True)
    with open_alignment(cram, reference) as af:
        contigs = list(zip(af.references, af.lengths))
    bed = out_dir / "t1k_input_regions.bed"
    with open(bed, "w", encoding="utf-8", newline="\n") as fh:
        for c, s, e in t1k_regions(coord_fa, contigs):
            fh.write(f"{c}\t{s}\t{e}\n")
    ref_args = ["-T", str(reference)] if reference is not None and Path(reference).exists() else []
    mapped, unmapped, final = out_dir / "t1k_mapped.bam", out_dir / "t1k_unmapped.bam", out_dir / "t1k_input.bam"
    steps = [
        [exe, "view", "-@", str(threads), *ref_args, "-b", "-M", "-L", str(bed), "-o", str(mapped), str(cram)],
        [exe, "view", "-@", str(threads), *ref_args, "-b", "-o", str(unmapped), str(cram), "*"],
        [exe, "cat", "-o", str(final), str(mapped), str(unmapped)],
        [exe, "index", str(final)],
    ]
    for cmd in steps:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"{' '.join(cmd[:2])} failed while preparing T1K input: {proc.stderr.strip()[-400:]}")
    mapped.unlink(missing_ok=True)
    unmapped.unlink(missing_ok=True)
    return final
