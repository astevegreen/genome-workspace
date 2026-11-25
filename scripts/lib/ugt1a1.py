"""UGT1A1 promoter TA repeat (rs3064744, *28) counted from BAM reads.

The GRCh38 reference carries A(TA)6TAA at chr2:233760233 (UGT1A1*1). Short-read VCFs often show a
reference block across the repeat, which says nothing about the repeat length. Here each read pair that
contains both flanking anchors is measured directly: the sequence between the anchors must be a pure
(AT)n run, and its length gives the TA count. Read pairs, not reads, are counted, so overlapping mates
are one observation.

The diplotype is the most likely pair of repeat lengths under a simple stutter model, and it is called
only when there are enough spanning fragments and it beats the next-best pair by MIN_LOD (log10).
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

UGT1A1_CHROM = "2"
UGT1A1_TA_POS = 233_760_233  # 1-based position of the C before the repeat (dbSNP rs3064744 REF 'CAT')
# Sequence either side of the A(TA)nTAA repeat in GRCh38: ...GGTTTTTGCC [AT]x(n+1) AAGTAGGAG...
LEFT_ANCHOR = "GGTTTTTGCC"
RIGHT_ANCHOR = "AAGTAGGAG"
FETCH_PAD = 20

MIN_MAPQ = 20
MIN_SPANNING = 10  # spanning read pairs needed before a call is made
MIN_LOD = 2.0  # best diplotype must be >= 100x more likely than the next best
STUTTER_ERROR = 0.02  # chance a fragment shows a length other than its allele's

TA_STAR = {5: "*36", 6: "*1", 7: "*28", 8: "*37"}
DECREASED_FUNCTION = {"*6", "*28", "*37"}  # CPIC; *1 and *36 are normal function
_AT_RUN = re.compile(r"(?:AT)+")


def ta_units(read_seq: str) -> int | None | str:
    """TA count from one read: an int, None if the read lacks either anchor, 'noncanonical' if the gap is not (AT)n."""
    i = read_seq.find(LEFT_ANCHOR)
    if i < 0:
        return None
    start = i + len(LEFT_ANCHOR)
    j = read_seq.find(RIGHT_ANCHOR, start)
    if j < 0:
        return None
    middle = read_seq[start:j]
    if not _AT_RUN.fullmatch(middle):
        return "noncanonical"
    return len(middle) // 2 - 1  # A(TA)nTAA reads as (AT)x(n+1) between the anchors


def count_ta_repeats(
    bam_path: Path, has_chr: bool = False, min_mapq: int = MIN_MAPQ, reference: Path | None = None
) -> dict[str, Any]:
    """Count TA repeat lengths over read pairs that fully span the repeat with both anchors (BAM or CRAM)."""
    from lib.alignment import open_alignment

    chrom = f"chr{UGT1A1_CHROM}" if has_chr else UGT1A1_CHROM
    per_fragment: dict[str, int | str] = {}
    partial = 0
    with open_alignment(bam_path, reference) as af:
        for r in af.fetch(chrom, UGT1A1_TA_POS - FETCH_PAD, UGT1A1_TA_POS + FETCH_PAD):
            if r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary or r.is_qcfail:
                continue
            if r.mapping_quality < min_mapq or not r.query_sequence:
                continue
            n = ta_units(r.query_sequence)
            if n is None:
                partial += 1
                continue
            prev = per_fragment.get(r.query_name)
            if prev is None:
                per_fragment[r.query_name] = n
            elif prev != n:
                per_fragment[r.query_name] = "mate_conflict"
    return summarize_fragments(list(per_fragment.values()), partial)


def summarize_fragments(values: list[int | str], partial_reads: int = 0) -> dict[str, Any]:
    counts: dict[int, int] = {}
    for v in values:
        if isinstance(v, int):
            counts[v] = counts.get(v, 0) + 1
    return {
        "counts": dict(sorted(counts.items())),
        "spanning_fragments": sum(counts.values()),
        "noncanonical_fragments": sum(1 for v in values if v == "noncanonical"),
        "mate_conflict_fragments": sum(1 for v in values if v == "mate_conflict"),
        "partial_reads": partial_reads,
    }


def _p_obs(observed: int, allele: int, error: float) -> float:
    if observed == allele:
        return 1.0 - error
    if abs(observed - allele) == 1:
        return error / 2
    return error / 100


def call_ta_diplotype(
    counts: dict[int, int],
    min_spanning: int = MIN_SPANNING,
    min_lod: float = MIN_LOD,
    error: float = STUTTER_ERROR,
) -> dict[str, Any]:
    """Most likely pair of TA lengths. status is OK, LOW_DEPTH or UNCERTAIN."""
    counts = {int(k): int(v) for k, v in counts.items() if int(v) > 0}
    total = sum(counts.values())
    out: dict[str, Any] = {"spanning_fragments": total, "min_spanning": min_spanning, "min_lod": min_lod}
    if total < min_spanning:
        return {**out, "status": "LOW_DEPTH", "alleles": None, "lod": None, "candidates": []}
    # neighbours too, so an all-TA6 pileup is weighed against TA5/TA6 and TA6/TA7 hets
    lengths = sorted({m for k in counts for m in (k - 1, k, k + 1) if m > 0})
    scored: list[tuple[float, tuple[int, int]]] = []
    for i, a in enumerate(lengths):
        for b in lengths[i:]:
            ll = sum(n * math.log(0.5 * _p_obs(k, a, error) + 0.5 * _p_obs(k, b, error)) for k, n in counts.items())
            scored.append((ll, (a, b)))
    scored.sort(key=lambda s: -s[0])
    best_ll, best = scored[0]
    lod = (best_ll - scored[1][0]) / math.log(10) if len(scored) > 1 else float("inf")
    candidates = [list(d) for _, d in scored[:2]]
    if lod < min_lod:
        return {**out, "status": "UNCERTAIN", "alleles": None, "lod": round(lod, 2), "candidates": candidates}
    return {**out, "status": "OK", "alleles": list(best), "lod": round(lod, 2) if math.isfinite(lod) else None,
            "candidates": candidates}


def star_for_ta(n: int) -> str:
    return TA_STAR.get(int(n), f"TA{int(n)}")


def ta_label(alleles: list[int] | None) -> str:
    """'TA6/TA7 (*1/*28)'."""
    if not alleles:
        return ""
    return "/".join(f"TA{a}" for a in alleles) + " (" + "/".join(star_for_ta(a) for a in alleles) + ")"


def counts_label(counts: dict[Any, Any]) -> str:
    """'TA6 x12, TA7 x10'."""
    return ", ".join(f"TA{k} x{v}" for k, v in sorted(((int(k), v) for k, v in (counts or {}).items())))


def analyze_ta_repeat(bam_path: Path, has_chr: bool = False, reference: Path | None = None) -> dict[str, Any]:
    """Step 10 entry point: read counts plus the diplotype call."""
    try:
        reads = count_ta_repeats(bam_path, has_chr, reference=reference)
    except (OSError, ValueError) as e:
        return {"status": "ERROR", "error": str(e)}
    call = call_ta_diplotype(reads["counts"])
    return {
        **reads,
        **call,
        "locus": f"{UGT1A1_CHROM}:{UGT1A1_TA_POS} (rs3064744)",
        "diplotype": ta_label(call.get("alleles")),
        "candidate_labels": [ta_label(c) for c in call.get("candidates") or []],
    }


def _star_key(s: str) -> tuple[int, str]:
    m = re.search(r"\*(\d+)", s)
    return (int(m.group(1)) if m else 999, s)


def interpret_ugt1a1(d6: int | None, ta: dict[str, Any]) -> dict[str, Any]:
    """Medication-card entry from the BAM TA-repeat call (ta['status'] == 'OK') and the VCF *6 dosage.

    *6 (c.211G>A) sits on a TA6 haplotype, so each *6 copy replaces a *1. A *6 that has no *1 left to
    replace must share a haplotype with a decreased-function TA allele, which does not change the count of
    decreased haplotypes. *6 and *28 are not phased; trans is assumed (they rarely share a haplotype).
    """
    haps = [star_for_ta(n) for n in ta["alleles"]]
    notes: list[str] = []
    for _ in range(d6 or 0):
        if "*1" in haps:
            haps[haps.index("*1")] = "*6"
        else:
            k = next((i for i, h in enumerate(haps) if h not in ("*1", "*6") and "+" not in h), None)
            if k is None:
                k = next((i for i, h in enumerate(haps) if "+" not in h), 0)
            haps[k] = f"{haps[k]}+*6"
            notes.append("*6 shares a haplotype with a TA-repeat allele (unphased)")
    if (d6 or 0) and any(h in ("*28", "*37") for h in haps):
        notes.append("*6 and the TA-repeat allele are assumed to be on different copies")
    haps.sort(key=_star_key)
    diplotype = "/".join(haps)
    if d6 is None:
        diplotype += " (*6 unassessed)"
    unknown = [h for h in haps if h.startswith("TA")]
    n_dec = sum(1 for h in haps if any(part in DECREASED_FUNCTION for part in h.split("+")))
    basis = (f"TA repeat counted from {ta['spanning_fragments']} BAM read pairs ({counts_label(ta.get('counts'))})")
    if unknown:
        return {
            "gene": "UGT1A1", "diplotype": diplotype, "phenotype": "Indeterminate (unnamed TA-repeat length)",
            "tier": "PARTIAL", "drugs": "Irinotecan, atazanavir",
            "guidance": f"{basis}; {', '.join(unknown)} has no CPIC function assignment, so confirm with a clinical UGT1A1 test.",
            "source": "BAM_TA_REPEAT",
        }
    if n_dec >= 2:
        pheno, tier = "Poor Metabolizer", "ACTIONABLE"
        guide = ("Irinotecan: start at a reduced dose (DPWG: 70% of the standard dose; FDA label: consider at least one "
                 "dose level lower), as severe neutropenia and diarrhoea are more likely. Atazanavir (CPIC): high "
                 "likelihood of jaundice leading to stopping it; consider an alternative. This is also the usual "
                 "Gilbert syndrome genotype (harmless mild jaundice).")
    elif n_dec == 1:
        pheno, tier = "Intermediate Metabolizer", "MODERATE"
        guide = ("Irinotecan: no starting-dose change (DPWG), though neutropenia risk is modestly higher. "
                 "Atazanavir (CPIC): low likelihood of stopping it for jaundice; no need to avoid.")
    else:
        pheno, tier = "Normal Metabolizer", "NORMAL"
        guide = "Irinotecan: standard starting dose. Atazanavir (CPIC): very low likelihood of stopping it for jaundice."
    if d6 is None:
        notes.append("the *6 position lacked PASS >=10x coverage")
        if tier == "NORMAL":
            tier = "PARTIAL"
    tail = f" ({'; '.join(notes)})" if notes else ""
    return {
        "gene": "UGT1A1", "diplotype": diplotype, "phenotype": pheno, "tier": tier,
        "drugs": "Irinotecan, atazanavir",
        "guidance": f"{guide} {basis}{tail}.",
        "source": "BAM_TA_REPEAT",
    }
