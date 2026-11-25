"""HLA disease associations from BAM-based T1K typing: celiac (HLA-DQ2/DQ8) and HLA-B27.

These are not diagnoses. Both are most useful as a negative: absent celiac HLA makes celiac
disease very unlikely, and absent HLA-B27 lowers the odds of ankylosing spondylitis. Carrying
either is common and usually means nothing on its own.
"""

from __future__ import annotations

from typing import Any

from lib.pgx import HLA_LOW_QUALITY, hla_slots

Slots = list[tuple[list[str], float | None]]

# Celiac heterodimers: (name, DQA1 allele prefix, DQB1 allele prefix). DQA1 and DQB1 sit next to each
# other but short reads cannot phase them; a heterodimer formed in trans confers the same risk, so
# presence of both chains counts either way.
CELIAC_HETERODIMERS: tuple[tuple[str, str, str], ...] = (
    ("DQ2.5", "HLA-DQA1*05:", "HLA-DQB1*02:01"),
    ("DQ2.2", "HLA-DQA1*02:01", "HLA-DQB1*02:02"),
    ("DQ8", "HLA-DQA1*03:", "HLA-DQB1*03:02"),
)
# Half-DQ2: one chain alone still occurs in a small share of people with celiac disease.
CELIAC_HALF_CHAINS: tuple[tuple[str, str], ...] = (
    ("DQB1*02 (DQ2 beta chain)", "HLA-DQB1*02:"),
    ("DQA1*05 (DQ2.5 alpha chain)", "HLA-DQA1*05:"),
)
# B*27 subtypes with little or no ankylosing spondylitis association.
B27_LOW_ASSOCIATION = ("HLA-B*27:06", "HLA-B*27:09")


def _alleles(slots: Slots, include_ties: bool) -> set[str]:
    """Alleles in the slots; tied candidates count only when include_ties is set."""
    return {a for cands, _ in slots if include_ties or len(cands) == 1 for a in cands}


def _has(alleles: set[str], prefix: str) -> bool:
    return any(a.startswith(prefix) for a in alleles)


def _display(gene_slots: dict[str, Slots]) -> str:
    parts = []
    for slots in gene_slots.values():
        if len(slots) == 1:
            slots = slots * 2  # homozygous (T1K reports one distinct allele)
        labels = []
        for cands, q in slots:
            text = cands[0].removeprefix("HLA-")
            if len(cands) > 1:
                text += " (tie: " + ", ".join(c.removeprefix("HLA-") for c in cands[1:]) + ")"
            if q is not None and q < HLA_LOW_QUALITY:
                text += "?"
            labels.append(text)
        parts.append("/".join(labels))
    return "; ".join(parts)


def _low_quality(gene_slots: dict[str, Slots]) -> list[str]:
    return sorted({c[0].removeprefix("HLA-") for slots in gene_slots.values() for c, q in slots if q is not None and q < HLA_LOW_QUALITY})


def _low_note(low_q: list[str], test: str) -> str:
    if not low_q:
        return ""
    return f" Low-confidence call for {', '.join(low_q)} (T1K quality < {HLA_LOW_QUALITY:g}): confirm with clinical {test} before relying on this."


def _celiac_molecules(alleles: set[str]) -> tuple[list[str], list[str]]:
    full = [name for name, a, b in CELIAC_HETERODIMERS if _has(alleles, a) and _has(alleles, b)]
    half = [] if full else [name for name, prefix in CELIAC_HALF_CHAINS if _has(alleles, prefix)]
    return full, half


def interpret_celiac_hla(hla: dict[str, Any]) -> dict[str, Any]:
    gene_slots = {g: hla_slots(hla, g) for g in ("HLA-DQA1", "HLA-DQB1")}
    entry: dict[str, Any] = {"trait": "Celiac disease (HLA-DQ2 / DQ8)", "genes": "HLA-DQA1, HLA-DQB1", "source": "BAM_T1K"}
    if not all(gene_slots.values()):
        missing = ", ".join(g for g, s in gene_slots.items() if not s)
        return {
            **entry,
            "genotype": _display({g: s for g, s in gene_slots.items() if s}) or "No confident call",
            "tier": "INCOMPLETE",
            "result": "Indeterminate",
            "meaning": f"T1K could not type {missing} confidently, so celiac HLA risk cannot be judged from this genome.",
            "next_step": "If celiac disease is ever in question, ask for clinical HLA-DQ2/DQ8 typing.",
        }
    low_note = _low_note(_low_quality(gene_slots), "HLA-DQ typing")
    full_c, half_c = _celiac_molecules(_alleles([s for v in gene_slots.values() for s in v], include_ties=False))
    full_a, half_a = _celiac_molecules(_alleles([s for v in gene_slots.values() for s in v], include_ties=True))
    entry["genotype"] = _display(gene_slots)
    if full_c or full_a:
        # A risk molecule that is only possible because of a tie must not read as absent.
        labels = full_c + [f"possible {m} (tied call)" for m in full_a if m not in full_c]
        dqb1 = gene_slots["HLA-DQB1"]
        dqb1_02_copies = sum(1 for cands, _ in (dqb1 * 2 if len(dqb1) == 1 else dqb1) if len(cands) == 1 and cands[0].startswith("HLA-DQB1*02:"))
        dose_note = (
            " Two copies of DQB1*02 (the DQ2 beta chain) is the highest-risk HLA pattern, several times the risk of one copy,"
            " though most people with it still never develop celiac disease."
            if dqb1_02_copies == 2
            else ""
        )
        return {
            **entry,
            "tier": "PRESENT",
            "result": "Carries " + ", ".join(labels) + (" (two DQB1*02 copies)" if dqb1_02_copies == 2 else ""),
            "meaning": (
                "About 95% of people with celiac disease carry DQ2 or DQ8, but so do roughly a third of people in "
                "European-ancestry populations, and only a few percent of carriers ever develop it. On its own this is not a diagnosis."
                + dose_note
                + low_note
            ),
            "next_step": "Only relevant alongside symptoms: celiac blood tests (tTG-IgA with total IgA) must be done while still eating gluten.",
        }
    if half_c or half_a:
        labels = half_c + [f"possible {m} (tied call)" for m in half_a if m not in half_c]
        return {
            **entry,
            "tier": "PARTIAL",
            "result": "Half of DQ2 only: " + ", ".join(labels),
            "meaning": (
                "No full DQ2 or DQ8 molecule. A small share of people with celiac disease carry only half of DQ2, "
                "so celiac disease is unlikely but not ruled out." + low_note
            ),
            "next_step": "If symptoms suggest celiac disease, the usual blood tests still apply.",
        }
    return {
        **entry,
        "tier": "ABSENT",
        "result": "No DQ2, DQ8 or half-DQ2",
        "meaning": (
            "Celiac disease is very unlikely: almost everyone with it (over 99%) carries DQ2, DQ8 or half of DQ2. "
            "This can rule celiac out without a gluten challenge." + low_note
        ),
        "next_step": "Mention this result if celiac testing is ever suggested.",
    }


def interpret_hla_b27(hla: dict[str, Any]) -> dict[str, Any]:
    gene_slots = {"HLA-B": hla_slots(hla, "HLA-B")}
    entry: dict[str, Any] = {"trait": "HLA-B27 (ankylosing spondylitis, uveitis)", "genes": "HLA-B", "source": "BAM_T1K"}
    if not gene_slots["HLA-B"]:
        return {
            **entry,
            "genotype": "No confident call",
            "tier": "INCOMPLETE",
            "result": "Indeterminate",
            "meaning": "T1K could not type HLA-B confidently, so HLA-B27 status is unknown.",
            "next_step": "An HLA-B27 blood test answers this if it is ever needed.",
        }
    low_note = _low_note(_low_quality(gene_slots), "HLA-B27 testing")
    certain = {a for a in _alleles(gene_slots["HLA-B"], include_ties=False) if a.startswith("HLA-B*27:")}
    every = {a for a in _alleles(gene_slots["HLA-B"], include_ties=True) if a.startswith("HLA-B*27:")}
    entry["genotype"] = _display(gene_slots)
    associated_c = sorted(a for a in certain if a not in B27_LOW_ASSOCIATION)
    associated_a = sorted(a for a in every if a not in B27_LOW_ASSOCIATION)
    if associated_c or associated_a:
        labels = [a.removeprefix("HLA-") for a in associated_c] + [
            f"possible {a.removeprefix('HLA-')} (tied call)" for a in associated_a if a not in associated_c
        ]
        return {
            **entry,
            "tier": "PRESENT",
            "result": "HLA-B27 positive: " + ", ".join(labels),
            "meaning": (
                "About 90% of people of European ancestry with ankylosing spondylitis carry HLA-B27, but only around 5% of "
                "carriers develop it (more with an affected relative). B27 is also linked to acute anterior uveitis and reactive arthritis."
                + low_note
            ),
            "next_step": "Worth mentioning if long-lasting inflammatory back pain (worse with rest, better with movement) or a painful red eye ever comes up.",
        }
    if every:
        labels = sorted(a.removeprefix("HLA-") for a in every)
        return {
            **entry,
            "tier": "PARTIAL",
            "result": "HLA-B27 low-association subtype: " + ", ".join(labels),
            "meaning": "B*27:06 and B*27:09 show little or no link to ankylosing spondylitis, unlike other B27 subtypes." + low_note,
            "next_step": "A standard HLA-B27 blood test reads these as positive; mention the subtype if that happens.",
        }
    return {
        **entry,
        "tier": "ABSENT",
        "result": "HLA-B27 negative",
        "meaning": "No HLA-B27 allele. Ankylosing spondylitis is less likely, though it still occurs without B27." + low_note,
        "next_step": "None.",
    }


def interpret_hla_associations(hla: dict[str, Any] | None) -> list[dict[str, Any]]:
    """Celiac (DQ2/DQ8) and HLA-B27 entries from T1K typing; empty when T1K did not run."""
    if not hla or hla.get("status") != "OK":
        return []
    return [interpret_celiac_hla(hla), interpret_hla_b27(hla)]
