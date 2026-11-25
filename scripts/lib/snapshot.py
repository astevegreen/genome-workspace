"""Annual run snapshots and the "what changed since last run" diff (Step 13).

A snapshot is one JSON file per sample per annual run, in data/snapshots/<label>/<run_date>.json. It records:
  - meta: run date, ClinVar / GenCC release dates, PharmCAT version, git commit, sample identity QC status
  - variants: every carried variant that ClinVar classed P, LP, P/LP, conflicting or VUS at the time, with its
    ClinVar VariationID, class, review status and stars, gene, condition, genotype and call QC
  - lof_rescues: Step 5's unreported ACMG loss-of-function candidates
  - medication_card, hla_associations, pharmcat: the per-gene / per-trait calls of Steps 7, 10 and 11

The snapshot holds each variant's class as ClinVar gave it that year, so ClinVar can be refreshed in place between
runs. Benign / likely benign variants are not tracked; one that turns P/LP appears as a new finding.
"""

from __future__ import annotations

import datetime
import json
from pathlib import Path
from typing import Any, Callable

from lib.genomics import parse_revstat_stars

SCHEMA_VERSION = 1

# Ordering of ClinVar germline classes, used to tell an upgrade (toward pathogenic) from a downgrade
CLASS_RANK: dict[str, float] = {
    "Benign": 0,
    "Benign/Likely benign": 0.5,
    "Likely benign": 1,
    "VUS": 2,
    "Conflicting": 2,
    "Likely pathogenic": 3,
    "Pathogenic/Likely pathogenic": 3.5,
    "Pathogenic": 4,
}
TRACKED_CLASSES = frozenset({"VUS", "Conflicting", "Likely pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic"})

SECTION_ORDER = (
    "reclassification",
    "new_finding",
    "dropped_finding",
    "review_status",
    "lof_rescue",
    "medication_card",
    "hla",
    "pharmcat",
)
SECTION_TITLES = {
    "reclassification": "ClinVar reclassifications of carried variants",
    "new_finding": "New findings",
    "dropped_finding": "Dropped findings",
    "review_status": "Review-status changes (same class)",
    "lof_rescue": "Unreported loss-of-function candidates (Step 5)",
    "medication_card": "Medication card",
    "hla": "HLA disease associations",
    "pharmcat": "PharmCAT calls",
}


def clinvar_class(clnsig: str) -> str:
    """Normalise a ClinVar CLNSIG value to one germline class.

    The part before the first '|' is the germline classification; '|risk_factor' style modifiers are ignored.
    Returns one of CLASS_RANK's keys, or 'Other: <value>' (risk factor, drug response, not provided, ...), or ''.
    """
    raw = (clnsig or "").strip()
    if not raw or raw == ".":
        return ""
    primary = raw.split("|")[0].strip().lower()
    if primary.startswith("conflicting"):
        return "Conflicting"
    tokens = [t.strip(" _") for t in primary.replace(",", "/").split("/") if t.strip(" _")]
    path = any(t.startswith("pathogenic") for t in tokens)
    lpath = any(t.startswith("likely_pathogenic") for t in tokens)
    if path and lpath:
        return "Pathogenic/Likely pathogenic"
    if path:
        return "Pathogenic"
    if lpath:
        return "Likely pathogenic"
    if any(t.startswith("uncertain_significance") for t in tokens):
        return "VUS"
    benign = "benign" in tokens
    lbenign = "likely_benign" in tokens
    if benign and lbenign:
        return "Benign/Likely benign"
    if benign:
        return "Benign"
    if lbenign:
        return "Likely benign"
    return f"Other: {raw.split('|')[0].strip()}"


def is_tracked(cls: str) -> bool:
    return cls in TRACKED_CLASSES


def class_direction(old_cls: str, new_cls: str) -> str:
    """'up' (toward pathogenic), 'down', or 'lateral' (same rank, e.g. VUS <-> conflicting, or unranked)."""
    a, b = CLASS_RANK.get(old_cls), CLASS_RANK.get(new_cls)
    if a is None or b is None or a == b:
        return "lateral"
    return "up" if b > a else "down"


def zygosity_from_gt(gt: str) -> str:
    alleles = [a for a in gt.replace("|", "/").split("/") if a not in ("", ".")]
    if not alleles:
        return "no call"
    alts = [a for a in alleles if a != "0"]
    if not alts:
        return "ref"
    if len(alleles) == 1:
        return "hemizygous"
    return "homozygous" if len(alts) == len(alleles) and len(set(alts)) == 1 else "heterozygous"


def variant_entry(
    clnsig: str,
    revstat: str,
    clinvar_id: str,
    geneinfo: str,
    condition: str,
    gt: str,
    filt: str = "",
    dp: str = "",
    gq: str = "",
) -> dict[str, Any]:
    """One tracked variant as stored in a snapshot."""
    revstat_s, stars = parse_revstat_stars(revstat)
    genes = [g.split(":")[0] for g in (geneinfo or "").split("|") if g and g != "."]
    return {
        "gene": ",".join(dict.fromkeys(genes)),
        "clinvar_id": clinvar_id if clinvar_id not in ("", ".") else "",
        "class": clinvar_class(clnsig),
        "clnsig": clnsig,
        "review_status": revstat_s.replace("_", " "),
        "stars": stars,
        "condition": (condition or "").replace("_", " ") if condition not in (None, ".") else "",
        "gt": gt,
        "zygosity": zygosity_from_gt(gt),
        "filter": filt,
        "dp": dp,
        "gq": gq,
        "qc_tier": "",
        "artifact_flag": "",
    }


def build_snapshot(
    label: str,
    run_date: str,
    meta: dict[str, Any],
    variants: dict[str, dict[str, Any]],
    lof_rescues: list[dict[str, Any]],
    card_rows: list[dict[str, str]],
    hla_rows: list[dict[str, str]],
    pharmcat_rows: list[dict[str, str]],
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "label": label,
        "run_date": run_date,
        "meta": meta,
        "variants": {k: variants[k] for k in sorted(variants)},
        "lof_rescues": {
            r["variant_id"]: {"gene": r.get("gene", ""), "consequence": r.get("consequence", r.get("most_severe_consequence", ""))}
            for r in lof_rescues if r.get("variant_id")
        },
        "medication_card": {
            r["gene"]: {"diplotype": r.get("diplotype", ""), "phenotype": r.get("phenotype", ""), "tier": r.get("tier", "")}
            for r in card_rows if r.get("gene")
        },
        "hla_associations": {
            r["trait"]: {"genotype": r.get("genotype", ""), "tier": r.get("tier", ""), "result": r.get("result", "")}
            for r in hla_rows if r.get("trait")
        },
        "pharmcat": {
            r["gene"]: {
                "diplotype": r.get("pharmcat_diplotype", ""),
                "phenotype": r.get("pharmcat_phenotype", ""),
                "status": r.get("status", ""),
            }
            for r in pharmcat_rows if r.get("gene")
        },
    }


def write_snapshot(path: Path, snap: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(json.dumps(snap, indent=1, sort_keys=True) + "\n")


def load_snapshot(path: Path) -> dict[str, Any]:
    snap = json.loads(Path(path).read_text(encoding="utf-8"))
    if snap.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"{path}: snapshot schema {snap.get('schema_version')} (this code reads {SCHEMA_VERSION})")
    return snap


def find_previous_snapshot(snap_dir: Path, run_date: str) -> Path | None:
    """Latest <date>.json in snap_dir dated strictly before run_date."""
    if not snap_dir.is_dir():
        return None
    earlier = []
    for p in snap_dir.glob("*.json"):
        try:
            datetime.date.fromisoformat(p.stem)
        except ValueError:
            continue
        if p.stem < run_date:
            earlier.append(p)
    return max(earlier, key=lambda p: p.stem) if earlier else None


def _row(section: str, key: str, gene: str, old: str, new: str, note: str = "", direction: str = "") -> dict[str, str]:
    return {"section": section, "key": key, "gene": gene, "old": old, "new": new, "direction": direction, "note": note}


def _stars(n: Any) -> str:
    return f"{n} star" if n == 1 else f"{n} stars"


def _fmt_variant(v: dict[str, Any] | None) -> str:
    if not v:
        return ""
    return f"{v['class']} ({_stars(v['stars'])})"


def _dict_changes(section: str, old: dict[str, dict], new: dict[str, dict], fields: tuple[str, ...]) -> list[dict[str, str]]:
    rows = []
    for key in sorted(set(old) | set(new)):
        a, b = old.get(key), new.get(key)
        fa = " / ".join(str(a.get(f, "")) for f in fields) if a else ""
        fb = " / ".join(str(b.get(f, "")) for f in fields) if b else ""
        if a is None:
            rows.append(_row(section, key, key, "", fb, "added"))
        elif b is None:
            rows.append(_row(section, key, key, fa, "", "removed"))
        elif fa != fb:
            changed = [f for f in fields if str(a.get(f, "")) != str(b.get(f, ""))]
            rows.append(_row(section, key, key, fa, fb, "changed: " + ", ".join(changed)))
    return rows


Lookup = Callable[[str], dict[str, Any]]


def diff_snapshots(old: dict[str, Any], new: dict[str, Any], lookup: Lookup | None = None) -> list[dict[str, str]]:
    """All changes from one snapshot to the next, as rows ordered by SECTION_ORDER.

    lookup(variant_key) -> {"class": current ClinVar class or "", "carried": bool} explains a dropped variant:
    still carried with a class below VUS is a downgrade; not carried means the VCF no longer calls it.
    """
    rows: list[dict[str, str]] = []
    ov, nv = old.get("variants", {}), new.get("variants", {})

    for key in sorted(set(ov) | set(nv)):
        a, b = ov.get(key), nv.get(key)
        gene = (b or a or {}).get("gene", "")
        if a and b:
            if a["class"] != b["class"]:
                rows.append(_row("reclassification", key, gene, _fmt_variant(a), _fmt_variant(b),
                                 b.get("condition", ""), class_direction(a["class"], b["class"])))
            elif a.get("stars") != b.get("stars") or a.get("review_status") != b.get("review_status"):
                rows.append(_row("review_status", key, gene, f"{a.get('review_status', '')} ({_stars(a.get('stars'))})",
                                 f"{b.get('review_status', '')} ({_stars(b.get('stars'))})", b["class"],
                                 "up" if b.get("stars", 0) > a.get("stars", 0) else "down"))
        elif b:
            rows.append(_row("new_finding", key, gene, "", _fmt_variant(b),
                             "Not tracked last run (then benign or likely benign, absent from ClinVar, or not called)."
                             + (f" {b['condition']}" if b.get("condition") else "")))
        else:
            info = lookup(key) if lookup else None
            if info and info.get("carried") and info.get("class"):
                cls = info["class"]
                rows.append(_row("reclassification", key, gene, _fmt_variant(a), f"{cls} (now untracked)",
                                 a.get("condition", ""), class_direction(a["class"], cls)))
            elif info and info.get("carried"):
                rows.append(_row("dropped_finding", key, gene, _fmt_variant(a), "",
                                 "Still carried, but no longer in ClinVar."))
            elif info is not None:
                rows.append(_row("dropped_finding", key, gene, _fmt_variant(a), "",
                                 "Not called in this run's VCF: check call quality or what changed in the pipeline."))
            else:
                rows.append(_row("dropped_finding", key, gene, _fmt_variant(a), "", "Not in this snapshot."))

    ol, nl = old.get("lof_rescues", {}), new.get("lof_rescues", {})
    for key in sorted(set(ol) ^ set(nl)):
        v = nl.get(key) or ol.get(key)
        rows.append(_row("lof_rescue", key, v.get("gene", ""), "" if key in nl else "present",
                         "present" if key in nl else "", "new candidate" if key in nl else "no longer a candidate"))

    rows += _dict_changes("medication_card", old.get("medication_card", {}), new.get("medication_card", {}),
                          ("diplotype", "phenotype", "tier"))
    rows += _dict_changes("hla", old.get("hla_associations", {}), new.get("hla_associations", {}),
                          ("genotype", "tier"))
    rows += _dict_changes("pharmcat", old.get("pharmcat", {}), new.get("pharmcat", {}),
                          ("diplotype", "phenotype"))

    dir_order = {"up": 0, "lateral": 1, "down": 2, "": 3}
    rows.sort(key=lambda r: (SECTION_ORDER.index(r["section"]), dir_order.get(r["direction"], 3), r["gene"], r["key"]))
    return rows


def _md_cell(s: Any) -> str:
    return str(s).replace("|", "\\|").replace("\n", " ") or "—"


def render_changes_report(label: str, old: dict[str, Any] | None, new: dict[str, Any], rows: list[dict[str, str]]) -> str:
    """Markdown 'what changed since the last run' report for one sample."""
    title = label.capitalize()
    nm = new.get("meta", {})
    md = [f"# What Changed Since the Last Annual Run — {title}", ""]
    if old is None:
        md += [
            f"**This run:** {new['run_date']} (ClinVar {nm.get('clinvar_date', '?')})",
            "",
            "This is the first snapshot for this sample, so there is nothing to compare yet. Next year's run will list "
            f"ClinVar reclassifications of the {len(new.get('variants', {}))} tracked variants (P/LP, conflicting and "
            "VUS), new or dropped findings, and medication-card, HLA and PharmCAT changes.",
            "",
        ]
        return "\n".join(md)

    om = old.get("meta", {})
    md += [
        "| | Last run | This run |",
        "|---|---|---|",
        f"| Run date | {old['run_date']} | {new['run_date']} |",
        f"| ClinVar release | {om.get('clinvar_date', '?')} | {nm.get('clinvar_date', '?')} |",
        f"| GenCC download | {om.get('gencc_date', '?')} | {nm.get('gencc_date', '?')} |",
        f"| PharmCAT | {om.get('pharmcat_version') or '—'} | {nm.get('pharmcat_version') or '—'} |",
        f"| Pipeline commit | {om.get('git_commit') or '—'} | {nm.get('git_commit') or '—'} |",
        f"| Tracked variants (P/LP, conflicting, VUS) | {len(old.get('variants', {}))} | {len(new.get('variants', {}))} |",
        "",
    ]
    if om.get("clinvar_date") and om.get("clinvar_date") == nm.get("clinvar_date"):
        md += ["> ClinVar was not refreshed between these runs, so no reclassifications are expected. "
               "Run `bash scripts/run_annual_checkup.sh --refresh-refs` for the annual comparison.", ""]

    counts = {s: sum(1 for r in rows if r["section"] == s) for s in SECTION_ORDER}
    up = sum(1 for r in rows if r["section"] == "reclassification" and r["direction"] == "up")
    md += [
        f"**Summary:** {counts['reclassification']} reclassification(s) ({up} toward pathogenic), "
        f"{counts['new_finding']} new and {counts['dropped_finding']} dropped finding(s), "
        f"{counts['medication_card']} medication-card, {counts['hla']} HLA and {counts['pharmcat']} PharmCAT change(s).",
        "",
        "> Reclassifications toward pathogenic are listed first. A new P/LP finding or an upgrade is worth raising "
        "with a genetic counselor or physician; confirm any finding in a CLIA-certified laboratory before acting on it.",
        "",
    ]
    for section in SECTION_ORDER:
        srows = [r for r in rows if r["section"] == section]
        md.append(f"## {SECTION_TITLES[section]}")
        md.append("")
        if not srows:
            md += ["No changes.", ""]
            continue
        if section in ("medication_card", "hla", "pharmcat"):
            md += ["| Gene / trait | Last run | This run | Change |", "|---|---|---|---|"]
            md += [f"| {_md_cell(r['key'])} | {_md_cell(r['old'])} | {_md_cell(r['new'])} | {_md_cell(r['note'])} |"
                   for r in srows]
        else:
            md += ["| Gene | Variant | Last run | This run | Direction | Note |", "|---|---|---|---|---|---|"]
            md += [f"| {_md_cell(r['gene'])} | `{r['key']}` | {_md_cell(r['old'])} | {_md_cell(r['new'])} | "
                   f"{_md_cell(r['direction'])} | {_md_cell(r['note'])} |" for r in srows]
        md.append("")
    return "\n".join(md)
