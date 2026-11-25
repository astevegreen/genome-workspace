"""PharmCAT cross-check of the medication card: input preparation, outside calls, and comparison.

PharmCAT treats any position absent from its input VCF as a no-call. A variants-only VCF therefore
leaves most CPIC positions "missing". The sample gVCF fills that gap: a PharmCAT position inside a
0/0 reference block with enough depth is written as an explicit 0/0 call, and every other absent
position stays absent (a no-call), which is safer than PharmCAT's --missing-to-ref flag.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pysam

from lib.pgx import hla_slots

REF_BLOCK_MIN_DP = 10  # same 10x floor the card uses for G6PD reference calls
HEMIZYGOUS_MIN_DP = 6  # male non-PAR chrX: one copy, so half the depth (Step 4 HEMI threshold)
GRCH38_X_PARS = ((10_001, 2_781_479), (155_701_383, 156_030_895))
SYMBOLIC_REF_ALTS = {"<*>", "<NON_REF>", "<X>"}

# Genes PharmCAT reports without a phenotype; compare their diplotypes instead.
DIPLOTYPE_ONLY_GENES = {"CYP4F2", "VKORC1"}
# Genes whose PharmCAT result comes from an outside call we supply (so only the phenotype translation is tested).
OUTSIDE_CALL_GENES = {"CYP2D6", "HLA-A", "HLA-B", "MT-RNR1"}


def _sample_contig(chrom: str, sample_contigs: set[str]) -> str | None:
    """Map a PharmCAT (chr-prefixed) contig to the sample's naming."""
    if chrom in sample_contigs:
        return chrom
    bare = chrom.removeprefix("chr")
    for cand in (bare, "MT" if bare == "M" else bare, "chrM" if bare in ("M", "MT") else None):
        if cand and cand in sample_contigs:
            return cand
    return None


def _is_ref_block(rec: pysam.VariantRecord) -> bool:
    if rec.alts and any(a not in SYMBOLIC_REF_ALTS for a in rec.alts):
        return False
    gt = rec.samples[0].get("GT")
    return bool(gt) and all(a == 0 for a in gt)


def _block_depth(rec: pysam.VariantRecord) -> int | None:
    """Minimum depth of a PASS reference block (INFO/MinDP, FORMAT/MIN_DP, then FORMAT/DP); None if filtered."""
    if any(f != "PASS" for f in rec.filter.keys()):
        return None
    if "MinDP" in rec.info:
        return int(rec.info["MinDP"])
    for key in ("MIN_DP", "DP"):
        if key in rec.format and rec.samples[0].get(key) is not None:
            return int(rec.samples[0][key])
    return None


def confirmed_reference_positions(
    gvcf_path: Path, positions_path: Path, min_dp: int = REF_BLOCK_MIN_DP, male: bool = False
) -> list[tuple[str, int, str, tuple[str, ...]]]:
    """PharmCAT positions the gVCF shows as homozygous reference across their full REF span.

    Returns (sample_contig, pos_1based, ref, alts) for each position where no variant record overlaps
    and 0/0 blocks of depth >= min_dp cover every base of the REF allele. For a male sample, non-PAR chrX
    positions use the hemizygous floor instead.
    """
    out: list[tuple[str, int, str, tuple[str, ...]]] = []
    with pysam.VariantFile(str(gvcf_path)) as gvcf, pysam.VariantFile(str(positions_path)) as positions:
        contigs = set(gvcf.header.contigs)
        for p in positions:
            contig = _sample_contig(p.chrom, contigs)
            if contig is None or not p.ref:
                continue
            start, end = p.start, p.start + len(p.ref)  # 0-based half-open span of the REF allele
            floor = min_dp
            if male and contig in ("X", "chrX") and not any(a <= p.pos <= b for a, b in GRCH38_X_PARS):
                floor = min(min_dp, HEMIZYGOUS_MIN_DP)
            try:
                overlapping = list(gvcf.fetch(contig, start, end))
            except ValueError:
                continue
            if not overlapping or any(not _is_ref_block(r) for r in overlapping):
                continue
            covered_to = start
            for r in sorted(overlapping, key=lambda r: r.start):
                dp = _block_depth(r)
                if r.start > covered_to or dp is None or dp < floor:
                    break
                covered_to = max(covered_to, r.stop)
            if covered_to >= end:
                out.append((contig, p.pos, p.ref, tuple(p.alts or ())))
    return out


def _minimal_header(src: pysam.VariantHeader) -> pysam.VariantHeader:
    """Contigs, FILTERs, GT and the sample only.

    PharmCAT's Java VCF parser rejects backslash-escaped quotes in header descriptions (vendor bcftools
    filter lines carry them), and it reads nothing beyond GT, so the source header is not copied.
    """
    header = pysam.VariantHeader()
    for name, contig in src.contigs.items():
        header.contigs.add(name, length=contig.length)
    for name, flt in src.filters.items():
        if name != "PASS":
            desc = (flt.description or name).replace("\\", "").replace('"', "'")
            header.filters.add(name, None, None, desc)
    header.formats.add("GT", 1, "String", "Genotype")
    header.add_sample(list(src.samples)[0])
    return header


def write_pharmcat_input(
    variants_vcf: Path,
    gvcf_path: Path,
    positions_path: Path,
    out_vcf: Path,
    min_dp: int = REF_BLOCK_MIN_DP,
    pad: int = 50,
    male: bool = False,
) -> dict[str, int]:
    """Write a bgzipped, indexed VCF of the sample's PGx-region variants plus gVCF-confirmed 0/0 calls.

    Only CHROM, POS, ID, REF, ALT, QUAL, FILTER and GT are kept; PharmCAT reads nothing else.
    """
    fills = confirmed_reference_positions(gvcf_path, positions_path, min_dp, male=male)
    with pysam.VariantFile(str(variants_vcf)) as vin, pysam.VariantFile(str(positions_path)) as positions:
        contigs = list(vin.header.contigs)
        contig_set = set(contigs)
        order = {c: i for i, c in enumerate(contigs)}
        spans: list[tuple[str, int, int]] = []
        for p in positions:
            c = _sample_contig(p.chrom, contig_set)
            if c:
                spans.append((c, max(0, p.start - pad), p.start + len(p.ref or "N") + pad))

        # (contig, 0-based start, alleles) -> (ID, QUAL, filters, GT)
        rows: dict[tuple[str, int, tuple[str, ...]], tuple[Any, ...]] = {}
        for c, s, e in spans:
            for r in vin.fetch(c, s, e):
                key = (r.chrom, r.start, tuple(r.alleles))
                rows.setdefault(key, (r.id, r.qual, list(r.filter.keys()), r.samples[0].get("GT")))
        variant_starts = {(c, s) for c, s, _ in rows}
        n_variants = len(rows)
        n_fills = 0
        for contig, pos, ref, alts in fills:
            if (contig, pos - 1) in variant_starts or not alts:
                continue
            rows[(contig, pos - 1, (ref, *alts))] = (None, None, [], (0, 0))
            n_fills += 1

        header = _minimal_header(vin.header)
        sample = list(header.samples)[0]
        tmp = out_vcf.parent / (out_vcf.name.removesuffix(".bgz").removesuffix(".gz") + ".tmp.vcf")
        with pysam.VariantFile(str(tmp), "w", header=header) as vout:
            for (contig, start, alleles), (vid, qual, filters, gt) in sorted(
                rows.items(), key=lambda kv: (order.get(kv[0][0], 10**6), kv[0][1], kv[0][2])
            ):
                rec = vout.new_record(contig=contig, start=start, alleles=alleles, id=vid, qual=qual)
                for f in filters:
                    rec.filter.add(f)
                rec.samples[sample]["GT"] = gt
                vout.write(rec)
    pysam.tabix_compress(str(tmp), str(out_vcf), force=True)
    tmp.unlink()
    pysam.tabix_index(str(out_vcf), preset="vcf", force=True, csi=True)
    return {"variant_records": n_variants, "reference_fills": n_fills, "positions_confirmed_ref": len(fills)}


def cyrius_to_outside(genotype: str | None) -> str | None:
    """Convert a Cyrius CYP2D6 call (e.g. '*1/*36+*10', '*1x3/*2') to PharmCAT outside-call syntax."""
    gt = (genotype or "").strip()
    if not gt or gt.lower() == "none" or ";" in gt or gt.count("/") != 1:
        return None
    haps = []
    for hap in gt.split("/"):
        hap = re.sub(r"x([3-9]|\d{2,})$", "≥3", hap)  # PharmCAT names 3+ copies '*1≥3'
        if "+" in hap:
            hap = "[" + " + ".join(h.strip() for h in hap.split("+")) + "]"
        haps.append(hap)
    return "/".join(haps)


def hla_to_outside(hla: dict[str, Any] | None, gene: str) -> str | None:
    """Two-field T1K diplotype for an outside call ('*57:01/*07:02'); None when a slot is tied or untyped."""
    if not hla or hla.get("status") != "OK":
        return None
    slots = hla_slots(hla, gene)
    if not slots or any(len(c) > 1 for c, _ in slots):
        return None
    alleles = [c[0].removeprefix(gene) for c, _ in slots]
    if len(alleles) == 1:
        alleles *= 2
    return "/".join(alleles[:2])


def build_outside_calls(bam: dict[str, Any] | None) -> dict[str, str]:
    """Outside calls from Step 10: Cyrius CYP2D6 and untied T1K HLA-A / HLA-B."""
    bam = bam or {}
    calls: dict[str, str] = {}
    cyr = bam.get("cyrius_cyp2d6") or {}
    if cyr.get("status") == "OK" and (d := cyrius_to_outside(cyr.get("genotype"))):
        calls["CYP2D6"] = d
    for gene in ("HLA-A", "HLA-B"):
        if d := hla_to_outside(bam.get("hla_t1k"), gene):
            calls[gene] = d
    return calls


def write_outside_calls(path: Path, calls: dict[str, str]) -> None:
    lines = ["# gene\tdiplotype (outside calls from the Step 10 BAM callers)"]
    lines += [f"{g}\t{d}" for g, d in sorted(calls.items())]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_calls_only_tsv(path: Path) -> dict[str, dict[str, str]]:
    """Parse PharmCAT's -reporterCallsOnlyTsv output into {gene: row}."""
    rows: dict[str, dict[str, str]] = {}
    header: list[str] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        cols = line.split("\t")
        if header is None:
            if cols[0].lstrip("#").strip() == "Gene":
                header = [c.lstrip("#").strip() for c in cols]
            continue
        row = dict(zip(header, cols + [""] * (len(header) - len(cols))))
        if row.get("Gene"):
            rows[row["Gene"]] = row
    return rows


def _fmt_score(v: Any) -> str:
    if v in (None, "", "n/a", "No Result"):
        return ""
    try:
        return f"{float(v):g}" if float(v) != int(float(v)) else f"{float(v):.1f}"
    except (TypeError, ValueError):
        return str(v)


def parse_phenotype_json(path: Path, tsv_rows: dict[str, dict[str, str]] | None = None) -> dict[str, dict[str, str]]:
    """Read PharmCAT's <name>.phenotype.json into {gene: row}, shaped like the calls-only TSV rows.

    The JSON is the fuller source: the calls-only TSV leaves out genes called from outside calls (HLA)
    and shows DPYD as its unphased variant list instead of the lowest-function pair CPIC dosing uses.
    The recommendation diplotypes are read; "Missing positions" (a yes/no flag) comes from the TSV.
    """
    data = json.loads(path.read_text(encoding="utf-8"))
    rows: dict[str, dict[str, str]] = {}
    for gene, rep in (data.get("geneReports") or {}).items():
        recs = rep.get("recommendationDiplotypes") or []
        labels = list(dict.fromkeys(r.get("label") or "" for r in recs if r.get("label")))
        # one diplotype can carry several phenotypes (HLA-B: one per risk allele); several diplotypes are alternatives
        per_dip = (", ".join(p for p in (r.get("phenotypes") or []) if p) for r in recs)
        phenos = [p for p in dict.fromkeys(per_dip) if p]
        scores = list(dict.fromkeys(s for r in recs if (s := _fmt_score(r.get("activityScore")))))
        rows[gene] = {
            "Gene": gene,
            "Source Diplotype": " OR ".join(labels),
            "Phenotype": " OR ".join(phenos),
            "Activity Score": " OR ".join(scores),
            "Call Source": str(rep.get("callSource") or ""),
            "Missing positions": ((tsv_rows or {}).get(gene) or {}).get("Missing positions", ""),
        }
    return rows


_PHENOTYPE_CLASSES: tuple[tuple[str, str], ...] = (
    ("ultrarapid", "Ultrarapid"),
    ("rapid", "Rapid"),
    ("poor metabolizer", "Poor metabolizer"),
    ("poor function", "Poor function"),
    ("intermediate", "Intermediate"),
    ("normal metabolizer", "Normal metabolizer"),
    ("decreased function", "Decreased function"),
    ("increased function", "Increased function"),
    ("normal function", "Normal function"),
    ("deficient", "Deficient"),
    ("variable", "Variable"),
    ("positive", "Positive"),
    ("negative", "Negative"),
    ("normal", "Normal"),
)
_NO_CALL = {"", "n/a", "no result", "unknown", "indeterminate", "no call", "not called"}


def phenotype_class(text: str | None) -> str | None:
    """Reduce a card or PharmCAT phenotype to a comparable class; None for a no-call."""
    t = (text or "").strip().lower()
    if t in _NO_CALL or t.startswith("indeterminate") or t.startswith("unknown"):
        return None
    for needle, label in _PHENOTYPE_CLASSES:
        if needle in t:
            return label
    return None


def _stars(diplotype: str | None) -> list[str]:
    return sorted(re.findall(r"\*\d+[A-Z]?(?:x\d|≥\d)?", diplotype or ""))


# CYP4F2 haplotypes carrying rs2108622 T (V433M), the variant CPIC's warfarin guideline uses. PharmCAT
# names V433M with W12G as *4 and with rs150579280 as *22, while the card checks rs2108622 alone (*3).
CYP4F2_V433M_ALLELES = {"*3", "*4", "*22"}


def _cyp4f2_v433m_count(diplotype: str | None) -> int | None:
    stars = _stars(diplotype)
    return None if len(stars) != 2 else sum(1 for s in stars if s in CYP4F2_V433M_ALLELES)


def _vkorc1_alt_count(text: str | None, pharmcat: bool) -> int | None:
    t = text or ""
    if pharmcat:
        return t.count("variant") if "rs9923231" in t else None
    m = re.search(r"-1639 ([CT])/([CT])", t)
    return None if not m else sum(1 for a in m.groups() if a == "T")


def compare_card_with_pharmcat(
    card_rows: list[dict[str, Any]], pc_rows: dict[str, dict[str, str]], outside_calls: dict[str, str]
) -> list[dict[str, Any]]:
    """One row per gene: card result, PharmCAT result, and a status."""
    card = {str(r["gene"]): r for r in card_rows}
    out: list[dict[str, Any]] = []
    for gene in sorted(set(card) | set(pc_rows)):
        c, p = card.get(gene), pc_rows.get(gene)
        row: dict[str, Any] = {
            "gene": gene,
            "card_diplotype": str((c or {}).get("diplotype") or ""),
            "card_phenotype": str((c or {}).get("phenotype") or ""),
            "pharmcat_diplotype": (p or {}).get("Source Diplotype", ""),
            "pharmcat_phenotype": (p or {}).get("Phenotype", ""),
            "pharmcat_activity_score": (p or {}).get("Activity Score", ""),
            "pharmcat_missing_positions": (p or {}).get("Missing positions", ""),
            "note": "",
        }
        if gene in outside_calls:
            row["note"] = f"PharmCAT used our BAM call ({outside_calls[gene]}) as input, so this checks the phenotype translation only."
        elif gene == "MT-RNR1" and p is not None:
            row["note"] = "PharmCAT reads MT-RNR1 only from an outside call; the card's own VCF call is not passed to it."
        elif gene in OUTSIDE_CALL_GENES and p is not None:
            row["note"] = "PharmCAT needs an outside call for this gene; none was passed (no confident BAM call)."
        if c is None:
            row["status"] = "PHARMCAT_ONLY"
        elif p is None:
            row["status"] = "CARD_ONLY"
        elif gene == "VKORC1":
            a = _vkorc1_alt_count(row["card_diplotype"], pharmcat=False)
            b = _vkorc1_alt_count(row["pharmcat_diplotype"], pharmcat=True)
            row["status"] = _diplotype_status(a, b)
        elif gene == "CYP4F2":
            a = _cyp4f2_v433m_count(row["card_diplotype"])
            b = _cyp4f2_v433m_count(row["pharmcat_diplotype"]) if " OR " not in row["pharmcat_diplotype"] else None
            row["status"] = _diplotype_status(a, b)
            if row["status"] == "AGREE" and _stars(row["card_diplotype"]) != _stars(row["pharmcat_diplotype"]):
                row["status"] = "DIPLOTYPE_DIFFERS"
                row["note"] = (f"Both show {a} copy(ies) of rs2108622 (V433M), the variant the warfarin guideline uses; "
                               "PharmCAT also reads other CYP4F2 positions the card does not check.")
        elif gene in DIPLOTYPE_ONLY_GENES:
            row["status"] = _diplotype_status(_stars(row["card_diplotype"]) or None, _stars(row["pharmcat_diplotype"]) or None)
        else:
            card_cls = phenotype_class(row["card_phenotype"])
            pc_cls = phenotype_class(row["pharmcat_phenotype"])
            if str(c.get("tier", "")).upper() in ("INCOMPLETE", "PARTIAL") or card_cls is None:
                row["status"] = "CARD_INCOMPLETE"
            elif pc_cls is None or " OR " in row["pharmcat_phenotype"]:
                row["status"] = "PHARMCAT_NO_CALL"
            elif card_cls != pc_cls:
                row["status"] = "DISAGREE"
            else:
                a, b = _stars(row["card_diplotype"]), _stars(row["pharmcat_diplotype"])
                row["status"] = "DIPLOTYPE_DIFFERS" if a and b and a != b else "AGREE"
        if row["status"] == "PHARMCAT_NO_CALL" and " OR " in row["pharmcat_diplotype"]:
            row["note"] = (row["note"] + " PharmCAT found more than one possible diplotype.").strip()
        out.append(row)
    return out


def _diplotype_status(card_value: Any, pharmcat_value: Any) -> str:
    if pharmcat_value is None:
        return "PHARMCAT_NO_CALL"
    if card_value is None:
        return "CARD_INCOMPLETE"
    return "AGREE" if card_value == pharmcat_value else "DISAGREE"


STATUS_ORDER = ("DISAGREE", "DIPLOTYPE_DIFFERS", "CARD_INCOMPLETE", "PHARMCAT_NO_CALL", "PHARMCAT_ONLY", "CARD_ONLY", "AGREE")
STATUS_TEXT = {
    "DISAGREE": "Disagree: check before relying on the card",
    "DIPLOTYPE_DIFFERS": "Same phenotype, different star alleles",
    "CARD_INCOMPLETE": "Card incomplete; PharmCAT may fill the gap",
    "PHARMCAT_NO_CALL": "PharmCAT made no call",
    "PHARMCAT_ONLY": "Only PharmCAT covers this gene",
    "CARD_ONLY": "Only the card covers this gene",
    "AGREE": "Agree",
}


def render_crosscheck_report(sample: str, rows: list[dict[str, Any]], meta: dict[str, Any]) -> str:
    counts = {s: sum(1 for r in rows if r["status"] == s) for s in STATUS_ORDER}
    md = [
        f"# PharmCAT Cross-Check of the Medication Card — {sample.title()}",
        "",
        f"**Run date:** {meta.get('run_date', '')} | **PharmCAT:** {meta.get('pharmcat_version', 'unknown')} | "
        f"**Input:** {meta.get('variant_records', 0)} PGx-region variant records + {meta.get('reference_fills', 0)} "
        f"gVCF-confirmed reference calls (0/0 blocks at >= {meta.get('min_dp', REF_BLOCK_MIN_DP)}x)",
        "",
        "PharmCAT is the reference CPIC caller: it builds full star-allele diplotypes from every defining position, "
        "while the card checks key variants one at a time. Where the two disagree, the card is the more likely to be wrong.",
        "",
        "**Summary:** " + ", ".join(f"{s.replace('_', ' ').lower()}: {counts[s]}" for s in STATUS_ORDER if counts[s]),
        "",
        "| Gene | Card | PharmCAT | Status | Note |",
        "|---|---|---|---|---|",
    ]
    for s in STATUS_ORDER:
        for r in (r for r in rows if r["status"] == s):
            card = " — ".join(b for b in (r["card_diplotype"], r["card_phenotype"]) if b) or "—"
            pc_bits = [b for b in (r["pharmcat_diplotype"], r["pharmcat_phenotype"]) if b]
            if r["pharmcat_activity_score"]:
                pc_bits.append(f"AS={r['pharmcat_activity_score']}")
            pc = " — ".join(pc_bits) or "—"
            note = r["note"]
            if r["pharmcat_missing_positions"].strip().lower() == "yes":
                note = (note + " Some defining positions had no call.").strip()
            md.append(f"| **{r['gene']}** | {card} | {pc} | **{STATUS_TEXT[s]}** | {note} |")
    md.extend([
        "",
        "_Research-grade cross-check. Positions with no gVCF evidence stay no-calls instead of being assumed reference, "
        "so PharmCAT may report fewer calls than a clinical assay. Confirm any actionable or disputed gene with clinical "
        "PGx testing. PharmCAT's own report is in the `pharmcat/` folder next to this file._",
    ])
    return "\n".join(md) + "\n"
