"""Shared clinical variant parsing, interval math, coverage auditing, and artifact flagging."""

from __future__ import annotations

import csv
import datetime
import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
import pysam

from .constants import (
    ACMG_SF_GENES,
    ACMG_SF_VERSION,
    AD_BS1_AF,
    BA1_COMMON_AF,
    CLINVAR_STARS,
)


def parse_revstat_stars(revstat_val: Any) -> tuple[str, int]:
    """Return (revstat_str, star_count) for a ClinVar CLNREVSTAT INFO value."""
    if revstat_val is None:
        return ("", 0)
    if isinstance(revstat_val, (tuple, list)):
        s = ",".join(str(x) for x in revstat_val)
    else:
        s = str(revstat_val)
    return (s, CLINVAR_STARS.get(s, 0))


def parse_info_str(val: Any) -> str:
    """Convert a pysam INFO value (scalar, tuple, list, or None) to string."""
    if val is None:
        return ""
    if isinstance(val, (tuple, list)):
        return ",".join(str(x) for x in val)
    return str(val)


def classify_clnsig(clnsig_str: str, clnsigconf_str: str = "") -> str | None:
    """Classify a ClinVar CLNSIG string into 'PLP', 'CONFLICTING_PLP', 'RISK_OR_DRUG', or None."""
    low = clnsig_str.lower()
    if not low or low == ".":
        return None
    if "conflicting" not in low and ("pathogenic" in low or "likely_pathogenic" in low):
        tokens = [t.strip() for t in low.replace("|", "/").replace(",", "/").split("/")]
        if any(t in ("pathogenic", "likely_pathogenic") for t in tokens):
            return "PLP"
    if "conflicting" in low and "pathogenic" in clnsigconf_str.lower():
        return "CONFLICTING_PLP"
    if any(k in low for k in ("drug_response", "risk_factor", "protective", "association")):
        return "RISK_OR_DRUG"
    return None


def extract_genes(geneinfo: str) -> list[str]:
    """Extract unique HGNC gene symbols from ClinVar GENEINFO string."""
    if not geneinfo:
        return []
    genes: list[str] = []
    for item in geneinfo.split("|"):
        sym = item.split(":")[0].strip()
        if sym and sym not in genes:
            genes.append(sym)
    return genes


def resolve_gene_moi(
    genes: list[str],
    gencc_moi: dict[str, dict[str, Any]],
    chrom: str = "",
) -> tuple[str, str, str]:
    """
    Resolve primary gene, Mode of Inheritance (MOI), and GenCC diseases across
    all genes listed in ClinVar GENEINFO. Prefers the first gene with a curated
    GenCC MOI rather than blindly stopping at genes[0].
    """
    primary_gene = genes[0] if genes else ""
    moi = "UNKNOWN"
    gencc_diseases = ""

    for g in genes:
        entry = gencc_moi.get(g)
        if entry and entry.get("moi", "UNKNOWN") != "UNKNOWN":
            primary_gene = g
            moi = entry["moi"]
            gencc_diseases = entry.get("gencc_diseases", "")
            break
    else:
        if primary_gene and primary_gene in gencc_moi:
            gencc_diseases = gencc_moi[primary_gene].get("gencc_diseases", "")

    if chrom == "MT" and moi == "UNKNOWN":
        moi = "MT"
    elif chrom == "X" and moi == "UNKNOWN":
        moi = "XL"

    return primary_gene, moi, gencc_diseases


def read_sample_sex(samples_tsv: Path, sample_label: str) -> str:
    """Look up sample sex ('M'/'F', or '' if unknown) from a config/samples.tsv registry."""
    if not samples_tsv.exists():
        return ""
    with open(samples_tsv, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("label") == sample_label:
                return (row.get("sex") or "").strip().upper()
    return ""


def normalize_gt_unphased(gt_str: str) -> str:
    """Normalize phased genotypes (e.g. '1|0' -> '0/1', '1|1' -> '1/1') for diplotype lookup."""
    s = (gt_str or "").strip()
    if "|" not in s:
        return s
    parts = s.split("|")
    if len(parts) == 2 and all(p.isdigit() for p in parts):
        return "/".join(sorted(parts, key=int))
    return s.replace("|", "/")


def derive_apoe_diplotype(gt_429: str, gt_741: str) -> str:
    """Derive APOE diplotype from rs429358 and rs7412 genotypes (phased or unphased)."""
    g429 = normalize_gt_unphased(gt_429)
    g741 = normalize_gt_unphased(gt_741)
    if g429 == "0/0":
        g429 = "0/0 (REF)"
    if g741 == "0/0":
        g741 = "0/0 (REF)"
    apoe_map = {
        ("0/0 (REF)", "0/0 (REF)"): "ε3/ε3 (most common / neutral reference)",
        ("0/0 (REF)", "0/1"): "ε2/ε3",
        ("0/0 (REF)", "1/1"): "ε2/ε2",
        ("0/1", "0/0 (REF)"): "ε3/ε4",
        ("1/1", "0/0 (REF)"): "ε4/ε4",
        ("0/1", "0/1"): "ε2/ε4 (unphased; rarely ε1/ε3)",
    }
    return apoe_map.get((g429, g741), f"Unresolved (rs429358={gt_429}, rs7412={gt_741})")


def is_strong_conflicting(row: dict[str, Any]) -> bool:
    """Return True if a Conflicting_PLP variant has meaningful P/LP submitter support."""
    conf = str(row.get("clnsigconf") or "")
    p_cnt = sum(int(m) for m in re.findall(r"(?:Pathogenic|Likely_pathogenic)_\((\d+)\)", conf))
    b_cnt = sum(int(m) for m in re.findall(r"(?:Benign|Likely_benign)_\((\d+)\)", conf))
    return p_cnt >= 2 and p_cnt >= b_cnt


def assess_artifact_or_common(row: dict[str, Any]) -> tuple[str, str]:
    """
    Apply automated ACMG BA1/BS1 & low-evidence checks to flag likely sequencing
    artifacts, common benign polymorphisms, blood-group traits, or 0-star legacy hits.

    Returns (flag, reason) where flag is 'OK' or 'LIKELY_ARTIFACT_OR_COMMON'.
    """
    af = row.get("gnomad_af")
    stars = int(row.get("clinvar_stars") or 0)
    moi = str(row.get("moi") or "UNKNOWN")
    clndn = str(row.get("clndn") or "").lower()
    reasons: list[str] = []

    if "blood group" in clndn or "blood_group" in clndn:
        reasons.append("Benign blood-group trait")

    if af is not None:
        if af >= BA1_COMMON_AF:
            reasons.append(f"BA1 common AF ({af:.2%} >= {BA1_COMMON_AF:.0%}; likely artifact or benign polymorphism)")
        elif moi == "AD" and af >= AD_BS1_AF:
            reasons.append(f"BS1 high AF for AD ({af:.2%} >= {AD_BS1_AF:.0%})")
        elif stars == 0 and af >= 0.005:
            reasons.append(f"0-star ClinVar entry with elevated AF ({af:.2%})")

    if stars == 0 and any(k in clndn for k in ("carcinoma", "melanoma", "leukemia", "lymphoma", "glioma")) and moi == "UNKNOWN":
        reasons.append("0-star legacy somatic tumor entry without curated germline MOI")

    if reasons:
        return ("LIKELY_ARTIFACT_OR_COMMON", "; ".join(reasons))
    return ("OK", "")


def load_gencc_moi(gencc_path: Path) -> dict[str, dict[str, Any]]:
    """Aggregate GenCC curated Mode of Inheritance per gene symbol."""
    if not gencc_path.exists():
        return {}
    df = pd.read_csv(gencc_path, sep="\t", dtype=str, low_memory=False)
    valid_strengths = {"Definitive", "Strong", "Moderate", "Supportive", "Limited"}
    gene_info: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"mois": set(), "diseases": set(), "classifications": set()}
    )
    for _, row in df.iterrows():
        gene = str(row.get("gene_symbol", "")).strip()
        cls = str(row.get("classification_title", "")).strip()
        moi = str(row.get("moi_title", "")).strip()
        dis = str(row.get("disease_title", "")).strip()
        if not gene or cls not in valid_strengths:
            continue
        gene_info[gene]["classifications"].add(cls)
        if dis and dis != "nan":
            gene_info[gene]["diseases"].add(dis)
        moi_low = moi.lower()
        if "autosomal recessive" in moi_low:
            gene_info[gene]["mois"].add("AR")
        elif "autosomal dominant" in moi_low:
            gene_info[gene]["mois"].add("AD")
        elif "x-linked" in moi_low:
            gene_info[gene]["mois"].add("XL")
        elif "mitochondrial" in moi_low:
            gene_info[gene]["mois"].add("MT")
        elif "semidominant" in moi_low:
            gene_info[gene]["mois"].add("AD/AR")
        elif "y-linked" in moi_low:
            gene_info[gene]["mois"].add("YL")

    summary: dict[str, dict[str, Any]] = {}
    for g, d in gene_info.items():
        mois = sorted(d["mois"])
        if "AD/AR" in mois or ("AD" in mois and "AR" in mois):
            moi_cat = "AD/AR"
        elif "AD" in mois:
            moi_cat = "AD"
        elif "AR" in mois:
            moi_cat = "AR"
        elif "XL" in mois:
            moi_cat = "XL"
        elif "MT" in mois:
            moi_cat = "MT"
        elif "YL" in mois:
            moi_cat = "YL"
        else:
            moi_cat = "UNKNOWN"
        summary[g] = {
            "moi": moi_cat,
            "mois_raw": ",".join(mois) if mois else "UNKNOWN",
            "gencc_diseases": "; ".join(sorted(d["diseases"])[:5]),
        }
    return summary


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or adjacent half-open [start, end) intervals."""
    if not intervals:
        return []
    s = sorted(intervals)
    merged = [s[0]]
    for start, end in s[1:]:
        p_start, p_end = merged[-1]
        if start <= p_end:
            merged[-1] = (p_start, max(p_end, end))
        else:
            merged.append((start, end))
    return merged


def reciprocal_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    """Return min(ov/len_a, ov/len_b) for two half-open intervals."""
    ov = max(0, min(a_end, b_end) - max(a_start, b_start))
    if ov == 0:
        return 0.0
    len_a = max(1, a_end - a_start)
    len_b = max(1, b_end - b_start)
    return min(ov / len_a, ov / len_b)


def load_gencode_cds(
    cds_bed_gz: Path,
) -> tuple[dict[str, list[tuple[int, int, str]]], dict[str, tuple[str, list[tuple[int, int]]]]]:
    """Load CDS intervals by chrom (for SV/CNV overlap) and merged per gene (for coverage audit)."""
    by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    by_gene_raw: dict[str, tuple[str, list[tuple[int, int]]]] = {}

    with gzip.open(cds_bed_gz, "rt", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start_s, end_s, gene = parts[0], int(parts[1]), int(parts[2]), parts[3]
            by_chrom[chrom].append((start_s, end_s, gene))
            if gene not in by_gene_raw:
                by_gene_raw[gene] = (chrom, [])
            by_gene_raw[gene][1].append((start_s, end_s))

    by_gene_merged = {
        g: (chrom, merge_intervals(ivs)) for g, (chrom, ivs) in by_gene_raw.items()
    }
    return by_chrom, by_gene_merged


def record_is_pass(rec: pysam.VariantRecord) -> bool:
    """Return True if VCF record FILTER is PASS or empty ('.')."""
    return len(rec.filter) == 0 or "PASS" in rec.filter


def record_effective_min_dp(rec: pysam.VariantRecord, sample: Any) -> int:
    """
    Extract conservative minimum depth for a gVCF block or variant record.
    Prefers MinDP (INFO) / MIN_DP (FORMAT) when present on reference blocks,
    falling back to sample DP or INFO DP.
    """
    min_dp = rec.info.get("MinDP")
    if min_dp is None:
        min_dp = sample.get("MIN_DP")
    dp = sample.get("DP")
    if dp is None:
        dp = rec.info.get("DP")
    if min_dp is not None and dp is not None:
        return int(min(min_dp, dp))
    if min_dp is not None:
        return int(min_dp)
    if dp is not None:
        return int(dp)
    return 0


def compute_gene_gvcf_coverage(
    gvcf: pysam.VariantFile,
    chrom: str,
    exons: list[tuple[int, int]],
    min_pass_dp: int = 10,
) -> dict[str, Any]:
    """
    Compute per-gene CDS coverage from an open gVCF:
      - pct_ge_10x: % of CDS bases covered by FILTER==PASS records with effective MinDP >= 10
      - pct_ge_1x:  % of CDS bases with any coverage (DP >= 1, regardless of FILTER)
      - pct_0x_nocall: % of CDS bases with 0x / no gVCF record
    """
    total_bp = sum(e - s for s, e in exons)
    bp_ge10_pass = 0
    bp_ge1 = 0

    for ex_s, ex_e in exons:
        ivs_10_pass: list[tuple[int, int]] = []
        ivs_1: list[tuple[int, int]] = []
        for rec in gvcf.fetch(chrom, ex_s, ex_e):
            r_start = rec.pos - 1
            r_end = rec.stop or rec.info.get("END") or rec.pos
            ov_s = max(ex_s, r_start)
            ov_e = min(ex_e, r_end)
            if ov_e <= ov_s:
                continue
            sample = next(iter(rec.samples.values()))
            eff_dp = record_effective_min_dp(rec, sample)
            if eff_dp >= 1:
                ivs_1.append((ov_s, ov_e))
            if eff_dp >= min_pass_dp and record_is_pass(rec):
                ivs_10_pass.append((ov_s, ov_e))
        bp_ge1 += sum(e - s for s, e in merge_intervals(ivs_1))
        bp_ge10_pass += sum(e - s for s, e in merge_intervals(ivs_10_pass))

    pct10 = round(100.0 * bp_ge10_pass / total_bp, 2) if total_bp else 0.0
    pct1 = round(100.0 * bp_ge1 / total_bp, 2) if total_bp else 0.0
    return {
        "cds_bp": total_bp,
        "pct_ge_10x": pct10,
        "pct_ge_1x": pct1,
        "pct_0x_nocall": round(max(0.0, 100.0 - pct1), 2),
    }


def get_reference_metadata(reference_dir: Path) -> dict[str, str]:
    """Read reference provenance (ClinVar fileDate, GenCC mtime, ACMG SF version, run date)."""
    meta = {
        "run_date": datetime.date.today().isoformat(),
        "clinvar_date": "unknown",
        "gencc_date": "unknown",
        "gencode_release": "v46",
        "acmg_sf_version": f"{ACMG_SF_VERSION} ({len(ACMG_SF_GENES)} genes)",
    }
    clinvar_vcf = reference_dir / "clinvar_GRCh38.norm.vcf.gz"
    if not clinvar_vcf.exists():
        clinvar_vcf = reference_dir / "clinvar_GRCh38.vcf.gz"
    if clinvar_vcf.exists():
        try:
            with pysam.VariantFile(str(clinvar_vcf)) as vf:
                for rec in vf.header.records:
                    if rec.key == "fileDate":
                        meta["clinvar_date"] = str(rec.value)
                        break
        except Exception:
            pass
        if meta["clinvar_date"] == "unknown":
            ts = clinvar_vcf.stat().st_mtime
            meta["clinvar_date"] = datetime.date.fromtimestamp(ts).isoformat()

    gencc_tsv = reference_dir / "gencc_submissions.tsv"
    if gencc_tsv.exists():
        ts = gencc_tsv.stat().st_mtime
        meta["gencc_date"] = datetime.date.fromtimestamp(ts).isoformat()

    return meta
