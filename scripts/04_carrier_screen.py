#!/usr/bin/env python3
"""
Step 4: Clinical Pathology & Carrier Screening Pipeline (Offline / Local)
Intersects normalized personal SNP/indel, CNV, and SV VCFs against normalized
NCBI ClinVar GRCh38, maps Mode of Inheritance via GenCC + ACMG SF v3.3 (84 genes),
computes call QC tiers, and stratifies findings into clinical buckets.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.constants import ACMG_SF_GENES, ACMG_SF_VERSION, PRIMARY_CHROMS_SET
from lib.genomics import (
    classify_clnsig,
    extract_genes,
    get_reference_metadata,
    load_gencc_moi,
    parse_info_str,
    parse_revstat_stars,
    resolve_gene_moi,
)

REPO_ROOT = SCRIPT_DIR.parent
REFERENCE_DIR = REPO_ROOT / "data" / "reference"

_CATEGORY_PRIORITY = {"PLP": 3, "CONFLICTING_PLP": 2, "RISK_OR_DRUG": 1}


def evaluate_personal_call(rec: pysam.VariantRecord) -> dict[str, Any]:
    """Compute zygosity, depth/GQ/allele-balance metrics, and transparent QC tier."""
    sample = next(iter(rec.samples.values()))
    gt = sample.get("GT")
    gt_str = "/".join("." if a is None else str(a) for a in gt) if gt else "./."
    non_ref_count = sum(1 for a in (gt or ()) if a is not None and a > 0)
    ploidy = len(gt) if gt else 0

    if ploidy == 1 and non_ref_count == 1:
        zygosity = "HEMI_ALT"
    elif ploidy == 2 and non_ref_count == 2:
        zygosity = "HOM_ALT"
    elif ploidy == 2 and non_ref_count == 1:
        zygosity = "HET"
    else:
        zygosity = f"OTHER({gt_str})"

    flt = ";".join(rec.filter.keys()) if len(rec.filter) > 0 else "PASS"
    dp = sample.get("DP")
    if dp is None:
        dp = rec.info.get("DP")
    gq = sample.get("GQ")
    ad = sample.get("AD")
    ad_ref = ad[0] if ad and len(ad) >= 1 else None
    ad_alt = ad[1] if ad and len(ad) >= 2 else None
    ab = None
    if ad_ref is not None and ad_alt is not None and (ad_ref + ad_alt) > 0:
        ab = round(ad_alt / (ad_ref + ad_alt), 3)
    mq = rec.info.get("MQ")

    reasons: list[str] = []
    if flt != "PASS":
        reasons.append(f"FILTER={flt}")
    min_dp = 6 if zygosity == "HEMI_ALT" else 10
    if dp is None or dp < min_dp:
        reasons.append(f"DP={dp}<{min_dp}")
    if gq is None or gq < 20:
        reasons.append(f"GQ={gq}<20")
    if ab is not None:
        if zygosity == "HET" and (ab < 0.25 or ab > 0.75):
            reasons.append(f"AB={ab}")
        elif zygosity in ("HOM_ALT", "HEMI_ALT") and ab < 0.85:
            reasons.append(f"AB={ab}")

    if flt != "PASS" or (ab is not None and ab < 0.20) or (dp is not None and dp < 6):
        qc_tier = "LOW_CONFIDENCE"
    elif reasons:
        qc_tier = "MODERATE"
    else:
        qc_tier = "HIGH"

    return {
        "gt": gt_str,
        "zygosity": zygosity,
        "filter": flt,
        "qual": round(float(rec.qual), 1) if rec.qual is not None else None,
        "dp": dp,
        "gq": gq,
        "ad_ref": ad_ref,
        "ad_alt": ad_alt,
        "ab": ab,
        "mq": mq,
        "qc_tier": qc_tier,
        "qc_notes": "; ".join(reasons) if reasons else "OK",
    }


def load_clinvar_index(
    clinvar_vcf: Path,
    custom_vcf: Path | None = None,
) -> tuple[
    dict[tuple[str, int, str, str], dict[str, Any]],
    dict[str, list[tuple[int, int, str, int]]],
    list[str],
]:
    """
    Single-pass loader for ClinVar (and optional custom clinical VCF) that builds:
      1. clin_index: {(chrom, pos, ref, alt): clinical metadata}
      2. gene_intervals: per-chrom P/LP (>=1*) gene bounds for CNV/SV overlap screening
      3. clinvar_plp_genes: sorted list of >=1* P/LP gene symbols (reused by Step 7)
    """
    print(f"[ClinVar] Single-pass indexing from {clinvar_vcf.name}...")
    clin_index: dict[tuple[str, int, str, str], dict[str, Any]] = {}
    clin_counts: dict[str, int] = defaultdict(int)
    gene_bounds: dict[tuple[str, str], list[int]] = {}
    gene_plp_counts: dict[tuple[str, str], int] = defaultdict(int)
    plp_genes_set: set[str] = set()

    def _ingest_vcf(vcf_path: Path, is_custom: bool = False) -> None:
        with pysam.VariantFile(str(vcf_path)) as cvf:
            for rec in cvf:
                if rec.chrom not in PRIMARY_CHROMS_SET or not rec.alts:
                    continue
                clnsig = parse_info_str(rec.info.get("CLNSIG"))
                clnsigconf = parse_info_str(rec.info.get("CLNSIGCONF"))
                if is_custom and not clnsig:
                    clnsig = "Pathogenic"
                cat = classify_clnsig(clnsig, clnsigconf)
                if cat is None:
                    continue
                clin_counts[cat] += 1
                revstat, stars = parse_revstat_stars(rec.info.get("CLNREVSTAT"))
                if is_custom and stars == 0 and not revstat:
                    revstat = "custom_clinical_panel"
                    stars = 2
                geneinfo = parse_info_str(rec.info.get("GENEINFO"))
                genes = extract_genes(geneinfo)
                clndn = parse_info_str(rec.info.get("CLNDN"))
                mc = parse_info_str(rec.info.get("MC"))
                rs = parse_info_str(rec.info.get("RS"))

                if cat == "PLP" and stars >= 1:
                    end = rec.stop or rec.pos
                    for g in genes:
                        plp_genes_set.add(g)
                        k = (rec.chrom, g)
                        if k not in gene_bounds:
                            gene_bounds[k] = [rec.pos, end]
                        else:
                            if rec.pos < gene_bounds[k][0]:
                                gene_bounds[k][0] = rec.pos
                            if end > gene_bounds[k][1]:
                                gene_bounds[k][1] = end
                        gene_plp_counts[k] += 1

                for alt in rec.alts:
                    key = (rec.chrom, rec.pos, rec.ref, alt)
                    candidate = {
                        "clinvar_id": rec.id,
                        "clinvar_category": cat,
                        "clnsig": clnsig,
                        "clnsigconf": clnsigconf,
                        "clnrevstat": revstat,
                        "clinvar_stars": stars,
                        "genes": genes,
                        "gene_primary": genes[0] if genes else "",
                        "clndn": clndn,
                        "mc": mc,
                        "rsid_clinvar": f"rs{rs}" if rs else "",
                    }
                    existing = clin_index.get(key)
                    if existing is not None:
                        old_rank = (_CATEGORY_PRIORITY.get(existing["clinvar_category"], 0), existing["clinvar_stars"])
                        new_rank = (_CATEGORY_PRIORITY.get(cat, 0), stars)
                        if new_rank <= old_rank:
                            continue
                    clin_index[key] = candidate

    _ingest_vcf(clinvar_vcf, is_custom=False)
    if custom_vcf is not None and custom_vcf.exists():
        print(f"[Custom] Merging optional custom clinical panel: {custom_vcf.name}...")
        _ingest_vcf(custom_vcf, is_custom=True)

    gene_intervals: dict[str, list[tuple[int, int, str, int]]] = defaultdict(list)
    for (chrom, gene), (start, end) in gene_bounds.items():
        gene_intervals[chrom].append((start, end, gene, gene_plp_counts[(chrom, gene)]))
    for chrom in gene_intervals:
        gene_intervals[chrom].sort()

    print(
        f"[ClinVar] Indexed {len(clin_index):,} clinical alleles "
        f"(PLP={clin_counts['PLP']:,}, CONFLICTING_PLP={clin_counts['CONFLICTING_PLP']:,}, "
        f"RISK_OR_DRUG={clin_counts['RISK_OR_DRUG']:,}; {len(plp_genes_set):,} P/LP >=1* genes)"
    )
    return clin_index, gene_intervals, sorted(plp_genes_set)


def _bcftools_record_count(vcf_path: Path) -> int | None:
    if not shutil.which("bcftools"):
        return None
    try:
        out = subprocess.check_output(["bcftools", "index", "-n", str(vcf_path)], text=True).strip()
        return int(out) if out.isdigit() else None
    except Exception:
        return None


def screen_snp_indels(
    personal_vcf: Path,
    clin_index: dict[tuple[str, int, str, str], dict[str, Any]],
    gencc_moi: dict[str, dict[str, Any]],
    targets_vcf: Path | None = None,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """
    Intersect normalized personal SNP/indel VCF against pre-loaded ClinVar index.
    When `bcftools` and `targets_vcf` are available, uses C-level `bcftools view -T`
    to pre-filter personal calls before Python inspection.
    """
    hits: list[dict[str, Any]] = []
    personal_stats: dict[str, int] = defaultdict(int)

    scan_vcf = personal_vcf
    tmp_path: Path | None = None
    used_prefilter = False

    if targets_vcf is not None and targets_vcf.exists() and shutil.which("bcftools"):
        total_records = _bcftools_record_count(personal_vcf)
        if total_records is not None:
            tmp_file = tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False)
            tmp_file.close()
            tmp_path = Path(tmp_file.name)
            try:
                subprocess.run(
                    [
                        "bcftools", "view", "--threads", "2",
                        "-T", str(targets_vcf),
                        "-Oz", "-o", str(tmp_path),
                        str(personal_vcf),
                    ],
                    check=True,
                    stderr=subprocess.DEVNULL,
                )
                scan_vcf = tmp_path
                personal_stats["total_variant_records"] = total_records
                personal_stats["non_ref_calls"] = total_records
                used_prefilter = True
            except Exception:
                scan_vcf = personal_vcf
                if tmp_path.exists():
                    tmp_path.unlink(missing_ok=True)
                tmp_path = None

    try:
        with pysam.VariantFile(str(scan_vcf)) as pvf:
            for rec in pvf:
                if not used_prefilter:
                    personal_stats["total_variant_records"] += 1
                if not rec.alts:
                    continue
                sample = next(iter(rec.samples.values()))
                gt = sample.get("GT")
                if not gt or all(a is None or a == 0 for a in gt):
                    continue
                if not used_prefilter:
                    personal_stats["non_ref_calls"] += 1

                present_alt_idxs = {a for a in gt if a is not None and a > 0}
                for idx in present_alt_idxs:
                    alt = rec.alts[idx - 1]
                    key = (rec.chrom, rec.pos, rec.ref, alt)
                    clin = clin_index.get(key)
                    if clin is None:
                        continue

                    qc = evaluate_personal_call(rec)
                    primary_gene, moi, gencc_diseases = resolve_gene_moi(
                        clin["genes"], gencc_moi, rec.chrom
                    )
                    is_acmg = any(g in ACMG_SF_GENES for g in clin["genes"])

                    hits.append({
                        "chrom": rec.chrom,
                        "pos": rec.pos,
                        "ref": rec.ref,
                        "alt": alt,
                        "variant_id": f"{rec.chrom}-{rec.pos}-{rec.ref}-{alt}",
                        "rsid": rec.id or clin["rsid_clinvar"] or ".",
                        "gene": primary_gene,
                        "all_genes": ",".join(clin["genes"]),
                        "acmg_sf": is_acmg,
                        "acmg_sf_version": ACMG_SF_VERSION,
                        "moi": moi,
                        "gencc_diseases": gencc_diseases,
                        **clin,
                        **qc,
                    })
    finally:
        if tmp_path is not None and tmp_path.exists():
            tmp_path.unlink(missing_ok=True)

    mode_str = "bcftools-prefiltered" if used_prefilter else "full-scan"
    print(
        f"[SNP/Indel] Scanned {personal_stats['total_variant_records']:,} personal variant records "
        f"({mode_str}) -> {len(hits):,} ClinVar matches."
    )
    return hits, dict(personal_stats)


def stratify_hits(hits: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """Stratify SNP/indel hits into clinical buckets."""
    plp_hits = [h for h in hits if h["clinvar_category"] == "PLP"]
    conf_hits = [h for h in hits if h["clinvar_category"] == "CONFLICTING_PLP"]
    risk_drug_hits = [h for h in hits if h["clinvar_category"] == "RISK_OR_DRUG"]

    het_by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for h in plp_hits:
        if h["zygosity"] == "HET" and h["gene"]:
            het_by_gene[h["gene"]].append(h)
    comphet_genes = {g for g, lst in het_by_gene.items() if len(lst) >= 2}

    buckets: dict[str, list[dict[str, Any]]] = {
        "incidental_ad_acmg": [],
        "hom_hemi_comphet": [],
        "carrier_ar_xl": [],
        "plp_unknown_moi": [],
        "conflicting_plp": conf_hits,
        "risk_or_drug": risk_drug_hits,
    }

    for h in plp_hits:
        h["potential_comphet"] = h["gene"] in comphet_genes and h["zygosity"] == "HET"
        if h["acmg_sf"] or h["moi"] in ("AD", "AD/AR", "MT"):
            buckets["incidental_ad_acmg"].append(h)

        if h["zygosity"] in ("HOM_ALT", "HEMI_ALT") or h["potential_comphet"]:
            buckets["hom_hemi_comphet"].append(h)
        elif h["moi"] in ("AR", "XL") and h["zygosity"] == "HET":
            buckets["carrier_ar_xl"].append(h)
        elif not (h["acmg_sf"] or h["moi"] in ("AD", "AD/AR", "MT")):
            buckets["plp_unknown_moi"].append(h)

    tier_rank = {"HIGH": 0, "MODERATE": 1, "LOW_CONFIDENCE": 2}
    for k in buckets:
        buckets[k].sort(
            key=lambda x: (
                -x["clinvar_stars"],
                tier_rank.get(x["qc_tier"], 9),
                x["gene"],
                x["chrom"],
                x["pos"],
            )
        )
    return buckets


def screen_cnv_sv(
    cnv_vcf: Path,
    sv_vcf: Path,
    gene_intervals: dict[str, list[tuple[int, int, str, int]]],
    gencc_moi: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    """Screen non-REF Canvas CNVs and Manta SVs against ClinVar P/LP gene intervals."""
    overlaps: list[dict[str, Any]] = []

    def find_overlapping_genes(chrom: str, start: int, end: int) -> list[tuple[str, int]]:
        hits = []
        for g_start, g_end, g_sym, n_plp in gene_intervals.get(chrom, []):
            if g_end < start:
                continue
            if g_start > end:
                break
            hits.append((g_sym, n_plp))
        return hits

    cnv_intervals_pass: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with pysam.VariantFile(str(cnv_vcf)) as vf:
        for rec in vf:
            if not rec.alts or rec.alts == (".",) or "Canvas:REF" in (rec.id or ""):
                continue
            flt = ";".join(rec.filter.keys()) if len(rec.filter) > 0 else "PASS"
            svtype = parse_info_str(rec.info.get("SVTYPE"))
            end = rec.stop or rec.info.get("END") or rec.pos
            length = abs(end - rec.pos)
            sample = next(iter(rec.samples.values()))
            cn = sample.get("CN")
            qs = sample.get("QS")
            ft = sample.get("FT")
            if flt == "PASS" and svtype != "LOH":
                cnv_intervals_pass[rec.chrom].append((rec.pos, end))

            genes = find_overlapping_genes(rec.chrom, rec.pos, end)
            if not genes:
                continue
            gene_syms = [g for g, _ in genes]
            acmg_hits = [g for g in gene_syms if g in ACMG_SF_GENES]
            if flt != "PASS" and not acmg_hits:
                continue
            overlaps.append({
                "source": "Canvas_CNV",
                "chrom": rec.chrom,
                "start": rec.pos,
                "end": end,
                "length_bp": length,
                "svtype": svtype or ("LOSS" if "LOSS" in (rec.id or "") else "GAIN"),
                "call_id": rec.id,
                "filter": flt,
                "sample_ft": str(ft),
                "cn_or_gt": f"CN={cn}",
                "qual_or_qs": round(float(qs), 1) if qs is not None else rec.qual,
                "n_clinvar_genes": len(gene_syms),
                "acmg_genes": ",".join(acmg_hits),
                "genes_preview": ",".join(gene_syms[:12]) + ("..." if len(gene_syms) > 12 else ""),
                "moi_summary": ",".join(sorted({gencc_moi.get(g, {}).get("moi", "UNKNOWN") for g in gene_syms[:12]})),
                "artifact_flag": "LOH_BLOCK" if svtype == "LOH" else "OK",
            })

    with pysam.VariantFile(str(sv_vcf)) as vf:
        for rec in vf:
            if rec.chrom not in PRIMARY_CHROMS_SET:
                continue
            flt = ";".join(rec.filter.keys()) if len(rec.filter) > 0 else "PASS"
            if flt != "PASS":
                continue
            svtype = parse_info_str(rec.info.get("SVTYPE"))
            if svtype == "BND":
                alt_str = ",".join(rec.alts or ())
                if "decoy" in alt_str or "random" in alt_str or "Un_" in alt_str:
                    continue
                end = rec.pos + 1
                length = 1
            else:
                end = rec.stop or rec.info.get("END") or rec.pos
                svlen_raw = rec.info.get("SVLEN")
                if isinstance(svlen_raw, (tuple, list)) and svlen_raw:
                    length = abs(int(svlen_raw[0]))
                elif svlen_raw is not None:
                    length = abs(int(svlen_raw))
                else:
                    length = abs(end - rec.pos)

            artifact_flag = "OK"
            if svtype in ("DEL", "DUP") and length > 5_000_000:
                has_cnv_support = any(
                    not (c_end < rec.pos or c_start > end)
                    for c_start, c_end in cnv_intervals_pass.get(rec.chrom, [])
                )
                if not has_cnv_support:
                    artifact_flag = "LIKELY_ARTIFACT_MEGA_SV_NO_CNV_SUPPORT"

            genes = find_overlapping_genes(rec.chrom, rec.pos, end)
            if not genes:
                continue
            gene_syms = [g for g, _ in genes]
            acmg_hits = [g for g in gene_syms if g in ACMG_SF_GENES]
            sample = next(iter(rec.samples.values()))
            gt = sample.get("GT")
            gt_str = "/".join("." if a is None else str(a) for a in gt) if gt else "./."
            pr = sample.get("PR")
            sr = sample.get("SR")

            overlaps.append({
                "source": "Manta_SV",
                "chrom": rec.chrom,
                "start": rec.pos,
                "end": end,
                "length_bp": length,
                "svtype": svtype,
                "call_id": rec.id,
                "filter": flt,
                "sample_ft": f"PR={pr};SR={sr}",
                "cn_or_gt": gt_str,
                "qual_or_qs": round(float(rec.qual), 1) if rec.qual is not None else None,
                "n_clinvar_genes": len(gene_syms),
                "acmg_genes": ",".join(acmg_hits),
                "genes_preview": ",".join(gene_syms[:12]) + ("..." if len(gene_syms) > 12 else ""),
                "moi_summary": ",".join(sorted({gencc_moi.get(g, {}).get("moi", "UNKNOWN") for g in gene_syms[:12]})),
                "artifact_flag": artifact_flag,
            })

    overlaps.sort(
        key=lambda x: (
            0 if x["artifact_flag"] == "OK" else 1,
            0 if x["acmg_genes"] else 1,
            x["source"],
            x["chrom"],
            x["start"],
        )
    )
    return overlaps


def main() -> None:
    parser = argparse.ArgumentParser(description="Run per-sample clinical carrier screen.")
    parser.add_argument("--sample", required=True, help="Sample label from config/samples.tsv (e.g. partner_a, partner_b)")
    args = parser.parse_args()

    processed_dir = REPO_ROOT / "data" / "processed" / args.sample
    reports_dir = REPO_ROOT / "results" / "reports" / args.sample
    qc_dir = REPO_ROOT / "results" / "qc" / args.sample
    reports_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    personal_vcf = processed_dir / "snp_indel.variants.norm.vcf.gz"
    clinvar_clinical_vcf = REFERENCE_DIR / "clinvar_clinical.norm.vcf.gz"
    clinvar_full_vcf = REFERENCE_DIR / "clinvar_GRCh38.norm.vcf.gz"
    clinvar_vcf = clinvar_clinical_vcf if clinvar_clinical_vcf.exists() else clinvar_full_vcf
    custom_vcf = REFERENCE_DIR / "custom_clinical.norm.vcf.gz"
    cnv_vcf = processed_dir / "cnv.vcf.gz"
    sv_vcf = processed_dir / "sv.vcf.gz"
    gencc_tsv = REFERENCE_DIR / "gencc_submissions.tsv"

    ref_meta = get_reference_metadata(REFERENCE_DIR)
    print(f"[Sample: {args.sample}] ClinVar={ref_meta['clinvar_date']}, GenCC={ref_meta['gencc_date']}, ACMG={ref_meta['acmg_sf_version']}")
    gencc_moi = load_gencc_moi(gencc_tsv)
    print(f"[MOI] Loaded curated MOI for {len(gencc_moi):,} genes.")

    clin_index, gene_intervals, clinvar_plp_genes = load_clinvar_index(
        clinvar_vcf, custom_vcf if custom_vcf.exists() else None
    )
    # Use fast bcftools C-level target pre-filtering when no separate custom VCF is active
    targets_vcf = clinvar_vcf if not custom_vcf.exists() else None
    hits, stats = screen_snp_indels(personal_vcf, clin_index, gencc_moi, targets_vcf=targets_vcf)
    buckets = stratify_hits(hits)

    cnv_sv_overlaps = screen_cnv_sv(cnv_vcf, sv_vcf, gene_intervals, gencc_moi)

    for name, rows in buckets.items():
        out_tsv = reports_dir / f"snp_indel_{name}.tsv"
        pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)

    pd.DataFrame(cnv_sv_overlaps).to_csv(
        reports_dir / "cnv_sv_clinical_overlaps.tsv", sep="\t", index=False
    )

    summary = {
        "sample": args.sample,
        "reference_metadata": ref_meta,
        "personal_stats": stats,
        "clinvar_plp_genes": clinvar_plp_genes,
        "counts": {k: len(v) for k, v in buckets.items()},
        "cnv_sv_total": len(cnv_sv_overlaps),
        "cnv_sv_non_artifact": sum(1 for x in cnv_sv_overlaps if x["artifact_flag"] == "OK"),
        "buckets": buckets,
        "cnv_sv_overlaps": cnv_sv_overlaps,
    }
    with open(reports_dir / "carrier_screen_raw.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"\n================ STEP 4 SCREENING SUMMARY ({args.sample}) ================")
    for k, v in buckets.items():
        high = sum(1 for x in v if x["qc_tier"] == "HIGH")
        print(f"  {k:22s}: {len(v):4d} total ({high:4d} HIGH QC)")
    print(
        f"  {'cnv_sv_overlaps':22s}: {len(cnv_sv_overlaps):4d} total "
        f"({summary['cnv_sv_non_artifact']:4d} non-artifact)"
    )
    print(f"Outputs written to {reports_dir}/")


if __name__ == "__main__":
    main()
