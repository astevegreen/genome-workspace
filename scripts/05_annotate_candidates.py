#!/usr/bin/env python3
"""
Step 5: Enrich Candidate Hits with Cached gnomAD v4 + Ensembl Annotations,
Flag Likely Artifacts/Common Polymorphisms (BA1/BS1), Scan ACMG SF v3.3 Genes
for Novel/Unreported High-Impact LoF Variants, and Render the Carrier Report.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.api_client import CachedAnnotationClient, format_af
from lib.constants import ACMG_SF_GENES, ACMG_SF_VERSION, NOVEL_LOF_MAX_AF
from lib.genomics import (
    assess_artifact_or_common,
    get_reference_metadata,
    load_gencode_cds,
    parse_info_str,
    parse_revstat_stars,
)

REPO_ROOT = SCRIPT_DIR.parent
REFERENCE_DIR = REPO_ROOT / "data" / "reference"


def _is_clinvar_benign(
    cvf: pysam.VariantFile | None,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> bool:
    """Return True if exact allele is ClinVar Benign / Likely_benign (>=1 star)."""
    if cvf is None:
        return False
    try:
        for rec in cvf.fetch(chrom, max(0, pos - 1), pos + len(ref)):
            if rec.pos == pos and rec.ref == ref and rec.alts and alt in rec.alts:
                clnsig = parse_info_str(rec.info.get("CLNSIG")).lower()
                _, stars = parse_revstat_stars(rec.info.get("CLNREVSTAT"))
                if stars >= 1 and "conflicting" not in clnsig and "pathogenic" not in clnsig:
                    tokens = [t.strip() for t in clnsig.replace("|", "/").replace(",", "/").split("/")]
                    if any(t in ("benign", "likely_benign") for t in tokens):
                        return True
    except ValueError:
        pass
    return False


def scan_novel_acmg_lof(
    personal_vcf: Path,
    cds_bed_gz: Path,
    known_plp_vids: set[str],
    api: CachedAnnotationClient,
    clinvar_vcf: Path | None = None,
) -> list[dict[str, Any]]:
    """
    Scan all 84 ACMG SF v3.3 genes (GENCODE v46 CDS ±2 bp canonical splice sites)
    for unreported (non-ClinVar P/LP and non-ClinVar Benign/LB) HIGH-impact
    loss-of-function variants with passing call QC and gnomAD v4 AF < 0.1% (or absent).
    """
    if not cds_bed_gz.exists() or not personal_vcf.exists():
        return []

    _, cds_by_gene = load_gencode_cds(cds_bed_gz)
    candidates: list[dict[str, Any]] = []
    vep_batch_inputs: list[tuple[str, str, str, int, str, str]] = []
    seen_vids: set[str] = set()

    cvf = pysam.VariantFile(str(clinvar_vcf)) if (clinvar_vcf and clinvar_vcf.exists()) else None
    try:
        with pysam.VariantFile(str(personal_vcf)) as vf:
            for gene in sorted(ACMG_SF_GENES):
                entry = cds_by_gene.get(gene)
                if not entry:
                    continue
                chrom, exons = entry
                for ex_s, ex_e in exons:
                    fetch_s = max(0, ex_s - 2)
                    fetch_e = ex_e + 2
                    for rec in vf.fetch(chrom, fetch_s, fetch_e):
                        if not rec.alts:
                            continue
                        sample = next(iter(rec.samples.values()))
                        gt = sample.get("GT")
                        if not gt or all(a is None or a == 0 for a in gt):
                            continue
                        flt = ";".join(rec.filter.keys()) if len(rec.filter) > 0 else "PASS"
                        if flt != "PASS":
                            continue

                        dp = sample.get("DP")
                        if dp is None:
                            dp = rec.info.get("DP")
                        gq = sample.get("GQ")
                        ploidy = len(gt)
                        min_dp = 6 if ploidy == 1 else 10
                        if dp is None or dp < min_dp or (gq is not None and gq < 20):
                            continue

                        ad = sample.get("AD")
                        present_alt_idxs = {a for a in gt if a is not None and a > 0}
                        for idx in present_alt_idxs:
                            alt = rec.alts[idx - 1]
                            vid = f"{rec.chrom}-{rec.pos}-{rec.ref}-{alt}"
                            if vid in known_plp_vids or vid in seen_vids:
                                continue
                            ab = (
                                round(ad[idx] / (ad[0] + ad[idx]), 3)
                                if ad and len(ad) > idx and (ad[0] + ad[idx]) > 0
                                else None
                            )
                            if ab is not None and ab < 0.20:
                                continue
                            if _is_clinvar_benign(cvf, rec.chrom, rec.pos, rec.ref, alt):
                                continue
                            seen_vids.add(vid)
                            gt_str = "/".join("." if a is None else str(a) for a in gt)
                            candidates.append({
                                "gene": gene,
                                "variant_id": vid,
                                "chrom": rec.chrom,
                                "pos": rec.pos,
                                "ref": rec.ref,
                                "alt": alt,
                                "rsid": rec.id or ".",
                                "gt": gt_str,
                                "dp": dp,
                                "gq": gq,
                                "ab": ab,
                                "filter": flt,
                            })
                            vep_batch_inputs.append((vid, gene, rec.chrom, rec.pos, rec.ref, alt))
    finally:
        if cvf is not None:
            cvf.close()

    if not candidates:
        return []

    print(f"[Novel ACMG LoF] Running VEP on {len(candidates)} coding/splice calls across {len(ACMG_SF_GENES)} ACMG SF {ACMG_SF_VERSION} genes...")
    vep_map = api.query_ensembl_vep_batch(vep_batch_inputs)

    novel_lofs: list[dict[str, Any]] = []
    for c in candidates:
        vinfo = vep_map.get(c["variant_id"], {})
        if vinfo.get("impact") != "HIGH":
            continue
        g_info = api.query_gnomad_variant(c["variant_id"])
        af = g_info.get("gnomad_af")
        if af is None or af < NOVEL_LOF_MAX_AF:
            novel_lofs.append({**c, **vinfo, **g_info})

    print(f"[Novel ACMG LoF] Identified {len(novel_lofs)} unreported rare HIGH-impact LoF candidate(s).")
    return novel_lofs


def render_variant_table(rows: list[dict[str, Any]], include_filter_reason: bool = False) -> str:
    if not rows:
        return "_None identified._\n"
    header = "| Gene | Variant (GRCh38) | rsID | GT / Zyg | QC (DP/GQ/AB) | ClinVar Sig (Stars) | MOI | gnomAD v4 AF (Popmax) | Condition(s) |"
    sep = "|---|---|---|---|---|---|---|---|---|"
    if include_filter_reason:
        header += " Filter Reason |"
        sep += "---|"
    lines = [header, sep]

    for r in rows:
        stars = "★" * r["clinvar_stars"] + "☆" * (4 - r["clinvar_stars"])
        acmg = " **[ACMG]**" if r.get("acmg_sf") else ""
        comp = " **[CompHet?]**" if r.get("potential_comphet") else ""
        qc_str = f"**{r['qc_tier']}** ({r['dp']}/{r['gq']}/{r['ab']})"
        af_str = format_af(r.get("gnomad_af"), r.get("gnomad_status"))
        if r.get("gnomad_popmax_af") is not None:
            af_str += f" ({r['gnomad_popmax_pop']}: {format_af(r['gnomad_popmax_af'], 'ok')})"
        cond = (r.get("clndn") or "").replace("|", "; ").replace("_", " ")
        if len(cond) > 80:
            cond = cond[:77] + "..."
        pos_str = f"`{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}`"
        if len(pos_str) > 36:
            pos_str = f"`{r['chrom']}:{r['pos']} {r['ref'][:8]}..>{r['alt'][:8]}..`"
        row_str = (
            f"| **{r['gene'] or '.'}**{acmg}{comp} | {pos_str} | {r['rsid']} | "
            f"`{r['gt']}` ({r['zygosity']}) | {qc_str} | {r['clnsig']} ({stars}) | "
            f"{r['moi']} | {af_str} | {cond} |"
        )
        if include_filter_reason:
            row_str += f" {r.get('artifact_reason', '')} |"
        lines.append(row_str)
    return "\n".join(lines) + "\n"


def render_cnv_sv_table(rows: list[dict[str, Any]], max_rows: int = 40) -> str:
    if not rows:
        return "_None identified._\n"
    lines = [
        "| Source | Type | Region (GRCh38) | Size | Call / GT | Filter / Support | ACMG Genes | ClinVar P/LP Genes Overlapped | Notes |",
        "|---|---|---|---:|---|---|---|---|---|",
    ]
    for r in rows[:max_rows]:
        size_kb = f"{r['length_bp'] / 1000:,.1f} kb" if r["length_bp"] >= 1000 else f"{r['length_bp']} bp"
        lines.append(
            f"| {r['source']} | {r['svtype']} | `{r['chrom']}:{r['start']:,}-{r['end']:,}` | "
            f"{size_kb} | `{r['cn_or_gt']}` | {r['filter']} ({r['sample_ft']}) | "
            f"**{r['acmg_genes'] or '-'}** | {r['genes_preview']} ({r['n_clinvar_genes']}) | {r['artifact_flag']} |"
        )
    if len(rows) > max_rows:
        lines.append(f"\n_Showing top {max_rows} of {len(rows)} entries; see `cnv_sv_clinical_overlaps.tsv` for full table._")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Annotate candidate P/LP & conflicting hits via cached gnomAD v4 & Ensembl.")
    parser.add_argument("--sample", required=True, help="Sample label from config/samples.tsv (e.g. partner_a, partner_b)")
    parser.add_argument("--refresh-cache", action="store_true", help="Ignore existing api_cache.json and re-query APIs")
    args = parser.parse_args()

    processed_dir = REPO_ROOT / "data" / "processed" / args.sample
    reports_dir = REPO_ROOT / "results" / "reports" / args.sample
    raw_path = reports_dir / "carrier_screen_raw.json"
    with open(raw_path, encoding="utf-8") as f:
        data = json.load(f)

    ref_meta = data.get("reference_metadata") or get_reference_metadata(REFERENCE_DIR)
    data["reference_metadata"] = ref_meta

    buckets = data["buckets"]
    primary_bucket_names = [
        "incidental_ad_acmg",
        "hom_hemi_comphet",
        "carrier_ar_xl",
        "plp_unknown_moi",
    ]
    annotate_bucket_names = primary_bucket_names + ["conflicting_plp"]

    unique_variants: set[str] = set()
    unique_genes: set[str] = set()
    for bname in annotate_bucket_names:
        for row in buckets.get(bname, []):
            unique_variants.add(row["variant_id"])
            if row["gene"]:
                unique_genes.add(row["gene"])

    api = CachedAnnotationClient(
        REFERENCE_DIR / "api_cache.json",
        refresh=True if args.refresh_cache else None,
    )
    print(
        f"[Annotate: {args.sample}] Enriching {len(unique_variants)} unique variants "
        f"(P/LP + Conflicting P/LP) across {len(unique_genes)} genes (with local cache)..."
    )

    var_cache: dict[str, dict[str, Any]] = {}
    for i, vid in enumerate(sorted(unique_variants), 1):
        var_cache[vid] = api.query_gnomad_variant(vid)
        if i % 15 == 0 or i == len(unique_variants):
            print(f"  variants: {i}/{len(unique_variants)}")

    gene_cache: dict[str, dict[str, Any]] = {}
    for i, g in enumerate(sorted(unique_genes), 1):
        gene_cache[g] = api.query_gnomad_gene(g)
        if i % 15 == 0 or i == len(unique_genes):
            print(f"  genes:    {i}/{len(unique_genes)}")

    # Merge annotations & apply automated BA1/BS1 artifact flagging
    enriched_rows: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    filtered_artifacts: list[dict[str, Any]] = []
    filtered_seen: set[str] = set()

    for bname in annotate_bucket_names:
        for row in buckets.get(bname, []):
            row.update(var_cache.get(row["variant_id"], {}))
            row.update(gene_cache.get(row["gene"], {}))
            flag, reason = assess_artifact_or_common(row)
            row["artifact_flag"] = flag
            row["artifact_reason"] = reason
            if bname in primary_bucket_names:
                if row["variant_id"] not in seen_ids:
                    seen_ids.add(row["variant_id"])
                    enriched_rows.append({"primary_bucket": bname, **row})
                if flag != "OK" and row["variant_id"] not in filtered_seen:
                    filtered_seen.add(row["variant_id"])
                    filtered_artifacts.append(row)

    # Scan for unreported HIGH-impact LoF variants across all 84 ACMG SF v3.3 genes
    known_plp_vids = {r["variant_id"] for r in enriched_rows}
    novel_acmg_lofs = scan_novel_acmg_lof(
        processed_dir / "snp_indel.variants.norm.vcf.gz",
        REFERENCE_DIR / "gencode.v46.cds.nochr.bed.gz",
        known_plp_vids,
        api,
        clinvar_vcf=REFERENCE_DIR / "clinvar_GRCh38.norm.vcf.gz",
    )
    data["novel_acmg_lof"] = novel_acmg_lofs
    api.save()

    pd.DataFrame(enriched_rows).to_csv(
        reports_dir / "carrier_screening_plp_enriched.tsv", sep="\t", index=False
    )

    # Split primary buckets into clean vs. filtered-artifact views for the Markdown report
    clean_incidental = [r for r in buckets["incidental_ad_acmg"] if r["artifact_flag"] == "OK"]
    clean_hom_hemi = [r for r in buckets["hom_hemi_comphet"] if r["artifact_flag"] == "OK"]
    clean_carrier = [r for r in buckets["carrier_ar_xl"] if r["artifact_flag"] == "OK"]
    clean_unknown = [r for r in buckets["plp_unknown_moi"] if r["artifact_flag"] == "OK"]

    cnv_sv_ok = [x for x in data["cnv_sv_overlaps"] if x["artifact_flag"] == "OK"]
    cnv_sv_flagged = [x for x in data["cnv_sv_overlaps"] if x["artifact_flag"] != "OK"]

    stats = data["personal_stats"]
    counts = data["counts"]
    md = [
        f"# Personal 30x WGS Carrier Screening & Clinical Variant Report ({args.sample.title()})",
        "",
        f"**Run Date:** {ref_meta['run_date']} | **ClinVar:** {ref_meta['clinvar_date']} | "
        f"**GenCC:** {ref_meta['gencc_date']} | **ACMG SF:** {ref_meta['acmg_sf_version']} | **GENCODE:** {ref_meta['gencode_release']}",
        "",
        "> **Disclaimer:** This is a research-grade bioinformatics analysis for informational and educational exploration only. "
        "It is **not** a clinical diagnostic test. Any actionable, incidental, or reproductive carrier finding must be confirmed "
        "by a CLIA-certified clinical laboratory and reviewed with a physician or certified genetic counselor.",
        "",
        "## 1. Run Summary",
        f"- **Personal normalized variant records scanned (GRCh38 primary contigs):** {stats.get('total_variant_records', 0):,} ({stats.get('non_ref_calls', 0):,} non-reference genotypes)",
        f"- **Incidental / Autosomal Dominant / ACMG SF {ACMG_SF_VERSION} P/LP hits (after artifact filter):** {len(clean_incidental)} *(raw: {counts['incidental_ad_acmg']})*",
        f"- **Unreported rare `HIGH`-impact LoF variants in ACMG SF {ACMG_SF_VERSION} genes:** {len(novel_acmg_lofs)}",
        f"- **Homozygous / Hemizygous / Potential Compound-Het P/LP hits:** {len(clean_hom_hemi)}",
        f"- **Autosomal Recessive / X-linked Carrier (Het `0/1`) P/LP hits (after artifact filter):** {len(clean_carrier)} *(raw: {counts['carrier_ar_xl']})*",
        f"- **Other P/LP hits (Uncertain/Unmapped MOI, after artifact filter):** {len(clean_unknown)} *(raw: {counts['plp_unknown_moi']})*",
        f"- **Filtered likely artifacts / common polymorphisms (BA1/BS1 or 0★ legacy):** {len(filtered_artifacts)} (listed in Section 6)",
        f"- **Conflicting-interpretation hits (with P/LP assertions, gnomAD-annotated):** {counts['conflicting_plp']} (see `snp_indel_conflicting_plp.tsv`)",
        f"- **Risk factor / Drug response hits:** {counts['risk_or_drug']} (see `snp_indel_risk_or_drug.tsv`)",
        f"- **CNV / SV gene-body overlaps with ClinVar P/LP genes:** {len(cnv_sv_ok)} non-artifact ({len(cnv_sv_flagged)} flagged LOH/mega-SV artifacts; see `deep_dive_report.md` for coding-exon filtered view)",
        "",
        f"## 2. Incidental / Autosomal Dominant / ACMG SF {ACMG_SF_VERSION} Findings",
        f"_ClinVar P/LP hits in ACMG SF {ACMG_SF_VERSION} ({len(ACMG_SF_GENES)} genes) or curated AD / AD-AR / MT genes, excluding automated BA1/BS1 artifacts._",
        "",
        render_variant_table(clean_incidental),
        f"### 2b. Unreported Rare `HIGH`-Impact LoF Candidates in ACMG SF {ACMG_SF_VERSION} Genes",
        "_Frameshift, stop-gained, or canonical splice variants (`PASS`, gnomAD AF < 0.1% or absent) not yet cataloged as ClinVar P/LP._",
        "",
    ]

    if novel_acmg_lofs:
        md.append("| Gene | Variant (GRCh38) | GT | DP/GQ/AB | Consequence | HGVS | gnomAD v4 AF |")
        md.append("|---|---|---|---|---|---|---|")
        for r in novel_acmg_lofs:
            md.append(
                f"| **{r['gene']}** | `{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}` | `{r['gt']}` | "
                f"{r['dp']}/{r['gq']}/{r['ab']} | {r.get('consequence', '')} | "
                f"{r.get('hgvsp') or r.get('hgvsc') or '-'} | {format_af(r.get('gnomad_af'), r.get('gnomad_status'))} |"
            )
        md.append("")
    else:
        md.append("_None identified._\n")

    md.extend([
        "## 3. Homozygous / Hemizygous / Potential Compound-Het Candidates",
        "_Includes homozygous alternate (`1/1`), hemizygous (`1` on non-PAR X/Y), and genes carrying >=2 heterozygous P/LP variants (unphased: could be in **cis** rather than **trans**)._",
        "",
        render_variant_table(clean_hom_hemi),
        "## 4. Autosomal Recessive / X-Linked Carrier Candidates (Heterozygous `0/1`)",
        "_Single heterozygous P/LP variants in curated Autosomal Recessive (AR) or X-linked (XL) genes._",
        "",
        render_variant_table(clean_carrier),
        "## 5. Other ClinVar P/LP Hits (Uncertain / Uncurated MOI in GenCC)",
        render_variant_table(clean_unknown),
        "## 6. Filtered Likely Artifacts & Common Polymorphisms (Audit Trail)",
        "_Automatically demoted from Sections 2–5 based on ACMG BA1/BS1 allele frequency rules (AF ≥ 5% overall or ≥ 1% for AD), benign blood-group traits, or 0★ legacy somatic entries._",
        "",
        render_variant_table(filtered_artifacts, include_filter_reason=True),
        "## 7. CNV & SV Clinical Gene Overlaps (Gene-Body Level)",
        "_See `deep_dive_report.md` for GENCODE v46 coding-exon (CDS) filtered CNV/SV results._",
        "",
        render_cnv_sv_table(cnv_sv_ok, max_rows=30),
        "## 8. Known 30x Short-Read WGS Blind Spots",
        "A negative result above does **not** rule out carrier status in regions poorly resolved by 150 bp short-read WGS: "
        "**SMN1/SMN2**, **FMR1**, **HBA1/HBA2**, **CYP21A2**, **GBA1**, **PMS2**, **PKD1**, **CYP2D6**, and **Repeat Expansions** "
        "(*HTT*, *DMPK*, *FXN*, *C9orf72*, *ATXN**). See `docs/wgs_blind_spots.md` for details.",
    ])

    report_path = reports_dir / "carrier_screening_report.md"
    report_path.write_text("\n".join(md), encoding="utf-8")
    with open(reports_dir / "carrier_screen_enriched.json", "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    print(f"[Annotate: {args.sample}] Report written to {report_path}")


if __name__ == "__main__":
    main()
