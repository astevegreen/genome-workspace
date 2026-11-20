#!/usr/bin/env python3
"""
Step 7: Deep-Dive Clinical & QC Analysis
  1. Exonic CNV/SV Filtering: intersect Canvas CNVs & Manta SVs against GENCODE v46
     coding exons (CDS) of ClinVar P/LP & ACMG SF v3.3 genes, requiring >=50%
     reciprocal Canvas overlap for Manta DEL/DUP > 1 Mb.
  2. gVCF Coding Coverage Audit: compute per-gene CDS coverage (>=10x PASS using
     conservative MinDP, >=1x, and 0x/no-call) from snp_indel.genome.vcf.gz across
     ACMG SF v3.3 (84 genes) + carrier core + blind spots + sample-specific carrier hits.
  3. Categorized Canonical Marker Audit: clinical risk/carrier markers, CPIC Tier-1
     pharmacogenomics (PGx), and common/non-actionable traits, plus ClinVar drug/risk hits.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.constants import (
    ACMG_SF_GENES,
    ACMG_SF_VERSION,
    CARRIER_AND_BLINDSPOT_GENES,
    KNOWN_POLYMORPHIC_SV_GENES,
    PRIMARY_CHROMS_SET,
    TARGETED_MARKERS,
)
from lib.genomics import (
    compute_gene_gvcf_coverage,
    derive_apoe_diplotype,
    get_reference_metadata,
    load_gencode_cds,
    reciprocal_overlap,
    record_effective_min_dp,
    record_is_pass,
)
from lib.pgx import interpret_pgx_profile, render_pgx_medication_card

REPO_ROOT = SCRIPT_DIR.parent
REFERENCE_DIR = REPO_ROOT / "data" / "reference"


def load_clinvar_plp_genes(clinvar_vcf: Path) -> set[str]:
    genes: set[str] = set()
    with pysam.VariantFile(str(clinvar_vcf)) as cvf:
        for rec in cvf:
            clnsig = str(rec.info.get("CLNSIG") or "").lower()
            if "conflicting" in clnsig or "pathogenic" not in clnsig:
                continue
            rev = str(rec.info.get("CLNREVSTAT") or "")
            if "no_assertion" in rev or "no_classification" in rev:
                continue
            gi = str(rec.info.get("GENEINFO") or "")
            for item in gi.split("|"):
                sym = item.split(":")[0].strip()
                if sym:
                    genes.add(sym)
    return genes


def _sv_annotation_note(genes: list[str], source: str) -> str:
    poly_notes = [KNOWN_POLYMORPHIC_SV_GENES[g] for g in genes if g in KNOWN_POLYMORPHIC_SV_GENES]
    if poly_notes:
        return poly_notes[0]
    if source == "Canvas_CNV":
        return "Exonic CNV overlap"
    return "Exonic SV overlap (verify in IGV/BAM if actionable)"


def run_exonic_cnv_sv_filter(
    cnv_vcf: Path,
    sv_vcf: Path,
    cds_by_chrom: dict[str, list[tuple[int, int, str]]],
    clinvar_genes: set[str],
) -> list[dict[str, Any]]:
    """Intersect Canvas CNVs & Manta SVs against GENCODE v46 CDS of clinical genes."""
    clinical_genes = clinvar_genes | ACMG_SF_GENES
    clin_cds: dict[str, list[tuple[int, int, str]]] = {
        chrom: [x for x in ivs if x[2] in clinical_genes]
        for chrom, ivs in cds_by_chrom.items()
    }

    def overlapping_cds_genes(chrom: str, start: int, end: int) -> list[str]:
        found: set[str] = set()
        for c_start, c_end, gene in clin_cds.get(chrom, []):
            if c_end <= start:
                continue
            if c_start >= end:
                break
            found.add(gene)
        return sorted(found)

    results: list[dict[str, Any]] = []
    canvas_pass_cnvs: dict[str, list[tuple[int, int, str]]] = defaultdict(list)

    # 1. Canvas PASS CNVs (LOSS / GAIN)
    with pysam.VariantFile(str(cnv_vcf)) as vf:
        for rec in vf:
            if not rec.alts or rec.alts == (".",) or "Canvas:REF" in (rec.id or ""):
                continue
            flt = ";".join(rec.filter.keys()) if len(rec.filter) > 0 else "PASS"
            svtype = str(rec.info.get("SVTYPE") or "")
            if flt != "PASS" or svtype == "LOH":
                continue
            start0 = rec.pos - 1
            end = rec.stop or rec.info.get("END") or rec.pos
            canvas_pass_cnvs[rec.chrom].append((start0, end, rec.id or ""))
            genes = overlapping_cds_genes(rec.chrom, start0, end)
            if not genes:
                continue
            sample = next(iter(rec.samples.values()))
            acmg = [g for g in genes if g in ACMG_SF_GENES]
            results.append({
                "source": "Canvas_CNV",
                "chrom": rec.chrom,
                "start": rec.pos,
                "end": end,
                "length_bp": end - rec.pos,
                "svtype": "LOSS" if "LOSS" in (rec.id or "") else "GAIN",
                "call_id": rec.id,
                "gt_or_cn": f"CN={sample.get('CN')}",
                "qual_or_qs": round(float(sample.get("QS") or rec.qual or 0), 1),
                "support": f"BC={sample.get('BC')};RC={sample.get('RC')}",
                "acmg_cds_genes": ",".join(acmg),
                "n_cds_genes": len(genes),
                "cds_genes": ",".join(genes),
                "note": _sv_annotation_note(genes, "Canvas_CNV"),
            })

    # 2. Manta PASS SVs overlapping CDS (require >=50% reciprocal Canvas support for DEL/DUP > 1 Mb)
    dropped_mega = 0
    with pysam.VariantFile(str(sv_vcf)) as vf:
        for rec in vf:
            if rec.chrom not in PRIMARY_CHROMS_SET:
                continue
            flt = ";".join(rec.filter.keys()) if len(rec.filter) > 0 else "PASS"
            if flt != "PASS":
                continue
            svtype = str(rec.info.get("SVTYPE") or "")
            if svtype == "BND":
                alt_str = ",".join(rec.alts or ())
                if any(k in alt_str for k in ("decoy", "random", "Un_")):
                    continue
                start0, end, length = rec.pos - 1, rec.pos, 1
            else:
                start0 = rec.pos - 1
                end = rec.stop or rec.info.get("END") or rec.pos
                svlen_raw = rec.info.get("SVLEN")
                if isinstance(svlen_raw, (tuple, list)) and svlen_raw:
                    length = abs(int(svlen_raw[0]))
                elif svlen_raw is not None:
                    length = abs(int(svlen_raw))
                else:
                    length = abs(end - rec.pos)

            if svtype in ("DEL", "DUP") and length > 1_000_000:
                max_rec_ov = max(
                    (
                        reciprocal_overlap(start0, end, c_s, c_e)
                        for c_s, c_e, _ in canvas_pass_cnvs.get(rec.chrom, [])
                    ),
                    default=0.0,
                )
                if max_rec_ov < 0.50:
                    dropped_mega += 1
                    continue

            genes = overlapping_cds_genes(rec.chrom, start0, end)
            if not genes:
                continue
            sample = next(iter(rec.samples.values()))
            gt = sample.get("GT")
            gt_str = "/".join("." if a is None else str(a) for a in gt) if gt else "./."
            acmg = [g for g in genes if g in ACMG_SF_GENES]
            results.append({
                "source": "Manta_SV",
                "chrom": rec.chrom,
                "start": rec.pos,
                "end": end,
                "length_bp": length,
                "svtype": svtype,
                "call_id": rec.id,
                "gt_or_cn": gt_str,
                "qual_or_qs": round(float(rec.qual), 1) if rec.qual is not None else None,
                "support": f"PR={sample.get('PR')};SR={sample.get('SR')}",
                "acmg_cds_genes": ",".join(acmg),
                "n_cds_genes": len(genes),
                "cds_genes": ",".join(genes),
                "note": _sv_annotation_note(genes, "Manta_SV"),
            })

    print(
        f"[Exonic CNV/SV] {len(results)} coding-exon clinical overlaps retained "
        f"({dropped_mega} unconfirmed Manta mega-SVs >1 Mb filtered out)."
    )
    return results


def run_gvcf_coverage_audit(
    gvcf_path: Path,
    cds_by_gene: dict[str, tuple[str, list[tuple[int, int]]]],
    sample_carrier_genes: set[str],
) -> list[dict[str, Any]]:
    """Compute per-gene CDS coverage >=10x (PASS, conservative MinDP), >=1x, and 0x."""
    audit_genes: dict[str, str] = {}
    for g in sorted(ACMG_SF_GENES):
        audit_genes[g] = f"ACMG_SF_{ACMG_SF_VERSION}"
    for g, cat in CARRIER_AND_BLINDSPOT_GENES.items():
        audit_genes[g] = f"ACMG+{cat}" if g in audit_genes else cat
    for g in sorted(sample_carrier_genes):
        if g not in audit_genes:
            audit_genes[g] = "Sample_Carrier_Hit"
        elif "Carrier_Hit" not in audit_genes[g]:
            audit_genes[g] = f"{audit_genes[g]}+Sample_Carrier_Hit"

    rows: list[dict[str, Any]] = []
    with pysam.VariantFile(str(gvcf_path)) as gvf:
        for gene, panel in sorted(audit_genes.items()):
            entry = cds_by_gene.get(gene)
            if not entry:
                rows.append({
                    "gene": gene,
                    "panel": panel,
                    "chrom": "NOT_IN_GTF",
                    "n_cds_exons": 0,
                    "cds_bp": 0,
                    "pct_ge_10x": 0.0,
                    "pct_ge_1x": 0.0,
                    "pct_0x_nocall": 100.0,
                })
                continue
            chrom, exons = entry
            cov = compute_gene_gvcf_coverage(gvf, chrom, exons, min_pass_dp=10)
            rows.append({
                "gene": gene,
                "panel": panel,
                "chrom": chrom,
                "n_cds_exons": len(exons),
                "cds_bp": cov["cds_bp"],
                "pct_ge_10x": cov["pct_ge_10x"],
                "pct_ge_1x": cov["pct_ge_1x"],
                "pct_0x_nocall": cov["pct_0x_nocall"],
            })

    rows.sort(key=lambda x: (x["pct_ge_10x"], x["gene"]))
    print(f"[Coverage Audit] Audited {len(rows)} ACMG/carrier/blind-spot genes (strict PASS + MinDP >= 10x).")
    return rows


def run_targeted_marker_genotyping(
    gvcf_path: Path,
    norm_vcf_path: Path,
) -> tuple[list[dict[str, Any]], str]:
    """Directly genotype canonical clinical, PGx, and common-trait markers in VCF + gVCF."""
    rows: list[dict[str, Any]] = []
    apoe_gts: dict[str, str] = {}

    with pysam.VariantFile(str(norm_vcf_path)) as nvf, pysam.VariantFile(str(gvcf_path)) as gvf:
        for category, gene, rsid, chrom, pos, ref, alt, desc in TARGETED_MARKERS:
            gt_str = "0/0 (REF)"
            dp = None
            gq = None
            ab = None
            flt = "PASS"
            status = "REF_CONFIRMED"

            found_var = False
            for rec in nvf.fetch(chrom, pos - 1, pos + len(ref)):
                if rec.pos == pos and rec.ref == ref and rec.alts and alt in rec.alts:
                    sample = next(iter(rec.samples.values()))
                    gt = sample.get("GT")
                    alt_idx = rec.alts.index(alt) + 1
                    if not gt or alt_idx not in gt:
                        continue
                    gt_str = "/".join("." if a is None else str(a) for a in gt)
                    dp = sample.get("DP")
                    if dp is None:
                        dp = rec.info.get("DP")
                    gq = sample.get("GQ")
                    ad = sample.get("AD")
                    if ad and len(ad) > alt_idx and (ad[0] + ad[alt_idx]) > 0:
                        ab = round(ad[alt_idx] / (ad[0] + ad[alt_idx]), 3)
                    flt = ";".join(rec.filter.keys()) if len(rec.filter) > 0 else "PASS"
                    status = "VARIANT_PRESENT" if record_is_pass(rec) else "VARIANT_LOW_CONFIDENCE"
                    found_var = True
                    break

            if not found_var:
                best_pass_dp = None
                best_any_dp = None
                best_any_flt = "NO_CALL"
                for rec in gvf.fetch(chrom, pos - 1, pos):
                    r_start = rec.pos
                    r_end = rec.stop or rec.info.get("END") or rec.pos
                    if r_start <= pos <= r_end:
                        sample = next(iter(rec.samples.values()))
                        d = record_effective_min_dp(rec, sample)
                        is_pass = record_is_pass(rec)
                        if is_pass and (best_pass_dp is None or d > best_pass_dp):
                            best_pass_dp = d
                        if best_any_dp is None or d > best_any_dp:
                            best_any_dp = d
                            best_any_flt = ";".join(rec.filter.keys()) if len(rec.filter) > 0 else "PASS"
                if best_pass_dp is not None:
                    dp = best_pass_dp
                    flt = "PASS"
                else:
                    dp = best_any_dp
                    flt = best_any_flt
                if dp is None:
                    gt_str = "./."
                    status = "NOT_COVERED"
                elif dp < 10 or flt != "PASS":
                    status = "REF_LOW_DP"

            if rsid in ("rs429358", "rs7412"):
                apoe_gts[rsid] = gt_str

            rows.append({
                "category": category,
                "gene": gene,
                "rsid": rsid,
                "locus": f"{chrom}:{pos} {ref}>{alt}",
                "status": status,
                "gt": gt_str,
                "dp": dp,
                "gq": gq,
                "ab": ab,
                "filter": flt,
                "description": desc,
            })

    apoe_diplotype = derive_apoe_diplotype(
        apoe_gts.get("rs429358", ""), apoe_gts.get("rs7412", "")
    )
    return rows, apoe_diplotype


def summarize_risk_and_drug(risk_tsv: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not risk_tsv.exists() or risk_tsv.stat().st_size <= 2:
        return pd.DataFrame(), pd.DataFrame()
    df = pd.read_csv(risk_tsv, sep="\t", dtype=str)
    if df.empty:
        return pd.DataFrame(), pd.DataFrame()
    df["clinvar_stars_num"] = pd.to_numeric(df["clinvar_stars"], errors="coerce").fillna(0).astype(int)
    drug_df = df[df["clnsig"].str.contains("drug_response", case=False, na=False)].copy()
    drug_df = drug_df.sort_values(["clinvar_stars_num", "gene"], ascending=[False, True])
    risk_df = df[
        ~df["clnsig"].str.contains("drug_response", case=False, na=False)
        & df["clnsig"].str.contains("risk_factor|protective", case=False, na=False)
    ].copy()
    risk_df = risk_df.sort_values(["clinvar_stars_num", "gene"], ascending=[False, True])
    return drug_df, risk_df


def _render_marker_table(rows: list[dict[str, Any]]) -> list[str]:
    lines = [
        "| Gene | rsID | Locus (GRCh38) | Status | Genotype | DP | GQ | AB | Description |",
        "|---|---|---|---|---|---:|---:|---:|---|",
    ]
    for r in rows:
        bold = "**" if r["status"] == "VARIANT_PRESENT" else ""
        lines.append(
            f"| {bold}{r['gene']}{bold} | {r['rsid']} | `{r['locus']}` | {bold}{r['status']}{bold} | "
            f"`{r['gt']}` | {r['dp'] if r['dp'] is not None else '-'} | "
            f"{r['gq'] if r['gq'] is not None else '-'} | {r['ab'] if r['ab'] is not None else '-'} | {r['description']} |"
        )
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description="Run deep-dive exonic CNV/SV, gVCF coverage & PGx audit.")
    parser.add_argument("--sample", required=True, help="Sample label from config/samples.tsv (e.g. partner_a, partner_b)")
    args = parser.parse_args()

    processed_dir = REPO_ROOT / "data" / "processed" / args.sample
    reports_dir = REPO_ROOT / "results" / "reports" / args.sample
    reports_dir.mkdir(parents=True, exist_ok=True)

    cds_bed = REFERENCE_DIR / "gencode.v46.cds.nochr.bed.gz"
    clinvar_vcf = REFERENCE_DIR / "clinvar_GRCh38.norm.vcf.gz"
    cnv_vcf = processed_dir / "cnv.vcf.gz"
    sv_vcf = processed_dir / "sv.vcf.gz"
    gvcf = processed_dir / "snp_indel.genome.vcf.gz"
    norm_vcf = processed_dir / "snp_indel.variants.norm.vcf.gz"
    risk_tsv = reports_dir / "snp_indel_risk_or_drug.tsv"
    enriched_json = reports_dir / "carrier_screen_enriched.json"

    ref_meta = get_reference_metadata(REFERENCE_DIR)
    raw_cnv_sv_non_artifact: int | None = None
    sample_carrier_genes: set[str] = set()
    clinvar_genes: set[str] = set()
    if enriched_json.exists():
        with open(enriched_json, encoding="utf-8") as f:
            edata = json.load(f)
        raw_cnv_sv_non_artifact = edata.get("cnv_sv_non_artifact")
        if isinstance(edata.get("clinvar_plp_genes"), list):
            clinvar_genes = set(edata["clinvar_plp_genes"])
        for bname in ("carrier_ar_xl", "incidental_ad_acmg", "hom_hemi_comphet", "plp_unknown_moi"):
            for r in edata.get("buckets", {}).get(bname, []):
                if r.get("gene") and r.get("artifact_flag", "OK") == "OK":
                    sample_carrier_genes.add(r["gene"])

    print(f"[Step 7: {args.sample}] Loading GENCODE v46 CDS intervals...")
    cds_by_chrom, cds_by_gene = load_gencode_cds(cds_bed)
    print(f"[Step 7] Loaded CDS for {len(cds_by_gene):,} genes.")

    if clinvar_genes:
        print(f"[Step 7] Reusing {len(clinvar_genes):,} ClinVar P/LP (>=1*) genes from Step 4.")
    else:
        print("[Step 7] Indexing ClinVar P/LP (>=1*) gene symbols...")
        clinvar_genes = load_clinvar_plp_genes(clinvar_vcf)

    exonic_svs = run_exonic_cnv_sv_filter(cnv_vcf, sv_vcf, cds_by_chrom, clinvar_genes)
    pd.DataFrame(exonic_svs).to_csv(
        reports_dir / "cnv_sv_exonic_clinical.tsv", sep="\t", index=False
    )

    cov_rows = run_gvcf_coverage_audit(gvcf, cds_by_gene, sample_carrier_genes)
    cov_df = pd.DataFrame(cov_rows)
    cov_df.to_csv(reports_dir / "gvcf_coding_coverage_audit.tsv", sep="\t", index=False)

    marker_rows, apoe_diplotype = run_targeted_marker_genotyping(gvcf, norm_vcf)
    pd.DataFrame(marker_rows).to_csv(
        reports_dir / "targeted_markers.tsv", sep="\t", index=False
    )
    pgx_entries = interpret_pgx_profile(marker_rows)
    pgx_card_md = render_pgx_medication_card(args.sample, ref_meta, pgx_entries)
    pgx_card_path = reports_dir / "pgx_medication_card.md"
    pgx_card_path.write_text(pgx_card_md, encoding="utf-8")
    pd.DataFrame(pgx_entries).to_csv(
        reports_dir / "pgx_medication_card.tsv", sep="\t", index=False
    )
    print(f"[Step 7: {args.sample}] CPIC Medication Card written to {pgx_card_path}")
    drug_df, risk_df = summarize_risk_and_drug(risk_tsv)

    acmg_cov = cov_df[cov_df["panel"].str.contains("ACMG")]
    low_cov = cov_df[cov_df["pct_ge_10x"] < 95.0]
    gene_body_note = (
        f" (down from {raw_cnv_sv_non_artifact} gene-body overlaps after restricting to CDS exons and filtering unconfirmed Manta mega-SVs >1 Mb lacking ≥50% reciprocal Canvas CNV support)"
        if raw_cnv_sv_non_artifact is not None
        else ""
    )

    md = [
        f"# Personal 30x WGS Deep-Dive Audit ({args.sample.title()}): Exonic CNV/SV, gVCF Coverage & PGx/Risk Markers",
        "",
        f"**Run Date:** {ref_meta['run_date']} | **ClinVar:** {ref_meta['clinvar_date']} | "
        f"**GenCC:** {ref_meta['gencc_date']} | **ACMG SF:** {ref_meta['acmg_sf_version']} | **GENCODE:** {ref_meta['gencode_release']}",
        "",
        f"## 1. Coding-Exon (GENCODE v46 CDS) CNV & SV Clinical Overlaps",
        f"- **Retained coding-exon overlaps in ClinVar P/LP / ACMG SF {ACMG_SF_VERSION} genes:** **{len(exonic_svs)}**{gene_body_note}.",
        "",
    ]

    if exonic_svs:
        md.append("| Source | Type | Region (GRCh38) | Size | GT/CN | Support | ACMG CDS Genes | CDS Genes Overlapped | Note |")
        md.append("|---|---|---|---:|---|---|---|---|---|")
        for r in exonic_svs:
            sz = f"{r['length_bp']/1000:,.1f} kb" if r["length_bp"] >= 1000 else f"{r['length_bp']} bp"
            gene_list = [g for g in str(r["cds_genes"]).split(",") if g]
            genes_preview = ",".join(gene_list[:15]) + ("..." if len(gene_list) > 15 else "")
            md.append(
                f"| {r['source']} | {r['svtype']} | `{r['chrom']}:{r['start']:,}-{r['end']:,}` | "
                f"{sz} | `{r['gt_or_cn']}` | {r['support']} | **{r['acmg_cds_genes'] or '-'}** | "
                f"{genes_preview} ({r['n_cds_genes']}) | {r['note']} |"
            )
    else:
        md.append(f"_No coding-exon CNV/SV overlaps detected in ClinVar P/LP or ACMG SF {ACMG_SF_VERSION} genes._")

    md.extend([
        "",
        f"## 2. gVCF Coding Coverage Audit (ACMG SF {ACMG_SF_VERSION} + Carrier & Blind-Spot Panel)",
        f"- **ACMG SF {ACMG_SF_VERSION} ({len(ACMG_SF_GENES)} genes):** median **{acmg_cov['pct_ge_10x'].median():.2f}%** of CDS bases at ≥10x (`FILTER == PASS`, conservative `MinDP`); **{int((acmg_cov['pct_ge_10x'] >= 98.0).sum())}/{len(ACMG_SF_GENES)}** genes have ≥98% CDS at ≥10x.",
        f"- **Genes with <95% CDS at ≥10x (`PASS`)** (confirms where short-read WGS has reduced sensitivity / blind spots):",
        "",
        "| Gene | Panel | Chrom | CDS Exons | CDS bp | ≥10x (`PASS`) % | ≥1x % | 0x / No-Call % |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ])
    for _, r in low_cov.iterrows():
        md.append(
            f"| **{r['gene']}** | {r['panel']} | {r['chrom']} | {r['n_cds_exons']} | "
            f"{r['cds_bp']:,} | **{r['pct_ge_10x']:.2f}%** | {r['pct_ge_1x']:.2f}% | {r['pct_0x_nocall']:.2f}% |"
        )

    md.extend([
        "",
        f"_Full {len(cov_df)}-gene table saved in `results/reports/{args.sample}/gvcf_coding_coverage_audit.tsv`._",
        "",
        "## 3. Targeted Canonical Clinical, Pharmacogenomic & Trait Markers",
        f"- **Derived *APOE* Diplotype:** **`{apoe_diplotype}`**",
        "",
        "### 3a. Clinical Risk & Core Carrier Markers",
        *_render_marker_table([r for r in marker_rows if r["category"] == "clinical"]),
        "",
        "### 3b. Pharmacogenomics (CPIC Tier-1 & Drug Safety Markers)",
        *_render_marker_table([r for r in marker_rows if r["category"] == "pgx"]),
        "",
        "### 3c. Common / Non-Actionable Traits (Informational Only)",
        *_render_marker_table([r for r in marker_rows if r["category"] == "common_trait"]),
        "",
        f"## 4. ClinVar Drug-Response Hits (`{len(drug_df)}` total in `snp_indel_risk_or_drug.tsv`)",
        "| Gene | rsID | GT | QC (DP/GQ/AB) | ClinVar Stars | Drug / Trait |",
        "|---|---|---|---|---:|---|",
    ])
    for _, r in drug_df.head(30).iterrows():
        cond = str(r.get("clndn", "")).replace("|", "; ").replace("_", " ")
        if len(cond) > 90:
            cond = cond[:87] + "..."
        md.append(
            f"| **{r.get('gene', '.')}** | {r.get('rsid', '.')} | `{r.get('gt', '.')}` | "
            f"{r.get('qc_tier')} ({r.get('dp')}/{r.get('gq')}/{r.get('ab')}) | {r['clinvar_stars_num']}★ | {cond} |"
        )

    md.extend([
        "",
        f"## 5. ClinVar Risk-Factor / Protective Hits (`{len(risk_df)}` total; showing ≥1★)",
        "| Gene | rsID | GT | ClinVar Sig | Stars | Condition / Trait |",
        "|---|---|---|---|---:|---|",
    ])
    if not risk_df.empty:
        for _, r in risk_df[risk_df["clinvar_stars_num"] >= 1].head(30).iterrows():
            cond = str(r.get("clndn", "")).replace("|", "; ").replace("_", " ")
            if len(cond) > 90:
                cond = cond[:87] + "..."
            md.append(
                f"| **{r.get('gene', '.')}** | {r.get('rsid', '.')} | `{r.get('gt', '.')}` | "
                f"{r.get('clnsig', '')} | {r['clinvar_stars_num']}★ | {cond} |"
            )

    out_md = reports_dir / "deep_dive_report.md"
    out_md.write_text("\n".join(md), encoding="utf-8")
    print(f"[Step 7: {args.sample}] Deep-dive report written to {out_md}")


if __name__ == "__main__":
    main()
