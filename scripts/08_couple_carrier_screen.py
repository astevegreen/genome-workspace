#!/usr/bin/env python3
"""
Step 8: Joint Couple / Reproductive Carrier Screening & Auto-Generated Plain-Language Summary
Compares two normalized 30x WGS datasets for potential offspring & personal health risk:
  1. Chromosomal sex concordance QC (non-PAR X het rate & Y variant count).
  2. Gene-level Autosomal Recessive (AR) matching across SNV/indels + exonic CNVs/SVs,
     explicitly reporting BOTH SNV/indel matches AND shared exonic SV+SV overlaps.
  3. Partner-gene rare/novel coding & splice variant rescue via cached Ensembl VEP + gnomAD v4.
  4. Reciprocal gVCF coding coverage verification (>=10x PASS using conservative MinDP).
  5. Maternal X-linked & Mitochondrial (MT) screening.
  6. Autosomal Dominant (AD) & ACMG SF v3.3 screening (including unreported rare LoF rescue).
  7. Automated separation of BA1/BS1 false alarms & common polymorphisms into an audit section.
  8. Auto-generation of results/reports/couple/plain_language_summary.md so annual re-runs
     never leave the human-friendly summary out of sync.
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

from lib.api_client import CachedAnnotationClient, format_af
from lib.constants import ACMG_SF_GENES, ACMG_SF_VERSION, KNOWN_POLYMORPHIC_SV_GENES
from lib.genomics import (
    compute_gene_gvcf_coverage,
    derive_apoe_diplotype,
    get_reference_metadata,
    is_strong_conflicting,
    load_gencode_cds,
)

REPO_ROOT = SCRIPT_DIR.parent
REFERENCE_DIR = REPO_ROOT / "data" / "reference"


def check_sex_qc(norm_vcf: Path, expected_sex: str) -> dict[str, Any]:
    x_het = 0
    x_hom_or_hemi = 0
    y_vars = 0
    with pysam.VariantFile(str(norm_vcf)) as vf:
        for rec in vf.fetch("X"):
            if rec.pos < 2_800_000 or rec.pos > 155_700_000:
                continue
            gt = next(iter(rec.samples.values())).get("GT")
            if not gt:
                continue
            if len(gt) == 2 and gt[0] != gt[1]:
                x_het += 1
            elif any(a is not None and a > 0 for a in gt):
                x_hom_or_hemi += 1
        try:
            for rec in vf.fetch("Y"):
                gt = next(iter(rec.samples.values())).get("GT")
                if gt and any(a is not None and a > 0 for a in gt):
                    y_vars += 1
        except ValueError:
            pass
    total_x = x_het + x_hom_or_hemi
    het_ratio = round(x_het / total_x, 4) if total_x else 0.0
    inferred = "F" if het_ratio > 0.25 else "M"
    return {
        "expected_sex": expected_sex,
        "inferred_sex": inferred,
        "concordant": inferred == expected_sex,
        "non_par_x_het": x_het,
        "non_par_x_hom_hemi": x_hom_or_hemi,
        "non_par_x_het_ratio": het_ratio,
        "y_variant_count": y_vars,
    }


def collect_partner_clinical_hits(
    enriched_json: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Return (clean_hits, filtered_artifact_hits, novel_acmg_lof_hits)."""
    with open(enriched_json, encoding="utf-8") as f:
        data = json.load(f)
    buckets = data["buckets"]
    seen: set[str] = set()
    clean_hits: list[dict[str, Any]] = []
    artifact_hits: list[dict[str, Any]] = []

    def _add(row: dict[str, Any], tier: str) -> None:
        vid = row["variant_id"]
        if vid in seen:
            return
        seen.add(vid)
        entry = {**row, "evidence_tier": tier}
        if row.get("artifact_flag", "OK") == "OK":
            clean_hits.append(entry)
        else:
            artifact_hits.append(entry)

    for bname in ("incidental_ad_acmg", "hom_hemi_comphet", "carrier_ar_xl", "plp_unknown_moi"):
        for r in buckets.get(bname, []):
            _add(r, "CLINVAR_PLP")

    for r in buckets.get("conflicting_plp", []):
        if is_strong_conflicting(r):
            _add(r, "STRONG_CONFLICTING_PLP")

    novel_lofs = data.get("novel_acmg_lof") or []
    return clean_hits, artifact_hits, novel_lofs


def scan_partner_gene_for_rare_coding_variants(
    partner_norm_vcf: Path,
    gene: str,
    chrom: str,
    exons: list[tuple[int, int]],
    api: CachedAnnotationClient,
    splice_pad: int = 10,
) -> list[dict[str, Any]]:
    """Extract all coding/splice variants in partner for `gene` and annotate via cached VEP + gnomAD."""
    candidates: list[dict[str, Any]] = []
    vep_inputs: list[tuple[str, str, str, int, str, str]] = []
    seen_vid: set[str] = set()

    with pysam.VariantFile(str(partner_norm_vcf)) as vf:
        for ex_s, ex_e in exons:
            fetch_s = max(0, ex_s - splice_pad)
            fetch_e = ex_e + splice_pad
            for rec in vf.fetch(chrom, fetch_s, fetch_e):
                if not rec.alts:
                    continue
                sample = next(iter(rec.samples.values()))
                gt = sample.get("GT")
                if not gt or all(a is None or a == 0 for a in gt):
                    continue
                ad = sample.get("AD")
                present_alt_idxs = {a for a in gt if a is not None and a > 0}
                gt_str = "/".join("." if a is None else str(a) for a in gt)
                for idx in present_alt_idxs:
                    alt = rec.alts[idx - 1]
                    vid = f"{rec.chrom}-{rec.pos}-{rec.ref}-{alt}"
                    if vid in seen_vid:
                        continue
                    seen_vid.add(vid)
                    ab = (
                        round(ad[idx] / (ad[0] + ad[idx]), 3)
                        if ad and len(ad) > idx and (ad[0] + ad[idx]) > 0
                        else None
                    )
                    candidates.append({
                        "gene": gene,
                        "variant_id": vid,
                        "chrom": rec.chrom,
                        "pos": rec.pos,
                        "ref": rec.ref,
                        "alt": alt,
                        "rsid": rec.id or ".",
                        "gt": gt_str,
                        "dp": sample.get("DP"),
                        "gq": sample.get("GQ"),
                        "ab": ab,
                        "filter": ";".join(rec.filter.keys()) if len(rec.filter) > 0 else "PASS",
                    })
                    vep_inputs.append((vid, gene, rec.chrom, rec.pos, rec.ref, alt))

    if not candidates:
        return []

    vep_map = api.query_ensembl_vep_batch(vep_inputs)
    rare_protein_altering: list[dict[str, Any]] = []
    for c in candidates:
        vinfo = vep_map.get(c["variant_id"], {})
        if vinfo.get("impact") not in ("HIGH", "MODERATE"):
            continue
        g_info = api.query_gnomad_variant(c["variant_id"])
        af = g_info.get("gnomad_af")
        if af is None or af < 0.01:
            rare_protein_altering.append({**c, **vinfo, **g_info})
    return rare_protein_altering


def _load_bam_summary(rep_dir: Path) -> dict[str, Any] | None:
    p = rep_dir / "bam_blind_spots.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _sv_bam_verdict_by_gene(bam_data: dict[str, Any] | None) -> dict[str, list[str]]:
    """Map gene -> list of concise BAM SV verification strings (e.g. 'DEL ratio=0.01 CONFIRMED_HOM_DEL')."""
    out: dict[str, list[str]] = defaultdict(list)
    if not bam_data:
        return out
    for s in bam_data.get("sv_verification") or []:
        desc = f"{s['svtype']} ratio={s['depth_ratio']:.2f} ({s['bam_verdict']})"
        for g in str(s.get("cds_genes") or "").split(","):
            g = g.strip()
            if g:
                out[g].append(desc)
    return out


def render_plain_language_summary(
    male_lbl: str,
    fem_lbl: str,
    ref_meta: dict[str, str],
    shared_snv_matches: list[dict[str, Any]],
    shared_sv_only: list[dict[str, Any]],
    partner_rescue_rows: list[dict[str, Any]],
    maternal_x_snv: list[dict[str, Any]],
    maternal_x_sv: list[dict[str, Any]],
    maternal_mt_snv: list[dict[str, Any]],
    ad_acmg_hits: list[dict[str, Any]],
    novel_acmg_lofs: list[dict[str, Any]],
    filtered_artifacts: list[dict[str, Any]],
    tm_merged: pd.DataFrame,
    male_bam: dict[str, Any] | None = None,
    fem_bam: dict[str, Any] | None = None,
) -> str:
    """Auto-generate a clear, jargon-light summary synced to the live pipeline outputs."""
    m_name, f_name = male_lbl.title(), fem_lbl.title()
    tm_by_rs = {str(r["rsid"]): r for _, r in tm_merged.iterrows()}

    def _gt(rsid: str, who: str) -> str:
        r = tm_by_rs.get(rsid)
        return str(r[f"gt_{who}"]) if r is not None else "N/A"

    m_apoe = derive_apoe_diplotype(_gt("rs429358", male_lbl), _gt("rs7412", male_lbl))
    f_apoe = derive_apoe_diplotype(_gt("rs429358", fem_lbl), _gt("rs7412", fem_lbl))

    # Count non-polymorphic shared SV-only overlaps
    non_poly_shared_svs = [s for s in shared_sv_only if not s["is_known_polymorphic"]]
    n_actionable_ad = len(ad_acmg_hits) + len(novel_acmg_lofs)

    md = [
        f"# Genetic Carrier Screening & Genome Comparison Summary ({m_name} & {f_name})",
        "",
        f"**Analysis Date:** {ref_meta['run_date']}  ",
        f"**Reference Databases:** NCBI ClinVar GRCh38 (`{ref_meta['clinvar_date']}`), GenCC (`{ref_meta['gencc_date']}`), "
        f"ACMG SF `{ref_meta['acmg_sf_version']}`, GENCODE `{ref_meta['gencode_release']}`, gnomAD v4, Ensembl VEP  ",
        "*(Auto-generated by `scripts/08_couple_carrier_screen.py` — always synced with the latest pipeline run.)*",
        "",
        "> **Medical Disclaimer:** This is a research-grade bioinformatics summary prepared for personal education and discussion with a healthcare provider. "
        "It is **not** a clinical diagnostic test. Before making reproductive or medical decisions, any findings—and especially the known blind spots of standard "
        "short-read whole-genome sequencing listed in Section 6—should be confirmed through a CLIA-certified clinical laboratory and reviewed with a physician or board-certified genetic counselor.",
        "",
        "---",
        "",
        "## Executive Summary (The Bottom Line)",
        "",
    ]

    if not shared_snv_matches and not maternal_x_snv and not maternal_mt_snv and not novel_acmg_lofs:
        md.append(f"Overall, the comparison between {m_name}'s and {f_name}'s whole-genome sequences is **very reassuring for potential offspring**:")
    else:
        md.append(f"Summary of the joint comparison between {m_name}'s and {f_name}'s whole-genome sequences:")

    has_bam = bool(male_bam or fem_bam)
    bam_exec_note = (
        "5. **BAM-Resolved Blind Spots & Remaining Clinical Follow-Up:** Specialized BAM callers (**ExpansionHunter**, **SMNCopyNumberCaller**, **Cyrius**) have now been run on both 30x WGS alignments for **Fragile X (*FMR1*)**, **SMA (*SMN1*/*SMN2*)**, **repeat expansions**, and **full *CYP2D6* star alleles/CNVs/hybrids** (see Section 3 & Section 5). Three high-homology blind spots still require a clinical carrier panel (**Alpha-Thalassemia**, **CAH**, and **Gaucher**), and any BAM-based reproductive/PGx finding should still be confirmed clinically."
        if has_bam
        else "5. **Recommended Clinical Follow-Up:** Even 30x whole-genome sequencing has five well-known short-read blind spots (**Fragile X**, **SMA**, **Alpha-Thalassemia**, **CAH**, and **Gaucher**). Standard clinical carrier screening covers these gaps; uploading BAM/CRAM alignment files in the future will also unlock local callers for several of them."
    )
    sv_exec_note = (
        f"*(Additionally, {len(shared_sv_only)} gene(s) had shared short-read structural variant calls—{len(shared_sv_only) - len(non_poly_shared_svs)} in known common/HLA regions and {len(non_poly_shared_svs)} in other loci; all have now been audited for read-level depth/breakpoint support in the technical couple report.)*"
        if has_bam
        else f"*(Additionally, {len(shared_sv_only)} gene(s) had shared short-read structural variant calls—{len(shared_sv_only) - len(non_poly_shared_svs)} in known common/HLA regions and {len(non_poly_shared_svs)} in other loci, documented transparently in the technical couple report for future BAM/IGV verification.)*"
    )

    md.extend([
        "",
        f"1. **Shared Recessive Disease Genes (25% Offspring Risk):** **{len(shared_snv_matches)}** shared SNV/indel or SNV+SV carrier gene(s) found. {sv_exec_note}",
        f"2. **Maternal X-Linked & Mitochondrial Screening ({f_name}):** **{len(maternal_x_snv)}** X-linked SNV/indel(s), **{len(maternal_x_sv)}** X-linked exonic SV(s), and **{len(maternal_mt_snv)}** mitochondrial disease variant(s) found.",
        f"3. **Medically Actionable ACMG SF {ACMG_SF_VERSION} ({len(ACMG_SF_GENES)} Genes) & Dominant Conditions:** **{len(ad_acmg_hits)}** ClinVar/conflicting hit(s) and **{len(novel_acmg_lofs)}** unreported rare loss-of-function candidate(s) found across both partners.",
        f"4. **Automatic False-Alarm Filtering:** **{len(filtered_artifacts)}** common sequencing artifact(s) / benign blood-group / legacy 0-star entries were automatically flagged and moved to Section 4 so they don't cause false alarm.",
        bam_exec_note,
        "",
        "---",
        "",
        "## 1. Joint Offspring Risk Snapshot",
        "",
        "| Risk Category | Result | Plain-Language Meaning |",
        "|---|---|---|",
        f"| **Shared Autosomal Recessive Genes (SNV/Indel)** | **{len(shared_snv_matches)} shared gene(s)** | A child has a 25% risk only if *both* parents carry a disease-causing variant in the same recessive gene. |",
        f"| **Maternal X-Linked Conditions ({f_name})** | **{len(maternal_x_snv)} SNV / {len(maternal_x_sv)} exonic SV** | Conditions passed from mother to son (50% risk to boys). |",
        f"| **Maternal Mitochondrial (`MT`) Conditions ({f_name})** | **{len(maternal_mt_snv)} variant(s)** | Mitochondrial DNA is inherited exclusively from the mother. |",
        f"| **ACMG SF {ACMG_SF_VERSION} ({len(ACMG_SF_GENES)} Actionable Genes) & AD** | **{n_actionable_ad} flagged item(s)** | Includes both ClinVar P/LP and unreported rare gene-breaking (LoF) variants. |",
        "",
        "---",
        "",
        "## 2. Individual Carrier Variants & Partner Verification",
        "",
        "Almost everyone is a healthy silent carrier of 2–7 recessive conditions. For **every** non-X/MT gene where one partner carried a variant, we verified that the other partner's DNA was read cleanly (`>=10x PASS`) AND scanned every coding exon and splice junction for any rare unreported mutation.",
        "",
    ])

    for who_lbl, who_name, partner_name in ((male_lbl, m_name, f_name), (fem_lbl, f_name, m_name)):
        rows = [r for r in partner_rescue_rows if r["carrier_partner"] == who_lbl]
        md.append(f"### {who_name}'s Carrier / Strong Conflicting Variants (Checked in {partner_name})")
        if not rows:
            md.append("_None identified._\n")
            continue
        md.extend([
            f"| Gene | MOI | Variant (rsID) | Evidence | {partner_name}'s Coverage (`>=10x PASS`) | {partner_name}'s Rare Coding Variants / Exonic SVs |",
            "|---|---|---|---|---:|---|",
        ])
        for r in rows:
            n_rv = r["partner_rare_coding_count"]
            n_sv = r["partner_exonic_sv_count"]
            status_txt = (
                "**Cleared** (0 rare coding, 0 exonic SVs)"
                if n_rv == 0 and n_sv == 0
                else f"**Review:** {n_rv} rare coding, {n_sv} exonic SV(s)"
            )
            md.append(
                f"| **{r['gene']}** | {r['moi']} | `{r['carrier_variant']}` ({r['carrier_rsid']}) | "
                f"{r['carrier_evidence']} | **{r['partner_pct_ge_10x']:.1f}%** | {status_txt} |"
            )
        md.append("")

    md.extend([
        "---",
        "",
        "## 3. Personal Health, Medication (PGx) & Common Traits",
        "",
        "### A. *APOE* (Cholesterol & Brain Aging Modifier)",
        f"- **{m_name}:** `{m_apoe}`",
        f"- **{f_name}:** `{f_apoe}`",
        "",
        "### B. Pharmacogenomics (Medication Processing — CPIC Tier-1 Markers)",
        f"| Gene | Marker (rsID) | {m_name} | {f_name} | Clinical Note |",
        "|---|---|---|---|---|",
    ])
    for _, r in tm_merged[tm_merged["category"] == "pgx"].iterrows():
        m_gt = str(r[f"gt_{male_lbl}"])
        f_gt = str(r[f"gt_{fem_lbl}"])
        # Show all PGx rows where either partner carries a variant, plus summarize key toxicity genes
        if "REF" not in m_gt or "REF" not in f_gt or r["gene"] in ("DPYD", "TPMT", "NUDT15", "CYP2D6"):
            md.append(f"| **{r['gene']}** | {r['rsid']} | `{m_gt}` | `{f_gt}` | {r['description']} |")

    md.extend([
        "",
        "### C. Other Canonical Clinical & Common Traits",
        f"| Category | Gene | Marker (rsID) | {m_name} | {f_name} | Note |",
        "|---|---|---|---|---|---|",
    ])
    for _, r in tm_merged[tm_merged["category"].isin(["clinical", "common_trait"])].iterrows():
        if r["gene"] == "APOE":
            continue
        cat_lbl = "Clinical" if r["category"] == "clinical" else "Common Trait"
        md.append(
            f"| {cat_lbl} | **{r['gene']}** | {r['rsid']} | `{r[f'gt_{male_lbl}']}` | `{r[f'gt_{fem_lbl}']}` | {r['description']} |"
        )

    md.extend([
        "",
        "---",
        "",
        "## 4. Technical False Alarms Automatically Ruled Out",
        "",
    ])
    if filtered_artifacts:
        md.append("| Partner | Gene | Variant (rsID) | gnomAD v4 AF | Why Filtered Out |")
        md.append("|---|---|---|---|---|")
        for a in filtered_artifacts:
            md.append(
                f"| **{a['partner'].title()}** | **{a['gene'] or '.'}** | `{a['variant_id']}` ({a['rsid']}) | "
                f"{format_af(a.get('gnomad_af'), a.get('gnomad_status'))} | {a.get('artifact_reason', '')} |"
            )
    else:
        md.append("_None._")

    md.extend([
        "",
        "---",
        "",
        "## 5. BAM-Resolved Blind Spots & Remaining Next Steps",
        "",
    ])
    if has_bam:
        def _fmt_smn(b: dict[str, Any] | None) -> str:
            if not b or b.get("smn", {}).get("status") != "OK":
                return "Not run"
            s = b["smn"]
            return f"SMN1 CN={s['smn1_cn']}, SMN2 CN={s['smn2_cn']}, isCarrier={s['is_carrier']} ({s['interpretation']})"

        def _fmt_fmr1(b: dict[str, Any] | None) -> str:
            if not b or b.get("expansion_hunter", {}).get("status") != "OK":
                return "Not run"
            fmr1 = b["expansion_hunter"].get("fmr1")
            if not fmr1:
                return "Not called"
            return f"`{fmr1['genotype']}` CGG repeats (CI `{fmr1['ci']}` — **{fmr1['classification']}**)"

        def _fmt_repeats(b: dict[str, Any] | None) -> str:
            if not b or b.get("expansion_hunter", {}).get("status") != "OK":
                return "Not run"
            eh = b["expansion_hunter"]
            if eh["n_flagged"] == 0:
                return f"All {eh['n_loci']} loci normal (*HTT*, *C9ORF72*, *DMPK*, *FXN*, *ATXN*)"
            fl = ", ".join(f"{r['locus']} `{r['genotype']}` ({r['classification']})" for r in eh["flagged"])
            return f"{eh['n_flagged']}/{eh['n_loci']} flagged: {fl}"

        def _fmt_cyp2d6(b: dict[str, Any] | None) -> str:
            if not b or b.get("cyrius_cyp2d6", {}).get("status") != "OK":
                return "Not run"
            c = b["cyrius_cyp2d6"]
            gt = str(c.get("genotype") or "").strip()
            if not gt or gt.lower() == "none":
                raw_cn = c.get("total_cn_raw")
                cn_txt = f"Total_CN_raw={raw_cn}" if raw_cn is not None else "borderline CN"
                return f"Cyrius no-call ({cn_txt}; confirm exact CN/hybrids clinically)"
            return f"`{gt}` ({c['filter']}) — **{c['cpic_phenotype']}**"

        md.extend([
            "The following former VCF blind spots have now been directly screened using the raw 30x WGS BAM alignments (research-grade; confirm clinically before reproductive or prescribing decisions):",
            "",
            f"| Blind-Spot Assay | {m_name} | {f_name} |",
            "|---|---|---|",
            f"| **Fragile X (*FMR1* CGG Repeats)** | {_fmt_fmr1(male_bam)} | {_fmt_fmr1(fem_bam)} |",
            f"| **Other Repeat Expansions (31 Loci)** | {_fmt_repeats(male_bam)} | {_fmt_repeats(fem_bam)} |",
            f"| **SMA (*SMN1*/*SMN2* Exon 7 & 2+0)** | {_fmt_smn(male_bam)} | {_fmt_smn(fem_bam)} |",
            f"| ***CYP2D6* Full Star Alleles / CNVs** | {_fmt_cyp2d6(male_bam)} | {_fmt_cyp2d6(fem_bam)} |",
            "",
            "### Still Unscreened by Short-Read 30x WGS (Recommend Clinical Carrier Panel)",
            "- **Alpha-Thalassemia (*HBA1/HBA2*):** Common 3.7 kb / 4.2 kb deletions — clinical test recommended for **both partners**.",
            "- **Congenital Adrenal Hyperplasia (*CYP21A2*) & Gaucher Disease (*GBA1*):** Pseudogene homology interferes with short-read calling — clinical panel recommended for **both partners**.",
            "- **Clinical Confirmation:** Even with BAM-based *FMR1*, *SMN1*, and *CYP2D6* callers, standard CLIA-certified carrier screening / PGx testing is recommended prior to reproductive or medication changes.",
        ])
    else:
        md.extend([
            "Short-read 30x WGS cannot reliably clear these five high-priority reproductive genes without clinical testing (or specialized BAM/CRAM callers once raw alignments are uploaded):",
            f"- **Fragile X (*FMR1*):** CGG repeat expansion — clinical test recommended for **{f_name}** (or ExpansionHunter once BAM/CRAM is available, with clinical confirmation).",
            "- **Spinal Muscular Atrophy (*SMN1*):** Exon 7 copy number & `2+0` silent carriers — clinical test recommended for **both partners** (or SMNCopyNumberCaller on BAM/CRAM).",
            "- **Alpha-Thalassemia (*HBA1/HBA2*):** Common 3.7 kb / 4.2 kb deletions — clinical test recommended for **both partners**.",
            "- **Congenital Adrenal Hyperplasia (*CYP21A2*) & Gaucher Disease (*GBA1*):** Pseudogene homology interferes with short-read calling — clinical panel recommended for **both partners**.",
            "- **Pharmacogenomics (*CYP2D6*):** Standard VCFs detect `*4` (`rs3892097`) but cannot detect whole-gene deletions/duplications or *CYP2D7* hybrids; once BAM/CRAM files are added, Cyrius/Aldy can resolve full *CYP2D6* star alleles.",
        ])

    md.extend([
        "",
        "### Annual Checkup Workflow",
        "To refresh this entire analysis once a year with the latest ClinVar & GenCC clinical databases, run:",
        "```bash",
        "bash scripts/run_annual_checkup.sh --refresh-refs",
        "```",
    ])
    return "\n".join(md) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Joint couple reproductive carrier screening.")
    parser.add_argument("--male", required=True, help="Male partner label from config/samples.tsv (e.g. partner_a)")
    parser.add_argument("--female", required=True, help="Female partner label from config/samples.tsv (e.g. partner_b)")
    parser.add_argument(
        "--from-existing-json",
        action="store_true",
        help="Re-render couple reports from existing couple_carrier_screen.json + fresh BAM blind-spot JSONs (when data/processed VCFs are not present).",
    )
    args = parser.parse_args()

    male_lbl, fem_lbl = args.male, args.female
    out_dir = REPO_ROOT / "results" / "reports" / "couple"
    out_dir.mkdir(parents=True, exist_ok=True)

    male_proc = REPO_ROOT / "data" / "processed" / male_lbl
    fem_proc = REPO_ROOT / "data" / "processed" / fem_lbl
    male_rep = REPO_ROOT / "results" / "reports" / male_lbl
    fem_rep = REPO_ROOT / "results" / "reports" / fem_lbl

    male_bam = _load_bam_summary(male_rep)
    fem_bam = _load_bam_summary(fem_rep)

    existing_json = out_dir / "couple_carrier_screen.json"
    if args.from_existing_json and existing_json.exists():
        print(f"[Couple] Loading existing couple JSON ({existing_json}) and merging BAM blind-spot findings...")
        prev = json.loads(existing_json.read_text(encoding="utf-8"))
        ref_meta = prev["reference_metadata"]
        import datetime as _dt
        ref_meta["run_date"] = _dt.date.today().isoformat()
        sex_qc = prev["sex_qc"]
        shared_clinical_matches = prev["shared_clinical_matches"]
        shared_sv_only_matches = prev["shared_sv_only_matches"]
        partner_rescue_rows = prev["partner_rescue"]
        maternal_x_snv = prev["maternal_x_snv"]
        maternal_x_sv = prev["maternal_x_sv"]
        maternal_mt_snv = prev["maternal_mt_snv"]
        ad_acmg_hits = prev["ad_acmg_hits"]
        novel_acmg_lofs = prev["novel_acmg_lofs"]
        filtered_artifacts = prev["filtered_artifacts"]
        tm_merged = pd.read_csv(out_dir / "couple_targeted_markers.tsv", sep="\t")
    else:
        ref_meta = get_reference_metadata(REFERENCE_DIR)
        api = CachedAnnotationClient(REFERENCE_DIR / "api_cache.json")

        print(f"[Couple] Running sex concordance QC ({male_lbl}=M, {fem_lbl}=F)...")
        sex_qc = {
            male_lbl: check_sex_qc(male_proc / "snp_indel.variants.norm.vcf.gz", "M"),
            fem_lbl: check_sex_qc(fem_proc / "snp_indel.variants.norm.vcf.gz", "F"),
        }

        print("[Couple] Loading per-partner clinical hits & exonic CNV/SV overlaps...")
        male_hits, male_artifacts, male_novel_lofs = collect_partner_clinical_hits(
            male_rep / "carrier_screen_enriched.json"
        )
        fem_hits, fem_artifacts, fem_novel_lofs = collect_partner_clinical_hits(
            fem_rep / "carrier_screen_enriched.json"
        )
        male_sv = pd.read_csv(male_rep / "cnv_sv_exonic_clinical.tsv", sep="\t")
        fem_sv = pd.read_csv(fem_rep / "cnv_sv_exonic_clinical.tsv", sep="\t")

        _, cds_by_gene = load_gencode_cds(REFERENCE_DIR / "gencode.v46.cds.nochr.bed.gz")

        # 1. Gene-level matching (explicitly separating SNV/indel matches from SV+SV-only overlaps!)
        male_by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
        fem_by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for h in male_hits:
            if h["gene"]:
                male_by_gene[h["gene"]].append(h)
        for h in fem_hits:
            if h["gene"]:
                fem_by_gene[h["gene"]].append(h)

        def sv_gene_map(df: pd.DataFrame) -> dict[str, list[dict[str, Any]]]:
            m: dict[str, list[dict[str, Any]]] = defaultdict(list)
            if df.empty:
                return m
            for _, r in df.iterrows():
                for g in str(r.get("cds_genes", "")).replace("...", "").split(","):
                    g = g.strip()
                    if g:
                        m[g].append(r.to_dict())
            return m

        male_sv_genes = sv_gene_map(male_sv)
        fem_sv_genes = sv_gene_map(fem_sv)

        all_shared_genes = sorted(
            (set(male_by_gene) | set(male_sv_genes)) & (set(fem_by_gene) | set(fem_sv_genes))
        )
        shared_clinical_matches: list[dict[str, Any]] = []
        shared_sv_only_matches: list[dict[str, Any]] = []
        for g in all_shared_genes:
            entry = {
                "gene": g,
                f"{male_lbl}_snv": male_by_gene.get(g, []),
                f"{fem_lbl}_snv": fem_by_gene.get(g, []),
                f"{male_lbl}_sv": male_sv_genes.get(g, []),
                f"{fem_lbl}_sv": fem_sv_genes.get(g, []),
            }
            if g in male_by_gene or g in fem_by_gene:
                shared_clinical_matches.append(entry)
            else:
                entry["is_known_polymorphic"] = g in KNOWN_POLYMORPHIC_SV_GENES
                entry["note"] = KNOWN_POLYMORPHIC_SV_GENES.get(
                    g, "Shared short-read exonic SV call; verify in BAM/IGV or population SV database (gnomAD-SV)"
                )
                shared_sv_only_matches.append(entry)

        # 2. Partner-Gene Rare/Novel Coding Rescue & Strict PASS+MinDP Reciprocal gVCF Coverage Audit
        print("[Couple] Scanning partner genomes for rare/novel coding variants & PASS gVCF coverage...")
        partner_rescue_rows: list[dict[str, Any]] = []
        carrier_gene_pairs: list[tuple[str, str, dict[str, Any]]] = []
        for h in male_hits:
            if h["gene"] and h["chrom"] not in ("X", "Y", "MT"):
                carrier_gene_pairs.append((male_lbl, fem_lbl, h))
        for h in fem_hits:
            if h["gene"] and h["chrom"] not in ("X", "Y", "MT"):
                carrier_gene_pairs.append((fem_lbl, male_lbl, h))

        checked_gene_partner: set[tuple[str, str]] = set()
        with pysam.VariantFile(str(male_proc / "snp_indel.genome.vcf.gz")) as male_gvf, \
             pysam.VariantFile(str(fem_proc / "snp_indel.genome.vcf.gz")) as fem_gvf:
            for carrier_lbl, partner_lbl, hit in carrier_gene_pairs:
                gene = hit["gene"]
                if (gene, partner_lbl) in checked_gene_partner:
                    continue
                checked_gene_partner.add((gene, partner_lbl))
                cds_entry = cds_by_gene.get(gene)
                if not cds_entry:
                    continue
                chrom, exons = cds_entry
                partner_proc = male_proc if partner_lbl == male_lbl else fem_proc
                partner_gvf = male_gvf if partner_lbl == male_lbl else fem_gvf
                cov = compute_gene_gvcf_coverage(partner_gvf, chrom, exons, min_pass_dp=10)
                rare_vars = scan_partner_gene_for_rare_coding_variants(
                    partner_proc / "snp_indel.variants.norm.vcf.gz", gene, chrom, exons, api
                )
                partner_svs = (male_sv_genes if partner_lbl == male_lbl else fem_sv_genes).get(gene, [])
                partner_rescue_rows.append({
                    "carrier_partner": carrier_lbl,
                    "checked_partner": partner_lbl,
                    "gene": gene,
                    "moi": hit["moi"],
                    "carrier_variant": hit["variant_id"],
                    "carrier_rsid": hit["rsid"],
                    "carrier_evidence": f"{hit['evidence_tier']} ({hit['clnsig']})",
                    "partner_cds_bp": cov["cds_bp"],
                    "partner_pct_ge_10x": cov["pct_ge_10x"],
                    "partner_rare_coding_count": len(rare_vars),
                    "partner_rare_coding_variants": rare_vars,
                    "partner_exonic_sv_count": len(partner_svs),
                })
        api.save()

        # 3. Maternal X-linked & Mitochondrial screening
        maternal_x_snv = [h for h in fem_hits if h["chrom"] == "X"]
        maternal_mt_snv = [h for h in fem_hits if h["chrom"] == "MT"]
        maternal_x_sv = fem_sv[fem_sv["chrom"].astype(str) == "X"].to_dict("records") if not fem_sv.empty else []

        # 4. Autosomal Dominant (AD / AD-AR) & ACMG SF v3.3 in either partner
        ad_acmg_hits = []
        for lbl, lst in ((male_lbl, male_hits), (fem_lbl, fem_hits)):
            for h in lst:
                if h.get("acmg_sf") or h["moi"] in ("AD", "AD/AR"):
                    ad_acmg_hits.append({"partner": lbl, **h})

        novel_acmg_lofs = [
            {"partner": male_lbl, **r} for r in male_novel_lofs
        ] + [
            {"partner": fem_lbl, **r} for r in fem_novel_lofs
        ]

        filtered_artifacts = [
            {"partner": male_lbl, **r} for r in male_artifacts
        ] + [
            {"partner": fem_lbl, **r} for r in fem_artifacts
        ]

        # 5. Side-by-side canonical markers
        male_tm = pd.read_csv(male_rep / "targeted_markers.tsv", sep="\t")
        fem_tm = pd.read_csv(fem_rep / "targeted_markers.tsv", sep="\t")
        tm_merged = male_tm[["category", "gene", "rsid", "locus", "description", "gt", "dp"]].merge(
            fem_tm[["rsid", "gt", "dp"]], on="rsid", suffixes=(f"_{male_lbl}", f"_{fem_lbl}")
        )

        pd.DataFrame(partner_rescue_rows).drop(columns=["partner_rare_coding_variants"], errors="ignore").to_csv(
            out_dir / "partner_gene_reciprocal_audit.tsv", sep="\t", index=False
        )
        tm_merged.to_csv(out_dir / "couple_targeted_markers.tsv", sep="\t", index=False)

    payload = {
        "reference_metadata": ref_meta,
        "sex_qc": sex_qc,
        "shared_clinical_matches": shared_clinical_matches,
        "shared_sv_only_matches": shared_sv_only_matches,
        "partner_rescue": partner_rescue_rows,
        "maternal_x_snv": maternal_x_snv,
        "maternal_x_sv": maternal_x_sv,
        "maternal_mt_snv": maternal_mt_snv,
        "ad_acmg_hits": ad_acmg_hits,
        "novel_acmg_lofs": novel_acmg_lofs,
        "filtered_artifacts": filtered_artifacts,
        "bam_blind_spots": {male_lbl: male_bam, fem_lbl: fem_bam},
    }
    with open(out_dir / "couple_carrier_screen.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    # Render technical Markdown Couple Report
    md = [
        f"# Couple Reproductive Carrier Screening Report ({male_lbl.title()} + {fem_lbl.title()})",
        "",
        f"**Run Date:** {ref_meta['run_date']} | **ClinVar:** {ref_meta['clinvar_date']} | "
        f"**GenCC:** {ref_meta['gencc_date']} | **ACMG SF:** {ref_meta['acmg_sf_version']} | **GENCODE:** {ref_meta['gencode_release']}",
        "",
        "> **Clinical Disclaimer:** This is a research-grade bioinformatics analysis for informational and educational exploration only. "
        "It is **not** a clinical reproductive carrier screening test. Before making reproductive decisions, any findings—and especially "
        "short-read WGS blind spots (*FMR1*, *SMN1*, *HBA1/HBA2*, *CYP21A2*, *GBA1*)—must be evaluated via a CLIA-certified clinical carrier "
        "screening panel and a board-certified genetic counselor.",
        "",
        "## 1. Sample & Sex-Chromosome Concordance QC",
        "| Partner | Expected Sex | Inferred Sex | Concordant | Non-PAR X Het | Non-PAR X Hom/Hemi | X Het Ratio | Y Variants |",
        "|---|---|---|---|---:|---:|---:|---:|",
    ]
    for lbl in (male_lbl, fem_lbl):
        q = sex_qc[lbl]
        md.append(
            f"| **{lbl.title()}** | {q['expected_sex']} | {q['inferred_sex']} | "
            f"{'PASS' if q['concordant'] else 'FAIL'} | {q['non_par_x_het']:,} | "
            f"{q['non_par_x_hom_hemi']:,} | {q['non_par_x_het_ratio']:.3f} | {q['y_variant_count']:,} |"
        )

    md.extend([
        "",
        "## 2. Shared Gene-Level Autosomal Recessive Matches (Potential 25% Offspring Risk)",
        "### 2a. Shared SNV/Indel or SNV/Indel + Exonic SV Matches",
        "_Genes where both partners carry a non-artifact ClinVar P/LP, strong Conflicting P/LP, or SNV + exonic SV call._",
        "",
    ])
    if shared_clinical_matches:
        for m in shared_clinical_matches:
            md.append(f"- **`{m['gene']}`**: {male_lbl}={len(m[f'{male_lbl}_snv'])} SNV / {len(m[f'{male_lbl}_sv'])} SV; {fem_lbl}={len(m[f'{fem_lbl}_snv'])} SNV / {len(m[f'{fem_lbl}_sv'])} SV")
    else:
        md.append("**None identified (0 shared AR carrier genes across non-artifact ClinVar P/LP, strong Conflicting P/LP, and SNV+SV pairs).**")

    md.extend([
        "",
        "### 2b. Shared Exonic CNV/SV-Only Overlaps (Transparency & BAM/IGV Audit)",
        "_Genes where both partners carry a short-read coding-exon CNV/SV call without an SNV/indel hit. Many represent common polymorphic SVs, HLA diversity, or segmental-duplication artifacts, but are listed transparently rather than silently dropped._",
        "",
    ])
    male_sv_bam = _sv_bam_verdict_by_gene(male_bam)
    fem_sv_bam = _sv_bam_verdict_by_gene(fem_bam)
    if shared_sv_only_matches:
        if male_bam or fem_bam:
            md.append(f"| Gene | {male_lbl.title()} Call(s) | {fem_lbl.title()} Call(s) | BAM Read-Depth Audit ({male_lbl.title()} / {fem_lbl.title()}) | Classification / Note |")
            md.append("|---|---|---|---|---|")
            for m in shared_sv_only_matches:
                m_sv_str = "; ".join(f"{s['svtype']} {s['length_bp']:,}bp ({s['gt_or_cn']})" for s in m[f"{male_lbl}_sv"])
                f_sv_str = "; ".join(f"{s['svtype']} {s['length_bp']:,}bp ({s['gt_or_cn']})" for s in m[f"{fem_lbl}_sv"])
                m_b = "; ".join(male_sv_bam.get(m["gene"], ["-"]))
                f_b = "; ".join(fem_sv_bam.get(m["gene"], ["-"]))
                tag = "Known polymorphic/complex locus" if m["is_known_polymorphic"] else "BAM audited"
                md.append(f"| **{m['gene']}** | {m_sv_str} | {f_sv_str} | {male_lbl.title()}: `{m_b}` / {fem_lbl.title()}: `{f_b}` | {tag}: {m['note']} |")
        else:
            md.append(f"| Gene | {male_lbl.title()} Call(s) | {fem_lbl.title()} Call(s) | Classification / Note |")
            md.append("|---|---|---|---|")
            for m in shared_sv_only_matches:
                m_sv_str = "; ".join(f"{s['svtype']} {s['length_bp']:,}bp ({s['gt_or_cn']})" for s in m[f"{male_lbl}_sv"])
                f_sv_str = "; ".join(f"{s['svtype']} {s['length_bp']:,}bp ({s['gt_or_cn']})" for s in m[f"{fem_lbl}_sv"])
                tag = "Known polymorphic/complex locus" if m["is_known_polymorphic"] else "Verify in BAM/gnomAD-SV"
                md.append(f"| **{m['gene']}** | {m_sv_str} | {f_sv_str} | {tag}: {m['note']} |")
    else:
        md.append("_None identified._")

    md.extend([
        "",
        "## 3. Reciprocal Partner-Gene Audit & Rare/Novel Coding Variant Rescue",
        "_For every non-X/MT gene where one partner carries a non-artifact ClinVar P/LP or strong Conflicting P/LP variant, we verified the other partner's strict gVCF coding coverage (`>=10x PASS`, conservative `MinDP`) AND scanned all GENCODE v46 CDS exons (±10 bp splice junctions) via Ensembl VEP + gnomAD v4 for any rare (`AF < 1%` or absent) `HIGH`/`MODERATE` impact variants._",
        "",
        "| Carrier | Gene | MOI | Carrier Variant (rsID) | Evidence | Partner Checked | Partner CDS ≥10x (`PASS`) | Partner Rare Coding/Splice Variants (`AF<1%`) | Partner Exonic SVs |",
        "|---|---|---|---|---|---|---:|---|---:|",
    ])
    for r in partner_rescue_rows:
        rv_str = "None"
        if r["partner_rare_coding_variants"]:
            parts = [
                f"`{rv['variant_id']}` ({rv['gt']}, {rv['impact']} {rv['consequence']}, "
                f"{rv.get('hgvsp') or rv.get('hgvsc') or ''}, gnomAD={format_af(rv.get('gnomad_af'), rv.get('gnomad_status'))})"
                for rv in r["partner_rare_coding_variants"]
            ]
            rv_str = "; ".join(parts)
        md.append(
            f"| **{r['carrier_partner'].title()}** | **{r['gene']}** | {r['moi']} | "
            f"`{r['carrier_variant']}` ({r['carrier_rsid']}) | {r['carrier_evidence']} | "
            f"{r['checked_partner'].title()} | **{r['partner_pct_ge_10x']:.2f}%** | {rv_str} | {r['partner_exonic_sv_count']} |"
        )

    md.extend([
        "",
        f"## 4. Maternal X-Linked & Mitochondrial Screening ({fem_lbl.title()})",
        "_Maternal X-linked pathogenic variants carry a **50% risk to male offspring** independent of paternal genotype; maternal `MT` variants are transmitted to **100% of offspring**._",
        "",
        "### 4a. Maternal X-Linked SNV/Indels (Non-Artifact ClinVar P/LP or Strong Conflicting)",
    ])
    if maternal_x_snv:
        md.append("| Gene | Variant | rsID | GT | QC (DP/GQ/AB) | ClinVar (Stars) | gnomAD AF | Condition |")
        md.append("|---|---|---|---|---|---|---|---|")
        for h in maternal_x_snv:
            md.append(
                f"| **{h['gene']}** | `{h['variant_id']}` | {h['rsid']} | `{h['gt']}` | "
                f"{h['qc_tier']} ({h['dp']}/{h['gq']}/{h['ab']}) | {h['clnsig']} ({h['clinvar_stars']}★) | "
                f"{format_af(h.get('gnomad_af'), h.get('gnomad_status'))} | {str(h.get('clndn','')).replace('|','; ').replace('_',' ')[:75]} |"
            )
    else:
        md.append("_None identified._")

    md.extend(["", "### 4b. Maternal X-Linked Coding-Exon CNVs/SVs"])
    if maternal_x_sv:
        for s in maternal_x_sv:
            md.append(
                f"- `{s['source']}` {s['svtype']} `{s['chrom']}:{s['start']:,}-{s['end']:,}` ({s['length_bp']:,} bp, GT=`{s['gt_or_cn']}`, {s['support']}) — CDS genes: **{s['cds_genes']}**"
            )
    else:
        md.append("_None identified._")

    md.extend(["", "### 4c. Maternal Mitochondrial (`MT`) ClinVar P/LP or Strong Conflicting Hits"])
    if maternal_mt_snv:
        for h in maternal_mt_snv:
            md.append(f"- **{h['gene']}** `{h['variant_id']}` ({h['rsid']}, GT=`{h['gt']}`, DP={h['dp']}, AB={h['ab']}, {h['clnsig']})")
    else:
        md.append("_None identified._")

    md.extend([
        "",
        f"## 5. Autosomal Dominant (AD) & ACMG SF {ACMG_SF_VERSION} Findings in Either Partner (50% Transmission Risk)",
    ])
    if ad_acmg_hits:
        md.append("| Partner | Gene | Variant | rsID | GT | QC (DP/GQ/AB/MQ) | ClinVar (Stars) | gnomAD v4 AF | Condition |")
        md.append("|---|---|---|---|---|---|---|---|---|")
        for h in ad_acmg_hits:
            md.append(
                f"| **{h['partner'].title()}** | **{h['gene']}** | `{h['variant_id']}` | {h['rsid']} | `{h['gt']}` | "
                f"{h['qc_tier']} ({h['dp']}/{h['gq']}/{h['ab']}/{h['mq']}) | {h['clnsig']} ({h['clinvar_stars']}★) | "
                f"{format_af(h.get('gnomad_af'), h.get('gnomad_status'))} | {str(h.get('clndn','')).replace('|','; ').replace('_',' ')[:75]} |"
            )
    else:
        md.append("_None identified._")

    md.extend(["", f"### 5b. Unreported Rare `HIGH`-Impact LoF Candidates in ACMG SF {ACMG_SF_VERSION} Genes"])
    if novel_acmg_lofs:
        md.append("| Partner | Gene | Variant | GT | DP/GQ/AB | Consequence | HGVS | gnomAD v4 AF |")
        md.append("|---|---|---|---|---|---|---|---|")
        for r in novel_acmg_lofs:
            md.append(
                f"| **{r['partner'].title()}** | **{r['gene']}** | `{r['variant_id']}` | `{r['gt']}` | "
                f"{r['dp']}/{r['gq']}/{r['ab']} | {r.get('consequence','')} | {r.get('hgvsp') or r.get('hgvsc') or '-'} | "
                f"{format_af(r.get('gnomad_af'), r.get('gnomad_status'))} |"
            )
    else:
        md.append("_None identified._")

    md.extend([
        "",
        "## 6. Filtered False Alarms & Common Polymorphisms (Audit Trail)",
        "_Automatically excluded from Sections 2–5 based on ACMG BA1/BS1 allele frequency rules, benign blood-group traits, or 0★ legacy entries._",
        "",
    ])
    if filtered_artifacts:
        md.append("| Partner | Gene | Variant | rsID | ClinVar (Stars) | gnomAD v4 AF | Filter Reason |")
        md.append("|---|---|---|---|---|---|---|")
        for a in filtered_artifacts:
            md.append(
                f"| **{a['partner'].title()}** | **{a['gene'] or '.'}** | `{a['variant_id']}` | {a['rsid']} | "
                f"{a['clnsig']} ({a['clinvar_stars']}★) | {format_af(a.get('gnomad_af'), a.get('gnomad_status'))} | {a.get('artifact_reason', '')} |"
            )
    else:
        md.append("_None._")

    md.extend([
        "",
        "## 7. Side-by-Side Canonical Clinical, PGx & Trait Markers",
        f"| Category | Gene | rsID | {male_lbl.title()} GT (DP) | {fem_lbl.title()} GT (DP) | Description |",
        "|---|---|---|---|---|---|",
    ])
    for _, r in tm_merged.iterrows():
        md.append(
            f"| {r['category']} | **{r['gene']}** | {r['rsid']} | `{r[f'gt_{male_lbl}']}` ({r[f'dp_{male_lbl}']}) | "
            f"`{r[f'gt_{fem_lbl}']}` ({r[f'dp_{fem_lbl}']}) | {r['description']} |"
        )

    if male_bam or fem_bam:
        def _smn_cell(b: dict[str, Any] | None) -> str:
            if not b or b.get("smn", {}).get("status") != "OK":
                return "Not run"
            s = b["smn"]
            return f"SMN1={s['smn1_cn']}, SMN2={s['smn2_cn']}, isCarrier={s['is_carrier']}, g.27134T>G={s['g27134tg']}"

        def _fmr1_cell(b: dict[str, Any] | None) -> str:
            if not b or b.get("expansion_hunter", {}).get("status") != "OK":
                return "Not run"
            fmr1 = b["expansion_hunter"].get("fmr1")
            return f"`{fmr1['genotype']}` (CI `{fmr1['ci']}`, {fmr1['classification']})" if fmr1 else "Not called"

        def _eh_cell(b: dict[str, Any] | None) -> str:
            if not b or b.get("expansion_hunter", {}).get("status") != "OK":
                return "Not run"
            eh = b["expansion_hunter"]
            if eh["n_flagged"] == 0:
                return f"0/{eh['n_loci']} expanded (all normal)"
            return "; ".join(f"{r['locus']} `{r['genotype']}` ({r['classification']})" for r in eh["flagged"])

        def _cyp_cell(b: dict[str, Any] | None) -> str:
            if not b or b.get("cyrius_cyp2d6", {}).get("status") != "OK":
                return "Not run"
            c = b["cyrius_cyp2d6"]
            gt = str(c.get("genotype") or "").strip()
            if not gt or gt.lower() == "none":
                raw_cn = c.get("total_cn_raw")
                cn_txt = f"Total_CN_raw={raw_cn}" if raw_cn is not None else "borderline CN"
                return f"No call ({cn_txt}; confirm CN clinically)"
            return f"`{gt}` ({c['filter']}; {c['cpic_phenotype']})"

        md.extend([
            "",
            "## 8. BAM-Based Blind-Spot Callers (ExpansionHunter, SMNCopyNumberCaller, Cyrius)",
            f"| Assay | {male_lbl.title()} | {fem_lbl.title()} |",
            "|---|---|---|",
            f"| **Fragile X (*FMR1* CGG Repeats)** | {_fmr1_cell(male_bam)} | {_fmr1_cell(fem_bam)} |",
            f"| **Other Repeat Expansions (31 GRCh38 Loci)** | {_eh_cell(male_bam)} | {_eh_cell(fem_bam)} |",
            f"| **SMA (*SMN1*/*SMN2* Exon 7 & 2+0)** | {_smn_cell(male_bam)} | {_smn_cell(fem_bam)} |",
            f"| ***CYP2D6* Star Alleles / CNVs / Hybrids** | {_cyp_cell(male_bam)} | {_cyp_cell(fem_bam)} |",
            "",
            "## 9. Remaining Reproductive Blind-Spot Checklist (Not Cleared by Short-Read 30x WGS)",
            "- **Both Partners *HBA1 / HBA2* (Alpha-Thalassemia):** 3.7 kb / 4.2 kb deletions.",
            "- **Both Partners *CYP21A2* (Congenital Adrenal Hyperplasia)** and ***GBA1* (Gaucher Disease):** Pseudogene homology interferes with short-read calling.",
            "- **Clinical Confirmation Note:** *FMR1*, *SMN1*, and *CYP2D6* were screened via research-grade BAM callers (`ExpansionHunter`, `SMNCopyNumberCaller`, `Cyrius`); confirm via a CLIA-certified clinical panel before reproductive or prescribing decisions.",
        ])
    else:
        md.extend([
            "",
            "## 8. Essential Reproductive Blind-Spot Checklist (Not Cleared by Short-Read 30x WGS)",
            "- **Maternal *FMR1* (Fragile X Syndrome):** CGG repeat expansion (invisible in short-read VCFs).",
            "- **Both Partners *SMN1* (Spinal Muscular Atrophy):** Exon 7 copy number & `2+0` silent carrier testing.",
            "- **Both Partners *HBA1 / HBA2* (Alpha-Thalassemia):** 3.7 kb / 4.2 kb deletions.",
            "- **Both Partners *CYP21A2* (Congenital Adrenal Hyperplasia)** and ***GBA1* (Gaucher Disease):** Pseudogene homology interferes with short-read calling.",
        ])

    report_path = out_dir / "couple_carrier_screening_report.md"
    report_path.write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"[Couple] Technical report written to {report_path}")

    plain_md = render_plain_language_summary(
        male_lbl=male_lbl,
        fem_lbl=fem_lbl,
        ref_meta=ref_meta,
        shared_snv_matches=shared_clinical_matches,
        shared_sv_only=shared_sv_only_matches,
        partner_rescue_rows=partner_rescue_rows,
        maternal_x_snv=maternal_x_snv,
        maternal_x_sv=maternal_x_sv,
        maternal_mt_snv=maternal_mt_snv,
        ad_acmg_hits=ad_acmg_hits,
        novel_acmg_lofs=novel_acmg_lofs,
        filtered_artifacts=filtered_artifacts,
        tm_merged=tm_merged,
        male_bam=male_bam,
        fem_bam=fem_bam,
    )
    plain_path = out_dir / "plain_language_summary.md"
    plain_path.write_text(plain_md, encoding="utf-8")
    print(f"[Couple] Auto-generated plain-language summary written to {plain_path}")


if __name__ == "__main__":
    main()
