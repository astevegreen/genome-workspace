"""Unit tests for core genomics classification, artifact flagging, interval math, and API formatting."""

from __future__ import annotations

import datetime
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pysam

SCRIPT_DIR = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.api_client import CachedAnnotationClient, format_af
from lib.constants import ACMG_SF_GENES, TARGETED_MARKERS
from lib.pgx import interpret_pgx_profile, render_pgx_medication_card
from lib.genomics import (
    assess_artifact_or_common,
    classify_clnsig,
    compute_gene_gvcf_coverage,
    derive_apoe_diplotype,
    extract_genes,
    is_strong_conflicting,
    merge_intervals,
    normalize_gt_unphased,
    parse_revstat_stars,
    reciprocal_overlap,
    resolve_gene_moi,
)


def _load_script_module(filename: str, mod_name: str):
    spec = importlib.util.spec_from_file_location(mod_name, SCRIPT_DIR / filename)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_acmg_sf_v33_gene_count_and_additions() -> None:
    assert len(ACMG_SF_GENES) == 84
    for added in ("ABCD1", "CYP27A1", "PLN"):
        assert added in ACMG_SF_GENES


def test_classify_clnsig() -> None:
    assert classify_clnsig("Pathogenic") == "PLP"
    assert classify_clnsig("Pathogenic/Likely_pathogenic") == "PLP"
    assert classify_clnsig("Likely_pathogenic|risk_factor") == "PLP"
    assert classify_clnsig("Conflicting_classifications_of_pathogenicity", "Pathogenic(3)|Uncertain_significance(1)") == "CONFLICTING_PLP"
    assert classify_clnsig("Conflicting_classifications_of_pathogenicity", "Benign(3)|Uncertain_significance(1)") is None
    assert classify_clnsig("drug_response") == "RISK_OR_DRUG"
    assert classify_clnsig("Benign") is None
    assert classify_clnsig(".") is None


def test_parse_revstat_stars() -> None:
    assert parse_revstat_stars("practice_guideline")[1] == 4
    assert parse_revstat_stars("reviewed_by_expert_panel")[1] == 3
    assert parse_revstat_stars("criteria_provided,_multiple_submitters,_no_conflicts")[1] == 2
    assert parse_revstat_stars("criteria_provided,_single_submitter")[1] == 1
    assert parse_revstat_stars("no_assertion_criteria_provided")[1] == 0


def test_extract_genes() -> None:
    assert extract_genes("BRCA1:672|NBR2:10230") == ["BRCA1", "NBR2"]
    assert extract_genes("") == []


def test_is_strong_conflicting() -> None:
    assert is_strong_conflicting({"clnsigconf": "Pathogenic_(18),Uncertain_significance_(1)"}) is True
    assert is_strong_conflicting({"clnsigconf": "Pathogenic_(1),Benign_(5)"}) is False
    assert is_strong_conflicting({"clnsigconf": "Likely_pathogenic_(2),Likely_benign_(2)"}) is True
    assert is_strong_conflicting({"clnsigconf": ""}) is False


def test_assess_artifact_or_common() -> None:
    # Common BA1 polymorphism / sequencing artifact (>= 5% AF)
    flag, reason = assess_artifact_or_common({"gnomad_af": 0.20, "clinvar_stars": 2, "moi": "AD", "clndn": "Example autosomal dominant disorder"})
    assert flag == "LIKELY_ARTIFACT_OR_COMMON"
    assert "BA1" in reason

    # Common benign blood-group trait (0-star)
    flag, reason = assess_artifact_or_common({"gnomad_af": 0.30, "clinvar_stars": 0, "moi": "XL", "clndn": "EXAMPLE BLOOD GROUP SYSTEM"})
    assert flag == "LIKELY_ARTIFACT_OR_COMMON"
    assert "blood-group" in reason.lower()

    # Legacy 0-star somatic entry (> 1% AF)
    flag, reason = assess_artifact_or_common({"gnomad_af": 0.02, "clinvar_stars": 0, "moi": "UNKNOWN", "clndn": "Example carcinoma"})
    assert flag == "LIKELY_ARTIFACT_OR_COMMON"

    # True rare AR carrier (textbook example, 2 stars)
    flag, _ = assess_artifact_or_common({"gnomad_af": 0.001, "clinvar_stars": 2, "moi": "AR", "clndn": "Cystic fibrosis"})
    assert flag == "OK"


def test_interval_math() -> None:
    assert merge_intervals([(10, 20), (15, 30), (40, 50), (50, 60)]) == [(10, 30), (40, 60)]
    assert merge_intervals([]) == []
    assert reciprocal_overlap(0, 100, 0, 100) == 1.0
    assert reciprocal_overlap(0, 100, 50, 150) == 0.5
    assert reciprocal_overlap(0, 100, 200, 300) == 0.0


def test_format_af_distinguishes_absent_vs_error_vs_unqueried() -> None:
    assert format_af(None, "absent") == "absent"
    assert format_af(None, "api_error") == "API error"
    assert format_af(None, "not_applicable") == "N/A (MT)"
    assert format_af(None, None) == "not queried"
    assert format_af(0.0025, "ok") == "0.2500%"
    assert format_af(1e-5, "ok") == "1.00e-05"


def test_targeted_markers_unique_rsids() -> None:
    rsids = [m[2] for m in TARGETED_MARKERS]
    assert len(rsids) == len(set(rsids))
    assert len(TARGETED_MARKERS) >= 30


def test_resolve_gene_moi_multigene_and_chrom_fallback() -> None:
    gencc = {
        "REALGENE": {"moi": "AR", "gencc_diseases": "Test AR disease"},
        "OTHER": {"moi": "UNKNOWN", "gencc_diseases": ""},
    }
    gene, moi, dis = resolve_gene_moi(["READTHROUGH", "REALGENE"], gencc, "1")
    assert gene == "REALGENE"
    assert moi == "AR"
    assert dis == "Test AR disease"

    gene_x, moi_x, _ = resolve_gene_moi(["UNCURATED_X"], gencc, "X")
    assert gene_x == "UNCURATED_X"
    assert moi_x == "XL"

    gene_mt, moi_mt, _ = resolve_gene_moi(["MT-ND1"], gencc, "MT")
    assert gene_mt == "MT-ND1"
    assert moi_mt == "MT"


def test_derive_apoe_diplotype_phased_and_unphased() -> None:
    assert normalize_gt_unphased("1|0") == "0/1"
    assert normalize_gt_unphased("1|1") == "1/1"
    assert derive_apoe_diplotype("0/0 (REF)", "0/0 (REF)").startswith("ε3/ε3")
    assert derive_apoe_diplotype("1|0", "0/0 (REF)") == "ε3/ε4"
    assert derive_apoe_diplotype("0/0", "0|1") == "ε2/ε3"
    assert derive_apoe_diplotype("1|1", "0|0") == "ε4/ε4"
    assert derive_apoe_diplotype("./.", "0/0 (REF)").startswith("Unresolved")


def test_api_cache_ttl_and_refresh(tmp_path: Path) -> None:
    cache_file = tmp_path / "api_cache.json"
    old_date = (datetime.datetime.now(datetime.timezone.utc) - datetime.timedelta(days=400)).isoformat()
    cache_file.write_text(
        json.dumps({
            "_meta": {"created_at": old_date, "updated_at": old_date},
            "gnomad_variants": {"1-100-A-G": {"gnomad_af": 0.001, "gnomad_status": "ok"}},
            "gnomad_genes": {},
            "ensembl_vep": {},
        }),
        encoding="utf-8",
    )
    # Expired (>365 days) -> ignored
    expired_client = CachedAnnotationClient(cache_file, max_age_days=365)
    assert expired_client._cache["gnomad_variants"] == {}

    # Non-expired when max_age_days=500 -> loaded
    fresh_client = CachedAnnotationClient(cache_file, max_age_days=500)
    assert "1-100-A-G" in fresh_client._cache["gnomad_variants"]

    # Explicit refresh=True -> cleared
    refreshed = CachedAnnotationClient(cache_file, max_age_days=500, refresh=True)
    assert refreshed._cache["gnomad_variants"] == {}


def test_evaluate_personal_call_and_stratify_hits() -> None:
    step4 = _load_script_module("04_carrier_screen.py", "step4_mod")

    class FakeFilter:
        def __init__(self, keys: list[str]) -> None:
            self._keys = keys

        def __len__(self) -> int:
            return len(self._keys)

        def keys(self) -> list[str]:
            return self._keys

    def make_rec(gt, dp, gq, ad, flt=None):
        return SimpleNamespace(
            samples={"S1": {"GT": gt, "DP": dp, "GQ": gq, "AD": ad}},
            filter=FakeFilter(flt or []),
            info={"MQ": 60},
            qual=100.0,
        )

    high = step4.evaluate_personal_call(make_rec((0, 1), 30, 99, (15, 15)))
    assert high["zygosity"] == "HET"
    assert high["qc_tier"] == "HIGH"

    hemi = step4.evaluate_personal_call(make_rec((1,), 8, 50, (0, 8)))
    assert hemi["zygosity"] == "HEMI_ALT"
    assert hemi["qc_tier"] == "HIGH"

    low_ab = step4.evaluate_personal_call(make_rec((0, 1), 30, 99, (26, 4)))
    assert low_ab["qc_tier"] == "LOW_CONFIDENCE"

    hits = [
        {
            "clinvar_category": "PLP", "clinvar_stars": 2, "qc_tier": "HIGH",
            "gene": "CFTR", "chrom": "7", "pos": 100, "zygosity": "HET",
            "acmg_sf": False, "moi": "AR",
        },
        {
            "clinvar_category": "PLP", "clinvar_stars": 2, "qc_tier": "HIGH",
            "gene": "CFTR", "chrom": "7", "pos": 200, "zygosity": "HET",
            "acmg_sf": False, "moi": "AR",
        },
        {
            "clinvar_category": "PLP", "clinvar_stars": 3, "qc_tier": "HIGH",
            "gene": "BRCA1", "chrom": "17", "pos": 300, "zygosity": "HET",
            "acmg_sf": True, "moi": "AD",
        },
    ]
    buckets = step4.stratify_hits(hits)
    assert len(buckets["hom_hemi_comphet"]) == 2
    assert all(r["potential_comphet"] for r in buckets["hom_hemi_comphet"])
    assert len(buckets["incidental_ad_acmg"]) == 1
    assert buckets["incidental_ad_acmg"][0]["gene"] == "BRCA1"


def test_compute_gene_gvcf_coverage_and_sex_qc_synthetic_vcf(tmp_path: Path) -> None:
    # Build a tiny bgzipped + tabix-indexed synthetic gVCF
    gvcf_raw = tmp_path / "test.gvcf"
    gvcf_raw.write_text(
        "\n".join([
            "##fileformat=VCFv4.2",
            '##FILTER=<ID=PASS,Description="All filters passed">',
            '##FILTER=<ID=LowQual,Description="Low quality">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
            '##INFO=<ID=MinDP,Number=1,Type=Integer,Description="Minimum DP">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
            '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP">',
            "##contig=<ID=1,length=10000>",
            "##contig=<ID=X,length=160000000>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
            # Bases 1..50 (0-based 0..50): PASS, MinDP=20 -> >=10x PASS
            "1\t1\t.\tA\t.\t.\tPASS\tEND=50;MinDP=20\tGT:DP\t0/0:25",
            # Bases 51..75 (0-based 50..75): LowQual, DP=20 -> >=1x only, NOT >=10x PASS
            "1\t51\t.\tC\t.\t.\tLowQual\tEND=75\tGT:DP:MIN_DP\t0/0:20:20",
            # Bases 76..90 (0-based 75..90): PASS, MIN_DP=5 -> >=1x only, NOT >=10x PASS
            "1\t76\t.\tG\t.\t.\tPASS\tEND=90\tGT:DP:MIN_DP\t0/0:15:5",
            # Bases 91..100: uncovered (0x)
            # Non-PAR X het variants for sex QC (female)
            "X\t10000000\t.\tA\tG\t99\tPASS\t.\tGT:DP\t0/1:30",
            "X\t20000000\t.\tC\tT\t99\tPASS\t.\tGT:DP\t0/1:30",
            "X\t30000000\t.\tG\tA\t99\tPASS\t.\tGT:DP\t1/1:30",
        ]) + "\n",
        encoding="utf-8",
    )
    gvcf_gz = tmp_path / "test.gvcf.gz"
    pysam.tabix_compress(str(gvcf_raw), str(gvcf_gz), force=True)
    pysam.tabix_index(str(gvcf_gz), preset="vcf", force=True)

    with pysam.VariantFile(str(gvcf_gz)) as vf:
        cov = compute_gene_gvcf_coverage(vf, "1", [(0, 100)], min_pass_dp=10)
    assert cov["cds_bp"] == 100
    assert cov["pct_ge_10x"] == 50.0
    assert cov["pct_ge_1x"] == 90.0
    assert cov["pct_0x_nocall"] == 10.0

    step8 = _load_script_module("08_couple_carrier_screen.py", "step8_mod")
    sex_res = step8.check_sex_qc(gvcf_gz, "F")
    assert sex_res["inferred_sex"] == "F"
    assert sex_res["concordant"] is True
    assert sex_res["non_par_x_het"] == 2
    assert sex_res["non_par_x_hom_hemi"] == 1


def test_cpic_pgx_medication_card_interpretation() -> None:
    markers = [
        # CYP2C19 *2/*17 -> Intermediate Metabolizer (ACTIONABLE)
        {"rsid": "rs4244285", "status": "VARIANT_PRESENT", "gt": "0/1"},
        {"rsid": "rs4986893", "status": "REF_CONFIRMED", "gt": "0/0 (REF)"},
        {"rsid": "rs28399504", "status": "REF_CONFIRMED", "gt": "0/0 (REF)"},
        {"rsid": "rs12248560", "status": "VARIANT_PRESENT", "gt": "0/1"},
        # SLCO1B1 *1/*5 -> Decreased Function (ACTIONABLE)
        {"rsid": "rs4149056", "status": "VARIANT_PRESENT", "gt": "0/1"},
        # DPYD c.2846A>T het -> AS=1.5 Intermediate Metabolizer (ACTIONABLE)
        {"rsid": "rs3918290", "status": "REF_CONFIRMED", "gt": "0/0 (REF)"},
        {"rsid": "rs55886062", "status": "REF_CONFIRMED", "gt": "0/0 (REF)"},
        {"rsid": "rs67376798", "status": "VARIANT_PRESENT", "gt": "0/1"},
        {"rsid": "rs75017182", "status": "REF_CONFIRMED", "gt": "0/0 (REF)"},
        # TPMT *3B + *3C het -> *1/*3A cis assumption (ACTIONABLE)
        {"rsid": "rs1800462", "status": "REF_CONFIRMED", "gt": "0/0 (REF)"},
        {"rsid": "rs1800460", "status": "VARIANT_PRESENT", "gt": "0/1"},
        {"rsid": "rs1142345", "status": "VARIANT_PRESENT", "gt": "0/1"},
        # CYP3A5 *3/*3 -> Non-Expressor (NORMAL standard tacrolimus dose)
        {"rsid": "rs776746", "status": "VARIANT_PRESENT", "gt": "1/1"},
        # CYP2D6 *4 negative -> PARTIAL
        {"rsid": "rs3892097", "status": "REF_CONFIRMED", "gt": "0/0 (REF)"},
    ]
    entries = {e["gene"]: e for e in interpret_pgx_profile(markers)}
    assert entries["CYP2C19"]["tier"] == "ACTIONABLE"
    assert "Intermediate" in entries["CYP2C19"]["phenotype"]
    assert entries["SLCO1B1"]["tier"] == "ACTIONABLE"
    assert "AS=1.5" in entries["DPYD"]["phenotype"]
    assert "*1/*3A" in entries["TPMT"]["diplotype"]
    assert entries["CYP3A5"]["diplotype"] == "*3/*3"
    assert entries["CYP2D6"]["tier"] == "PARTIAL"

    card_md = render_pgx_medication_card("test", {"run_date": "2026-11-20"}, list(entries.values()))
    assert "Personal Pharmacogenomics (PGx) Medication Card" in card_md
    assert "[ACTIONABLE]" in card_md


def test_prs_weight_orientation_and_exact_scoring(tmp_path: Path) -> None:
    step9 = _load_script_module("09_polygenic_risk.py", "step9_mod")

    # Effect allele == ALT: dosage(ALT) * w
    assert step9.orient_weight_to_alt("G", "A", 0.5, "A", "G") == (0.5, 0.0)
    # Effect allele == REF: 2*w + dosage(ALT)*(-w)
    assert step9.orient_weight_to_alt("A", "G", 0.5, "A", "G") == (-0.5, 1.0)
    # Incompatible alleles -> None
    assert step9.orient_weight_to_alt("T", "C", 0.5, "A", "G") is None

    # Empirical percentile
    assert step9.compute_empirical_percentile(5.0, [1.0, 2.0, 3.0, 4.0, 6.0]) == 80.0
    assert step9.risk_tier_from_percentile(96.0).startswith("High")

    # End-to-end exact raw score verification across 3 SNPs:
    # SNP1: 1:100 A>G, eff=G (ALT), w=+1.0, personal=0/1 -> contributes 1 * 1.0 = +1.0
    # SNP2: 1:200 C>T, eff=C (REF), w=+2.0, personal=0/1 -> contributes 1 * 2.0 = +2.0 (offset 4.0 + 1*(-2.0) = 2.0)
    # SNP3: 1:300 G>A, eff=G (REF), w=+0.5, personal=0/0 (absent from VCF) -> contributes 2 * 0.5 = +1.0
    # Expected exact raw score = 1.0 + 2.0 + 1.0 = 4.0
    vcf_raw = tmp_path / "p.vcf"
    vcf_raw.write_text(
        "\n".join([
            "##fileformat=VCFv4.2",
            '##FILTER=<ID=PASS,Description="All filters passed">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "##contig=<ID=1,length=10000>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
            "1\t100\t.\tA\tG\t99\tPASS\t.\tGT\t0/1",
            "1\t200\t.\tC\tT\t99\tPASS\t.\tGT\t0/1",
        ]) + "\n",
        encoding="utf-8",
    )
    vcf_gz = tmp_path / "p.vcf.gz"
    pysam.tabix_compress(str(vcf_raw), str(vcf_gz), force=True)
    pysam.tabix_index(str(vcf_gz), preset="vcf", force=True)

    score_tsv = tmp_path / "scores.tsv"
    score_tsv.write_text(
        "ID\tALT\tPGS_TEST\n"
        "1:100:A:G\tG\t1.0\n"
        "1:200:C:T\tT\t-2.0\n"
        "1:300:G:A\tA\t-0.5\n",
        encoding="utf-8",
    )
    dists = {"scores": {"PGS_TEST": {"ref_offset": 4.0 + 1.0}}}
    raw_scores, hits = step9.score_personal_wgs(vcf_gz, score_tsv, dists)
    assert abs(raw_scores["PGS_TEST"] - 4.0) < 1e-9
    assert hits["PGS_TEST"] == 2
