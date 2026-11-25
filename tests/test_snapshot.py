"""Unit tests for Step 13 annual snapshots and the 'what changed' diff (synthetic snapshots and VCFs only)."""

from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

SCRIPT_DIR = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.snapshot import (
    SCHEMA_VERSION,
    build_snapshot,
    class_direction,
    clinvar_class,
    diff_snapshots,
    find_previous_snapshot,
    is_tracked,
    load_snapshot,
    render_changes_report,
    variant_entry,
    write_snapshot,
    zygosity_from_gt,
)

_spec = importlib.util.spec_from_file_location("step13", SCRIPT_DIR / "13_annual_snapshot.py")
step13 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(step13)


@pytest.mark.parametrize("clnsig,expected", [
    ("Pathogenic", "Pathogenic"),
    ("Likely_pathogenic", "Likely pathogenic"),
    ("Pathogenic/Likely_pathogenic", "Pathogenic/Likely pathogenic"),
    ("Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance", "Pathogenic/Likely pathogenic"),
    ("Likely_pathogenic,_low_penetrance", "Likely pathogenic"),
    ("Pathogenic|risk_factor", "Pathogenic"),
    ("Uncertain_significance", "VUS"),
    ("Uncertain_significance|drug_response", "VUS"),
    ("Conflicting_classifications_of_pathogenicity", "Conflicting"),
    ("Conflicting_interpretations_of_pathogenicity|other", "Conflicting"),
    ("Benign/Likely_benign", "Benign/Likely benign"),
    ("Likely_benign", "Likely benign"),
    ("Benign", "Benign"),
    ("drug_response", "Other: drug_response"),
    ("", ""),
    (".", ""),
])
def test_clinvar_class(clnsig, expected):
    assert clinvar_class(clnsig) == expected


def test_tracked_classes_and_direction():
    assert all(is_tracked(c) for c in ("VUS", "Conflicting", "Likely pathogenic", "Pathogenic"))
    assert not any(is_tracked(c) for c in ("Benign", "Likely benign", "Other: risk_factor", ""))
    assert class_direction("VUS", "Likely pathogenic") == "up"
    assert class_direction("Pathogenic", "Likely benign") == "down"
    assert class_direction("VUS", "Conflicting") == "lateral"
    assert class_direction("VUS", "Other: drug_response") == "lateral"


@pytest.mark.parametrize("gt,zyg", [("0/1", "heterozygous"), ("1|0", "heterozygous"), ("1/1", "homozygous"),
                                    ("1/2", "heterozygous"), ("1", "hemizygous"), ("0/0", "ref"), ("./.", "no call")])
def test_zygosity_from_gt(gt, zyg):
    assert zygosity_from_gt(gt) == zyg


def test_variant_entry_parses_stars_genes_and_condition():
    v = variant_entry("Uncertain_significance", "reviewed_by_expert_panel", "999001", "GENEA:111|GENEB:222",
                      "Synthetic_disorder_type_1", "0/1", "PASS", "31", "99")
    assert v["class"] == "VUS" and v["stars"] == 3
    assert v["gene"] == "GENEA,GENEB"
    assert v["condition"] == "Synthetic disorder type 1"
    assert v["review_status"] == "reviewed by expert panel"
    assert v["zygosity"] == "heterozygous"


def _var(cls_clnsig: str, stars_rev: str = "criteria_provided,_single_submitter", gene: str = "GENEA") -> dict:
    return variant_entry(cls_clnsig, stars_rev, "900000", f"{gene}:1", "Synthetic_condition", "0/1", "PASS", "30", "99")


def _snap(run_date: str, clinvar_date: str, variants: dict, card=None, hla=None, pc=None, lof=None) -> dict:
    return build_snapshot(
        "synthetic", run_date,
        {"clinvar_date": clinvar_date, "gencc_date": clinvar_date, "pharmcat_version": "9.9.9", "git_commit": "abc1234"},
        variants, lof or [],
        card_rows=card or [], hla_rows=hla or [], pharmcat_rows=pc or [],
    )


CARD_OLD = [
    {"gene": "GENEP", "diplotype": "*90/*90", "phenotype": "Normal Metabolizer", "tier": "NORMAL"},
    {"gene": "GENEQ", "diplotype": "*90/*91", "phenotype": "Intermediate Metabolizer", "tier": "MODERATE"},
]
CARD_NEW = [
    {"gene": "GENEP", "diplotype": "*90/*92", "phenotype": "Poor Metabolizer", "tier": "ACTIONABLE"},
    {"gene": "GENEQ", "diplotype": "*90/*91", "phenotype": "Intermediate Metabolizer", "tier": "MODERATE"},
    {"gene": "GENER", "diplotype": "*90/*90", "phenotype": "Normal Metabolizer", "tier": "NORMAL"},
]
HLA_OLD = [{"trait": "Synthetic trait", "genotype": "X*90:01/X*91:01", "tier": "ABSENT", "result": "absent"}]
HLA_NEW = [{"trait": "Synthetic trait", "genotype": "X*90:01/X*91:01", "tier": "PRESENT", "result": "present"}]
PC_OLD = [{"gene": "GENEP", "pharmcat_diplotype": "*90/*90", "pharmcat_phenotype": "Normal Metabolizer", "status": "AGREE"}]
PC_NEW = [{"gene": "GENEP", "pharmcat_diplotype": "*90/*92", "pharmcat_phenotype": "Poor Metabolizer", "status": "AGREE"}]


@pytest.fixture
def pair():
    old = _snap("2025-11-20", "2025-11-15", {
        "1-1000-A-G": _var("Uncertain_significance"),                     # -> LP (up)
        "1-2000-C-T": _var("Pathogenic", gene="GENEB"),                    # -> VUS (down)
        "2-3000-G-A": _var("Uncertain_significance", gene="GENEC"),       # -> conflicting (lateral)
        "2-4000-T-C": _var("Likely_pathogenic", gene="GENED"),            # stars 1 -> 3
        "3-5000-A-T": _var("Uncertain_significance", gene="GENEE"),       # dropped: now likely benign
        "3-6000-G-C": _var("Uncertain_significance", gene="GENEF"),       # dropped: no longer in ClinVar
        "4-7000-C-G": _var("Conflicting_classifications_of_pathogenicity", gene="GENEG"),  # dropped: not called
    }, CARD_OLD, HLA_OLD, PC_OLD, lof=[{"variant_id": "5-8000-AT-A", "gene": "GENEL"}])
    new = _snap("2026-11-25", "2026-11-20", {
        "1-1000-A-G": _var("Likely_pathogenic"),
        "1-2000-C-T": _var("Uncertain_significance", gene="GENEB"),
        "2-3000-G-A": _var("Conflicting_classifications_of_pathogenicity", gene="GENEC"),
        "2-4000-T-C": _var("Likely_pathogenic", "reviewed_by_expert_panel", gene="GENED"),
        "6-9000-T-G": _var("Pathogenic", gene="GENEH"),                   # new finding
    }, CARD_NEW, HLA_NEW, PC_NEW, lof=[{"variant_id": "5-8100-G-GA", "gene": "GENEL"}])
    return old, new


def _lookup(key: str) -> dict:
    return {
        "3-5000-A-T": {"class": "Likely benign", "carried": True},
        "3-6000-G-C": {"class": "", "carried": True},
        "4-7000-C-G": {"class": "Conflicting", "carried": False},
    }[key]


def test_diff_classifies_every_change(pair):
    old, new = pair
    rows = diff_snapshots(old, new, _lookup)
    by_key = {(r["section"], r["key"]): r for r in rows}

    rec = [r for r in rows if r["section"] == "reclassification"]
    assert [r["key"] for r in rec] == ["1-1000-A-G", "2-3000-G-A", "1-2000-C-T", "3-5000-A-T"]  # up, lateral, down
    assert [r["direction"] for r in rec] == ["up", "lateral", "down", "down"]
    assert by_key[("reclassification", "1-1000-A-G")]["old"] == "VUS (1 star)"
    assert by_key[("reclassification", "1-1000-A-G")]["new"] == "Likely pathogenic (1 star)"
    assert by_key[("reclassification", "3-5000-A-T")]["new"] == "Likely benign (now untracked)"

    assert by_key[("review_status", "2-4000-T-C")]["direction"] == "up"
    assert "Not tracked last run" in by_key[("new_finding", "6-9000-T-G")]["note"]
    assert "no longer in ClinVar" in by_key[("dropped_finding", "3-6000-G-C")]["note"]
    assert "Not called" in by_key[("dropped_finding", "4-7000-C-G")]["note"]

    assert by_key[("lof_rescue", "5-8100-G-GA")]["note"] == "new candidate"
    assert by_key[("lof_rescue", "5-8000-AT-A")]["note"] == "no longer a candidate"

    card = {r["key"]: r for r in rows if r["section"] == "medication_card"}
    assert set(card) == {"GENEP", "GENER"}  # GENEQ unchanged
    assert card["GENEP"]["note"] == "changed: diplotype, phenotype, tier"
    assert card["GENER"]["note"] == "added"
    assert by_key[("hla", "Synthetic trait")]["note"] == "changed: tier"
    assert by_key[("pharmcat", "GENEP")]["new"] == "*90/*92 / Poor Metabolizer"

    sections = [r["section"] for r in rows]
    assert sections == sorted(sections, key=["reclassification", "new_finding", "dropped_finding", "review_status",
                                             "lof_rescue", "medication_card", "hla", "pharmcat"].index)


def test_diff_without_lookup_and_identical_snapshots(pair):
    old, new = pair
    dropped = [r for r in diff_snapshots(old, new) if r["section"] == "dropped_finding"]
    assert len(dropped) == 3 and all(r["note"] == "Not in this snapshot." for r in dropped)
    assert diff_snapshots(new, json.loads(json.dumps(new))) == []


def test_render_report_sections_and_first_run(pair):
    old, new = pair
    md = render_changes_report("synthetic", old, new, diff_snapshots(old, new, _lookup))
    assert md.startswith("# What Changed Since the Last Annual Run — Synthetic")
    assert "| ClinVar release | 2025-11-15 | 2026-11-20 |" in md
    assert "4 reclassification(s) (1 toward pathogenic)" in md
    assert "## Review-status changes (same class)" in md
    assert "not refreshed" not in md
    first = render_changes_report("synthetic", None, new, [])
    assert "first snapshot" in first and "5 tracked variants" in first


def test_render_report_warns_when_clinvar_not_refreshed_and_escapes_pipes(pair):
    old, new = pair
    new["meta"]["clinvar_date"] = old["meta"]["clinvar_date"]
    new["variants"]["1-1000-A-G"]["condition"] = "Type A|Type B"
    md = render_changes_report("synthetic", old, new, diff_snapshots(old, new, _lookup))
    assert "ClinVar was not refreshed" in md
    assert "Type A\\|Type B" in md


def test_snapshot_round_trip_and_schema_check(tmp_path, pair):
    _, new = pair
    p = tmp_path / "snaps" / "2026-11-25.json"
    write_snapshot(p, new)
    assert b"\r\n" not in p.read_bytes()
    assert load_snapshot(p) == new
    bad = dict(new, schema_version=SCHEMA_VERSION + 1)
    write_snapshot(tmp_path / "bad.json", bad)
    with pytest.raises(ValueError, match="schema"):
        load_snapshot(tmp_path / "bad.json")


def test_find_previous_snapshot(tmp_path):
    for name in ("2024-11-01.json", "2025-11-20.json", "2026-11-25.json", "notes.json", "2027-01-01.json"):
        (tmp_path / name).write_text("{}", encoding="utf-8")
    assert find_previous_snapshot(tmp_path, "2026-11-25").name == "2025-11-20.json"
    assert find_previous_snapshot(tmp_path, "2024-11-01") is None
    assert find_previous_snapshot(tmp_path / "missing", "2026-11-25") is None


def test_cli_diff_mode_writes_md_and_tsv(tmp_path, pair):
    old, new = pair
    write_snapshot(tmp_path / "old.json", old)
    write_snapshot(tmp_path / "new.json", new)
    out = tmp_path / "report"
    proc = subprocess.run([sys.executable, str(SCRIPT_DIR / "13_annual_snapshot.py"), "--sample", "synthetic",
                           "--old", str(tmp_path / "old.json"), "--new", str(tmp_path / "new.json"),
                           "--out-dir", str(out)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    md = (out / "changes_since_last_run.md").read_bytes()
    tsv = (out / "changes_since_last_run.tsv").read_text(encoding="utf-8").splitlines()
    assert b"\r\n" not in md and b"New findings" in md
    assert tsv[0] == "section\tkey\tgene\told\tnew\tdirection\tnote"
    assert len(tsv) == 1 + len(diff_snapshots(old, new))


# --- synthetic VCFs for the ClinVar lookup -------------------------------------------------------------------------

CLINVAR_VCF = """##fileformat=VCFv4.2
##fileDate=2026-11-20
##contig=<ID=1,length=100000>
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="sig">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="rev">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="gene">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="dn">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t1000\t900001\tA\tG\t.\t.\tCLNSIG=Uncertain_significance;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=GENEA:1;CLNDN=Synthetic_condition
1\t2000\t900002\tC\tT\t.\t.\tCLNSIG=Likely_benign;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=GENEB:2;CLNDN=not_provided
1\t3000\t900003\tG\tA\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;GENEINFO=GENEC:3;CLNDN=Synthetic_disorder
1\t4000\t900004\tT\tC\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;GENEINFO=GENED:4;CLNDN=Synthetic_disorder
"""

SAMPLE_VCF = """##fileformat=VCFv4.2
##contig=<ID=1,length=100000>
##FILTER=<ID=PASS,Description="pass">
##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="dp">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="gq">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t1000\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99
1\t2000\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ\t1/1:28:90
1\t3000\t.\tG\tC\t50\tPASS\t.\tGT:DP:GQ\t0/1:25:80
1\t4000\t.\tT\tC\t50\tPASS\t.\tGT:DP:GQ\t0/0:33:99
"""


def _bgzip_vcf(tmp_path: Path, name: str, text: str) -> Path:
    plain = tmp_path / name
    plain.write_text(text, encoding="utf-8")
    pysam.tabix_index(str(plain), preset="vcf", force=True)
    return tmp_path / f"{name}.gz"


def test_make_lookup_reports_current_class_and_carrier_status(tmp_path):
    clinvar = _bgzip_vcf(tmp_path, "clinvar.vcf", CLINVAR_VCF)
    sample = _bgzip_vcf(tmp_path, "sample.vcf", SAMPLE_VCF)
    lookup = step13.make_lookup(sample, clinvar)
    assert lookup("1-1000-A-G") == {"class": "VUS", "carried": True}
    assert lookup("1-2000-C-T") == {"class": "Likely benign", "carried": True}
    assert lookup("1-3000-G-A") == {"class": "Pathogenic", "carried": False}  # sample has G>C, not G>A
    assert lookup("1-4000-T-C") == {"class": "Pathogenic", "carried": False}  # 0/0
    assert lookup("1-5000-A-C") == {"class": "", "carried": False}


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not installed")
def test_extract_tracked_variants_keeps_carried_plp_and_vus_only(tmp_path):
    clinvar = _bgzip_vcf(tmp_path, "clinvar.vcf", CLINVAR_VCF)
    sample = _bgzip_vcf(tmp_path, "sample.vcf", SAMPLE_VCF)
    got = step13.extract_tracked_variants(sample, clinvar, tmp_path)
    assert list(got) == ["1-1000-A-G"]  # LB not tracked, 1:3000 allele mismatch, 1:4000 not carried
    v = got["1-1000-A-G"]
    assert (v["class"], v["clinvar_id"], v["gene"], v["dp"], v["zygosity"]) == ("VUS", "900001", "GENEA", "30", "heterozygous")
