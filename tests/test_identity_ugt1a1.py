"""Unit tests for Step 3 sample identity QC and the Step 10 UGT1A1 TA-repeat call (synthetic data only)."""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import sys
from pathlib import Path

import pysam
import pytest

SCRIPT_DIR = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.identity import (
    bam_results_usable,
    bam_verified,
    classify_identity,
    compare_genotypes,
    coverage_ratios,
    genotype_bam_at_sites,
    load_bam_results,
    overall_status,
    render_identity_report,
    select_identity_sites,
    sex_from_coverage,
    sex_from_vcf,
    sex_status,
    vcf_chrx_metrics,
)
from lib.pgx import interpret_pgx_profile, render_pgx_medication_card
from lib.pharmcat import compare_card_with_pharmcat, describe_input_at
from lib.ugt1a1 import (
    LEFT_ANCHOR,
    RIGHT_ANCHOR,
    UGT1A1_TA_POS,
    call_ta_diplotype,
    count_ta_repeats,
    interpret_ugt1a1,
    ta_units,
)


def _load_script_module(filename: str, mod_name: str):
    spec = importlib.util.spec_from_file_location(mod_name, SCRIPT_DIR / filename)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _vcf(path: Path, body: list[str], contigs: dict[str, int], fmt: str = "GT:DP:GQ") -> Path:
    head = ["##fileformat=VCFv4.2", '##FILTER=<ID=PASS,Description="All filters passed">',
            '##FILTER=<ID=LowQual,Description="Low quality">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">']
    head += [f"##contig=<ID={c},length={n}>" for c, n in contigs.items()]
    head.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1")
    raw = path.with_suffix("")
    raw.write_text("\n".join(head + body) + "\n", encoding="utf-8")
    pysam.tabix_compress(str(raw), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    return path


# ------------------------------------------------------------------------------------------ UGT1A1


def _ta_read(n: int, left_pad: str = "CAGTCAAACATTAACTTGGTGTATCGATT", right_pad: str = "AGAGGGCGAACCTCTGGCAGGAGCAAAGG") -> str:
    return left_pad + LEFT_ANCHOR + "AT" * (n + 1) + RIGHT_ANCHOR + right_pad


def test_ta_units_from_read_sequence() -> None:
    assert ta_units(_ta_read(6)) == 6
    assert ta_units(_ta_read(7)) == 7
    assert ta_units(_ta_read(5)) == 5
    assert ta_units(_ta_read(6)[:45]) is None  # read ends inside the repeat: no right anchor
    assert ta_units(LEFT_ANCHOR + "ATATGTATAT" + RIGHT_ANCHOR) == "noncanonical"


def test_count_ta_repeats_from_bam(tmp_path: Path) -> None:
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": "2", "LN": 242_193_529}]})
    unsorted = tmp_path / "u.bam"
    reads = [("p1", 6), ("p1", 6), ("p2", 7), ("p3", 6), ("p4", 7), ("p5", 6), ("p5", 7), ("p6", None)]
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for i, (name, n) in enumerate(reads):
            a = pysam.AlignedSegment(header)
            seq = _ta_read(n) if n is not None else _ta_read(6)[:40]
            a.query_name, a.query_sequence = name, seq
            a.reference_id, a.reference_start = 0, UGT1A1_TA_POS - 40 + i
            a.mapping_quality, a.cigarstring = 60, f"{len(seq)}M"
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            if name in ("p1", "p5"):
                a.flag = 0x1 | (0x40 if i % 2 == 0 else 0x80)
            out.write(a)
    bam = tmp_path / "s.bam"
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    res = count_ta_repeats(bam)
    assert res["counts"] == {6: 2, 7: 2}  # p1 mates count once; p5 mates disagree
    assert res["mate_conflict_fragments"] == 1
    assert res["partial_reads"] == 1


def test_call_ta_diplotype() -> None:
    hom = call_ta_diplotype({6: 20})
    assert hom["status"] == "OK" and hom["alleles"] == [6, 6] and hom["lod"] > 2
    het = call_ta_diplotype({6: 10, 7: 9})
    assert het["status"] == "OK" and het["alleles"] == [6, 7]
    assert call_ta_diplotype({7: 15})["alleles"] == [7, 7]
    assert call_ta_diplotype({6: 5})["status"] == "LOW_DEPTH"
    unsure = call_ta_diplotype({6: 12, 7: 2})  # two TA7 reads: stutter or a real allele?
    assert unsure["status"] == "UNCERTAIN" and unsure["alleles"] is None
    assert sorted(map(tuple, unsure["candidates"])) == [(6, 6), (6, 7)]
    assert call_ta_diplotype({"6": 20})["alleles"] == [6, 6]  # keys come back as strings from JSON


def _ta(alleles: list[int], counts: dict[int, int] | None = None) -> dict:
    return {"status": "OK", "alleles": alleles, "spanning_fragments": 20, "counts": counts or {alleles[0]: 20}}


def test_interpret_ugt1a1() -> None:
    nm = interpret_ugt1a1(0, _ta([6, 6]))
    assert (nm["diplotype"], nm["phenotype"], nm["tier"]) == ("*1/*1", "Normal Metabolizer", "NORMAL")
    im = interpret_ugt1a1(0, _ta([6, 7]))
    assert (im["diplotype"], im["phenotype"], im["tier"]) == ("*1/*28", "Intermediate Metabolizer", "MODERATE")
    pm = interpret_ugt1a1(0, _ta([7, 7]))
    assert (pm["diplotype"], pm["tier"]) == ("*28/*28", "ACTIONABLE") and "70%" in pm["guidance"]
    assert interpret_ugt1a1(0, _ta([5, 6]))["phenotype"] == "Normal Metabolizer"  # *36 is normal function
    assert interpret_ugt1a1(0, _ta([6, 8]))["diplotype"] == "*1/*37"
    six_28 = interpret_ugt1a1(1, _ta([6, 7]))
    assert six_28["diplotype"] == "*6/*28" and six_28["phenotype"] == "Poor Metabolizer"
    assert interpret_ugt1a1(1, _ta([6, 6]))["diplotype"] == "*1/*6"
    assert interpret_ugt1a1(2, _ta([6, 7]))["diplotype"] == "*6/*28+*6"
    no6 = interpret_ugt1a1(None, _ta([6, 6]))
    assert "*6 unassessed" in no6["diplotype"] and no6["tier"] == "PARTIAL"
    odd = interpret_ugt1a1(0, _ta([6, 9]))
    assert odd["tier"] == "PARTIAL" and "TA9" in odd["diplotype"]


def test_pgx_card_uses_bam_ta_repeat() -> None:
    markers = [{"rsid": "rs4148323", "status": "REF_CONFIRMED", "gt": "0/0 (REF)"}]
    old = {e["gene"]: e for e in interpret_pgx_profile(markers)}["UGT1A1"]
    assert old["diplotype"] == "No *6 (*28 unassessed)"
    new = {e["gene"]: e for e in interpret_pgx_profile(markers, ugt1a1_ta=_ta([6, 7], {6: 11, 7: 9}))}["UGT1A1"]
    assert new["diplotype"] == "*1/*28" and new["source"] == "BAM_TA_REPEAT" and "TA6 x11, TA7 x9" in new["guidance"]
    unsure = {"status": "UNCERTAIN", "alleles": None, "spanning_fragments": 14, "counts": {6: 12, 7: 2}}
    partial = {e["gene"]: e for e in interpret_pgx_profile(markers, ugt1a1_ta=unsure)}["UGT1A1"]
    assert partial["tier"] == "PARTIAL" and "not decisive" in partial["guidance"]
    md = render_pgx_medication_card("t", {"run_date": "2026-11-25"}, interpret_pgx_profile(markers, ugt1a1_ta=_ta([6, 6])))
    assert "BAM read count" in md


def test_pharmcat_ugt1a1_comparison_and_input(tmp_path: Path) -> None:
    pc = {"UGT1A1": {"Source Diplotype": "*1/*1", "Phenotype": "Normal Metabolizer"}}
    card_im = {"gene": "UGT1A1", "diplotype": "*1/*28", "phenotype": "Intermediate Metabolizer", "tier": "MODERATE",
               "source": "BAM_TA_REPEAT"}
    ta = _ta([6, 7], {6: 10, 7: 10})
    row = compare_card_with_pharmcat([card_im], pc, {}, ta, "0/0 from a gVCF reference block")[0]
    assert row["status"] == "DISAGREE"
    assert "trust it over PharmCAT" in row["note"] and "does not measure the repeat length" in row["note"]
    card_nm = {**card_im, "diplotype": "*1/*1", "phenotype": "Normal Metabolizer", "tier": "NORMAL"}
    assert compare_card_with_pharmcat([card_nm], pc, {}, _ta([6, 6]), "0/0 from a gVCF reference block")[0]["status"] == "AGREE"
    card_partial = {"gene": "UGT1A1", "diplotype": "No *6 (*28 unassessed)", "phenotype": "Partial — *6 Negative", "tier": "PARTIAL"}
    unsure = {"status": "UNCERTAIN", "counts": {6: 12, 7: 2}, "spanning_fragments": 14}
    row = compare_card_with_pharmcat([card_partial], pc, {}, unsure, "no call")[0]
    assert row["status"] == "CARD_INCOMPLETE" and "made no call" in row["note"]

    vcf = _vcf(tmp_path / "in.vcf.gz", [f"2\t{UGT1A1_TA_POS}\t.\tCAT\tCATAT\t.\tPASS\t.\tGT\t0/0"], {"2": 242_193_529}, "GT")
    assert describe_input_at(vcf) == "0/0 from a gVCF reference block"
    vcf2 = _vcf(tmp_path / "in2.vcf.gz", [f"2\t{UGT1A1_TA_POS}\t.\tCAT\tCATAT\t.\tPASS\t.\tGT\t0/1"], {"2": 242_193_529}, "GT")
    assert describe_input_at(vcf2) == "a VCF variant call (0/1)"
    vcf3 = _vcf(tmp_path / "in3.vcf.gz", ["2\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1"], {"2": 242_193_529}, "GT")
    assert describe_input_at(vcf3) == "no call"


def test_bam_report_ugt1a1_section(tmp_path: Path) -> None:
    step10 = _load_script_module("10_bam_blind_spots.py", "step10_ugt_mod")
    ugt = {"status": "OK", "alleles": [6, 7], "diplotype": "TA6/TA7 (*1/*28)", "counts": {6: 9, 7: 8},
           "spanning_fragments": 17, "min_spanning": 10, "min_lod": 2.0, "lod": 4.1, "partial_reads": 3,
           "noncanonical_fragments": 0, "mate_conflict_fragments": 0, "candidate_labels": []}
    qc = {"status": "PASS", "summary": "BAM matches the VCF."}
    md = step10.render_bam_report("t", tmp_path / "t.bam", {}, {}, {}, [], [], None, ugt, qc)
    assert "## 8. *UGT1A1* Promoter TA Repeat" in md and "TA6/TA7 (*1/*28)" in md
    assert "Sample identity QC: **PASS**" in md
    assert "Unavailable" in "\n".join(step10.render_ugt1a1_section(None))


# ------------------------------------------------------------------------------------------ identity


def test_compare_genotypes_and_classify() -> None:
    sites = [{"chrom": "1", "pos": i, "ref": "A", "alt": "G", "gt": ("A", "G")} for i in range(1, 1001)]
    same = {("1", i): {"gt": ("A", "G"), "dp": 30} for i in range(1, 1001)}
    same[("1", 5)] = {"gt": ("A", "A"), "dp": 30}  # one discordant site
    same[("1", 6)] = {"gt": ("A", "G"), "dp": 3}  # too shallow to compare
    res = compare_genotypes(sites, same)
    assert (res["compared"], res["matching"], res["low_depth"]) == (999, 998, 1)
    assert classify_identity(res["discordance"], res["compared"]) == "MATCH"
    other = {("1", i): {"gt": ("A", "A") if i % 2 else ("A", "G"), "dp": 30} for i in range(1, 1001)}
    res = compare_genotypes(sites, other)
    assert res["discordance"] == 0.5 and classify_identity(res["discordance"], res["compared"]) == "MISMATCH"
    assert classify_identity(0.1, 1000) == "INCONCLUSIVE"
    assert classify_identity(0.0, 100) == "INCONCLUSIVE"  # too few sites


def test_select_identity_sites(tmp_path: Path) -> None:
    body = [
        # windows are [500, 900) and [1500, 1900)
        "1\t600\t.\tA\tG\t50\tLowQual\t.\tGT:DP:GQ\t0/1:30:99",  # filtered
        "1\t650\t.\tAT\tA\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99",  # indel
        "1\t700\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99",  # first good site in window 1
        "1\t800\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ\t1/1:30:99",
        "1\t1600\t.\tG\tA\t50\tPASS\t.\tGT:DP:GQ\t0/1:5:99",  # shallow
        "1\t1650\t.\tG\tA\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99",  # two records at one position (split multiallelic)
        "1\t1650\t.\tG\tC\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99",
        "1\t1680\t.\tG\tA\t50\tPASS\t.\tGT:DP:GQ\t1/1:30:99",
        "X\t200\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99",  # not autosomal
    ]
    vcf = _vcf(tmp_path / "v.vcf.gz", body, {"1": 3000, "X": 3000})
    sites = select_identity_sites(vcf, spacing=1000, window=400)
    assert [(s["pos"], s["gt"]) for s in sites] == [(700, ("C", "T")), (1680, ("A", "A"))]


def test_sex_inference(tmp_path: Path) -> None:
    autos = [(str(i), 1_000_000, 300_000) for i in range(1, 23)]
    male = coverage_ratios(autos + [("X", 1_000_000, 150_000), ("Y", 1_000_000, 120_000)])
    female = coverage_ratios(autos + [("X", 1_000_000, 300_000), ("Y", 1_000_000, 9_000)])
    assert sex_from_coverage(male["x_ratio"], male["y_ratio"]) == "M"
    assert sex_from_coverage(female["x_ratio"], female["y_ratio"]) == "F"
    assert sex_from_coverage(0.5, 0.0) == "UNCLEAR"  # one X, no Y: not a plain M or F

    hap = _vcf(tmp_path / "m.vcf.gz", ["X\t50000\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99",  # PAR1: ignored
                                       "X\t3000000\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ\t1:30:99",
                                       "X\t3000100\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ\t1:30:99"], {"X": 156_040_895})
    m = vcf_chrx_metrics(hap)
    assert (m["x_haploid"], m["x_het"]) == (2, 0) and sex_from_vcf(m) == "M"
    dip = _vcf(tmp_path / "f.vcf.gz", ["X\t3000000\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ\t0/1:30:99",
                                       "X\t3000100\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ\t1/1:30:99"], {"X": 156_040_895})
    assert sex_from_vcf(vcf_chrx_metrics(dip)) == "F"
    assert sex_status("F", "M") == "MISMATCH" and sex_status("M", "M") == "MATCH" and sex_status("", "M") == "UNCLEAR"

    assert overall_status("MATCH", ["MATCH", "MATCH"]) == "PASS"
    assert overall_status("MISMATCH", ["MATCH"]) == "FAIL"
    assert overall_status("INCONCLUSIVE", ["MATCH"]) == "FAIL"
    assert overall_status("NOT_CHECKED", ["MISMATCH"]) == "FAIL"
    assert overall_status("NOT_CHECKED", ["UNCLEAR"]) == "WARN"
    assert overall_status("NOT_CHECKED", ["MATCH"]) == "PASS"


def test_bam_results_blocked_after_identity_failure(tmp_path: Path) -> None:
    rep = tmp_path / "rep"
    rep.mkdir()
    (rep / "bam_blind_spots.json").write_text(json.dumps({"cyrius_cyp2d6": {"status": "OK"}}), encoding="utf-8")
    assert load_bam_results(rep) is not None  # no identity QC yet: earlier results stay usable
    bam = tmp_path / "s.bam"
    bam.write_bytes(b"x")
    ok, why = bam_verified(rep, bam)
    assert not ok and "no sample identity QC" in why

    qc = {"status": "PASS", "summary": "ok", "bam_check": {"status": "MATCH", "alignment_target": os.path.realpath(bam)}}
    (rep / "sample_identity_qc.json").write_text(json.dumps(qc), encoding="utf-8")
    assert bam_verified(rep, bam) == (True, "ok") and bam_results_usable(rep)
    other = tmp_path / "other.bam"
    other.write_bytes(b"y")
    assert bam_verified(rep, other)[0] is False  # link now points at a different file

    qc["status"] = "FAIL"
    qc["bam_check"]["status"] = "MISMATCH"
    (rep / "sample_identity_qc.json").write_text(json.dumps(qc), encoding="utf-8")
    assert load_bam_results(rep) is None and not bam_results_usable(rep)
    assert bam_verified(rep, bam)[0] is False

    md = render_identity_report("t", {**qc, "run_date": "2026-11-25",
                                      "sex_check": {"declared": "F", "vcf": {"inferred": "F", "status": "MATCH"}}})
    assert "**MISMATCH**" in md and "Declared sex" in md


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not installed")
def test_genotype_bam_at_sites_with_bcftools(tmp_path: Path) -> None:
    import random

    rng = random.Random(7)
    ref = "".join(rng.choice("ACGT") for _ in range(2000))
    fa = tmp_path / "ref.fa"
    fa.write_text(">1\n" + "\n".join(ref[i:i + 60] for i in range(0, len(ref), 60)) + "\n", encoding="utf-8")
    pysam.faidx(str(fa))
    het_pos, hom_pos = 800, 1200  # 1-based
    alt = {b: ("A" if b != "A" else "C") for b in "ACGT"}
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": len(ref)}]})
    unsorted = tmp_path / "u.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for i in range(60):
            start = 700 + i * 10
            seq = list(ref[start:start + 100])
            for pos, carries in ((het_pos, i % 2 == 0), (hom_pos, True)):
                if start < pos <= start + 100 and carries:
                    seq[pos - 1 - start] = alt[ref[pos - 1]]
            a = pysam.AlignedSegment(header)
            a.query_name, a.query_sequence = f"r{i}", "".join(seq)
            a.reference_id, a.reference_start, a.mapping_quality, a.cigarstring = 0, start, 60, "100M"
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            out.write(a)
    bam = tmp_path / "s.bam"
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    sites = [{"chrom": "1", "pos": het_pos, "gt": tuple(sorted((ref[het_pos - 1], alt[ref[het_pos - 1]])))},
             {"chrom": "1", "pos": hom_pos, "gt": (alt[ref[hom_pos - 1]],) * 2}]
    calls = genotype_bam_at_sites(bam, fa, sites, tmp_path / "work", has_chr=False)
    res = compare_genotypes(sites, calls)
    assert res["compared"] == 2 and res["matching"] == 2


def test_step10_drops_caller_cache_when_alignment_changes(tmp_path: Path) -> None:
    step10 = _load_script_module("10_bam_blind_spots.py", "step10_cache_mod")
    work = tmp_path / "bam_callers"
    bam_a, bam_b = tmp_path / "a.bam", tmp_path / "b.bam"
    bam_a.write_bytes(b"aaaa")
    bam_b.write_bytes(b"bbbbbb")
    step10.reset_caller_cache_if_alignment_changed(work, bam_a, "t")
    (work / "cyrius").mkdir()
    step10.reset_caller_cache_if_alignment_changed(work, bam_a, "t")
    assert (work / "cyrius").exists()  # same alignment: cache kept
    step10.reset_caller_cache_if_alignment_changed(work, bam_b, "t")
    assert not (work / "cyrius").exists() and (work / "alignment.json").exists()
