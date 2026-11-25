#!/usr/bin/env python3
"""
Step 3: Sample identity QC.
  1. BAM vs VCF: re-genotypes ~3,000 of the VCF's own SNVs from the BAM (bcftools mpileup/call) and counts
     how many agree. Same person: < 5% discordant. Another person: >= 20%. On anything but a match, the BAM
     is also tried against every other sample's VCF in config/samples.tsv, to name a swap.
  2. Genetic sex: chrX/chrY depth from the BAM index and the chrX het/haploid pattern of the VCF, each
     against the sex in config/samples.tsv.
Writes results/reports/<label>/sample_identity_qc.{md,json} and results/qc/<label>/bam.idxstats.tsv.
Exits 1 on FAIL, so run_sample.sh can stop the BAM steps; Step 10 refuses to run without a MATCH.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import shutil
import sys
from pathlib import Path
from typing import Any

import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.identity import (
    QC_JSON,
    QC_MD,
    alignment_fingerprint,
    classify_identity,
    compare_genotypes,
    coverage_ratios,
    genotype_bam_at_sites,
    load_identity_qc,
    overall_status,
    render_identity_report,
    select_identity_sites,
    sex_from_coverage,
    sex_from_vcf,
    sex_status,
    vcf_chrx_metrics,
)

REPO_ROOT = SCRIPT_DIR.parent
SAMPLES_TSV = REPO_ROOT / "config" / "samples.tsv"


def read_samples() -> list[dict[str, str]]:
    with open(SAMPLES_TSV, encoding="utf-8") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r.get("label")]


def find_alignment(label: str) -> Path | None:
    d = REPO_ROOT / "data" / "alignments" / label
    for ext in ("*.bam", "*.cram"):
        found = sorted(d.glob(ext)) if d.exists() else []
        if found:
            return found[0]
    return None


def variants_vcf(label: str) -> Path:
    return REPO_ROOT / "data" / "processed" / label / "snp_indel.variants.norm.vcf.gz"


def check_bam_against(bam: Path, vcf: Path, fasta: Path, work: Path, has_chr: bool, bcftools: str) -> dict[str, Any]:
    sites = select_identity_sites(vcf)
    calls = genotype_bam_at_sites(bam, fasta, sites, work, has_chr, bcftools)
    res = compare_genotypes(sites, calls)
    return {**res, "status": classify_identity(res["discordance"], res["compared"])}


def main() -> None:
    parser = argparse.ArgumentParser(description="Check that the BAM matches the VCF and that genetic sex matches samples.tsv")
    parser.add_argument("--sample", required=True, help="Sample label from config/samples.tsv")
    parser.add_argument("--reference", default=str(REPO_ROOT / "data" / "reference" / "GRCh38.primary_assembly.fa.bgz"))
    parser.add_argument("--force", action="store_true", help="Re-run the BAM genotype check even if the alignment is unchanged")
    args = parser.parse_args()
    label = args.sample

    samples = {r["label"]: r for r in read_samples()}
    if label not in samples:
        sys.exit(f"[{label}] Not in {SAMPLES_TSV}")
    declared = (samples[label].get("sex") or "").strip().upper()
    vcf = variants_vcf(label)
    if not vcf.exists():
        sys.exit(f"[{label}] Step 3 needs {vcf} (run the indexing/normalisation steps first)")

    rep_dir = REPO_ROOT / "results" / "reports" / label
    qc_dir = REPO_ROOT / "results" / "qc" / label
    work = REPO_ROOT / "data" / "processed" / label / "identity_qc"
    rep_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    previous = load_identity_qc(rep_dir) or {}
    bcftools = shutil.which("bcftools") or "bcftools"

    print(f"[{label}] Step 3: genetic sex from the VCF's chrX calls...")
    vm = vcf_chrx_metrics(vcf)
    vcf_sex = sex_from_vcf(vm)
    sex_check: dict[str, Any] = {"declared": declared, "vcf": {**vm, "inferred": vcf_sex, "status": sex_status(declared, vcf_sex)}}

    bam = find_alignment(label)
    if bam is None or not bam.exists():
        reason = "no BAM/CRAM in data/alignments" if bam is None else f"{bam.name} points to a missing file (drive not mounted?)"
        prev_bam = previous.get("bam_check") or {}
        if bam is not None and prev_bam.get("status") not in (None, "NOT_CHECKED"):
            bam_check = {**prev_bam, "carried_over": True, "reason": reason}
            if (previous.get("sex_check") or {}).get("bam"):
                sex_check["bam"] = previous["sex_check"]["bam"]
        else:
            bam_check = {"status": "NOT_CHECKED", "reason": reason}
        print(f"[{label}] BAM not checked: {reason}.")
    else:
        fp = alignment_fingerprint(bam)
        vcf_mtime = int(vcf.stat().st_mtime)
        with pysam.AlignmentFile(str(bam)) as af:
            has_chr = "chr1" in set(af.references[:30])
            stats = [(s.contig, af.get_reference_length(s.contig), s.mapped) for s in af.get_index_statistics()]
            rg_samples = sorted({rg.get("SM", "") for rg in af.header.to_dict().get("RG", []) if rg.get("SM")})
        with open(qc_dir / "bam.idxstats.tsv", "w", encoding="utf-8", newline="\n") as fh:
            for contig, length, mapped in stats:
                fh.write(f"{contig}\t{length}\t{mapped}\t0\n")
        ratios = coverage_ratios(stats)
        bam_sex = sex_from_coverage(ratios["x_ratio"], ratios["y_ratio"])
        sex_check["bam"] = {**ratios, "inferred": bam_sex, "status": sex_status(declared, bam_sex)}

        prev_bam = previous.get("bam_check") or {}
        same = all(prev_bam.get(k) == v for k, v in fp.items()) and prev_bam.get("vcf_mtime") == vcf_mtime
        if same and not args.force and prev_bam.get("status") in ("MATCH", "MISMATCH", "INCONCLUSIVE"):
            bam_check = {**prev_bam, "carried_over": False}
            print(f"[{label}] Alignment and VCF unchanged since {previous.get('run_date')}: reusing BAM check ({prev_bam['status']}).")
        else:
            print(f"[{label}] Re-genotyping the VCF's SNVs from {bam.name} (bcftools mpileup/call)...")
            res = check_bam_against(bam, vcf, Path(args.reference), work, has_chr, bcftools)
            bam_check = {**fp, "vcf_mtime": vcf_mtime, **res, "carried_over": False}
            if res["status"] != "MATCH":
                others = []
                for other in samples:
                    if other == label or not variants_vcf(other).exists():
                        continue
                    print(f"[{label}] Trying the BAM against {other}'s VCF...")
                    o = check_bam_against(bam, variants_vcf(other), Path(args.reference), work / other, has_chr, bcftools)
                    others.append({"label": other, **o})
                bam_check["other_samples"] = others
        sid = samples[label].get("sample_id", "")
        bam_check["read_group_samples"] = rg_samples
        bam_check["read_group_matches"] = bool(sid) and sid in rg_samples

    statuses = [sex_check[k]["status"] for k in ("bam", "vcf") if k in sex_check]
    status = overall_status(bam_check["status"], statuses)

    parts = []
    if bam_check["status"] == "NOT_CHECKED":
        parts.append("BAM not checked")
    else:
        parts.append(f"BAM {'matches' if bam_check['status'] == 'MATCH' else 'DOES NOT match' if bam_check['status'] == 'MISMATCH' else 'is inconclusive against'} the VCF")
        swapped = [o["label"] for o in bam_check.get("other_samples") or [] if o["status"] == "MATCH"]
        if swapped:
            parts.append(f"the BAM matches {', '.join(swapped)}'s VCF instead (swapped link?)")
    inferred = {sex_check[k]["inferred"] for k in ("bam", "vcf") if k in sex_check}
    if "MISMATCH" in statuses:
        parts.append(f"genetic sex {'/'.join(sorted(inferred))} does NOT match samples.tsv ({declared or '?'})")
    elif "UNCLEAR" in statuses:
        parts.append(f"genetic sex unclear ({'/'.join(sorted(inferred))}; samples.tsv says {declared or '?'})")
    else:
        parts.append(f"genetic sex {declared} matches samples.tsv")

    qc = {
        "sample": label,
        "sample_id": samples[label].get("sample_id", ""),
        "run_date": datetime.date.today().isoformat(),
        "status": status,
        "summary": "; ".join(parts) + ".",
        "bam_check": bam_check,
        "sex_check": sex_check,
    }
    (rep_dir / QC_JSON).write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    (rep_dir / QC_MD).write_text(render_identity_report(label, qc), encoding="utf-8")
    print(f"[{label}] Sample identity QC: {status} — {qc['summary']} ({rep_dir / QC_MD})")
    if status == "FAIL":
        bar = "!" * 78
        print(f"\n{bar}\n[{label}] SAMPLE IDENTITY QC FAILED: {qc['summary']}\n"
              f"[{label}] BAM-based results (Step 10, and their use in Steps 7, 8 and 11) are blocked for this sample.\n"
              f"[{label}] Check the data/alignments/{label}/ link and the sex column in config/samples.tsv.\n{bar}\n",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
