#!/usr/bin/env python3
"""
Step 13: annual snapshot and "what changed since the last run" report.
  1. Looks up every carried variant in the full normalized ClinVar and keeps those classed P, LP, P/LP,
     conflicting or VUS (with VariationID, review status, stars, condition, genotype and call QC).
  2. Adds Step 5's LoF candidates, the Step 7 medication card, Step 10 HLA associations and Step 11 PharmCAT calls,
     and writes data/snapshots/<label>/<run_date>.json (gitignored; one per run date, --force to overwrite).
  3. Diffs it against the latest earlier snapshot and writes results/reports/<label>/changes_since_last_run.{md,tsv}.
Run by run_annual_checkup.sh after every sample. To compare two snapshot files directly:
  13_annual_snapshot.py --sample <label> --old A.json --new B.json [--out-dir DIR]
"""

from __future__ import annotations

import argparse
import datetime
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.genomics import get_reference_metadata, parse_info_str
from lib.identity import load_identity_qc
from lib.snapshot import (
    build_snapshot,
    clinvar_class,
    diff_snapshots,
    find_previous_snapshot,
    is_tracked,
    load_snapshot,
    render_changes_report,
    variant_entry,
    write_snapshot,
)

REPO_ROOT = SCRIPT_DIR.parent
REPORT_STEM = "changes_since_last_run"


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size < 2:
        return []
    return pd.read_csv(path, sep="\t", dtype=str).fillna("").to_dict("records")


def clinvar_at(cv: pysam.VariantFile, chrom: str, pos: int, ref: str, alt: str) -> pysam.VariantRecord | None:
    try:
        for rec in cv.fetch(chrom, pos - 1, pos):
            if rec.pos == pos and rec.ref == ref and rec.alts and alt in rec.alts:
                return rec
    except ValueError:  # contig absent from ClinVar
        return None
    return None


def extract_tracked_variants(var_vcf: Path, clinvar_vcf: Path, tmp: Path) -> dict[str, dict]:
    """Carried variants that ClinVar classes P, LP, P/LP, conflicting or VUS, keyed chrom-pos-ref-alt."""
    # Filtered in Python rather than with a bcftools '~' expression: some bcftools builds (1.22) print a debug
    # line to stderr for every record a regex filter tests, millions of lines over the full ClinVar.
    targets = tmp / "clinvar_tracked_positions.tsv"
    with open(targets, "w", encoding="utf-8", newline="\n") as fh:
        q = subprocess.Popen(["bcftools", "query", "-f", "%CHROM\t%POS\t%INFO/CLNSIG\n", str(clinvar_vcf)],
                             stdout=subprocess.PIPE, text=True)
        last = None
        for line in q.stdout:
            chrom, pos, clnsig = line.rstrip("\n").split("\t", 2)
            if (chrom, pos) != last and is_tracked(clinvar_class(clnsig)):
                fh.write(f"{chrom}\t{pos}\n")
                last = (chrom, pos)
        if q.wait() != 0:
            raise RuntimeError("bcftools query on ClinVar failed")
    view = subprocess.run(
        ["bcftools", "view", "-H", "-T", str(targets), "-i", 'GT="alt"', str(var_vcf)],
        capture_output=True, text=True,
    )
    if view.returncode != 0:
        raise RuntimeError(f"bcftools view on the sample VCF failed: {view.stderr.strip()[-500:]}")
    out: dict[str, dict] = {}
    with pysam.VariantFile(str(clinvar_vcf)) as cv:
        for line in view.stdout.splitlines():
            f = line.split("\t")
            chrom, pos, ref, alts, filt = f[0], int(f[1]), f[3], f[4], f[6]
            fmt = f[8].split(":")
            sample = dict(zip(fmt, f[9].split(":"))) if len(f) > 9 else {}
            for alt in alts.split(","):
                rec = clinvar_at(cv, chrom, pos, ref, alt)
                if rec is None:
                    continue
                clnsig = parse_info_str(rec.info.get("CLNSIG"))
                if not is_tracked(clinvar_class(clnsig)):
                    continue
                out[f"{chrom}-{pos}-{ref}-{alt}"] = variant_entry(
                    clnsig=clnsig,
                    revstat=rec.info.get("CLNREVSTAT"),
                    clinvar_id=str(rec.id or ""),
                    geneinfo=parse_info_str(rec.info.get("GENEINFO")),
                    condition=parse_info_str(rec.info.get("CLNDN")),
                    gt=sample.get("GT", ""),
                    filt=filt,
                    dp=sample.get("DP", ""),
                    gq=sample.get("GQ", ""),
                )
    return out


def make_lookup(var_vcf: Path, clinvar_vcf: Path):
    """lookup(key) for diff_snapshots: current ClinVar class and whether this run's VCF still calls the variant."""
    cv = pysam.VariantFile(str(clinvar_vcf))
    vf = pysam.VariantFile(str(var_vcf))

    def lookup(key: str) -> dict:
        chrom, pos, ref, alt = key.split("-", 3)
        rec = clinvar_at(cv, chrom, int(pos), ref, alt)
        carried = False
        try:
            for r in vf.fetch(chrom, int(pos) - 1, int(pos)):
                if r.pos == int(pos) and r.ref == ref and r.alts and alt in r.alts:
                    gt = next(iter(r.samples.values())).get("GT") or ()
                    idx = r.alts.index(alt) + 1
                    carried = idx in gt
        except ValueError:
            pass
        return {"class": clinvar_class(parse_info_str(rec.info.get("CLNSIG"))) if rec else "", "carried": carried}

    return lookup


def git_commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT, capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception:
        return ""


def pharmcat_version(rep_dir: Path) -> str:
    md = rep_dir / "pharmcat_crosscheck.md"
    if md.exists():
        m = re.search(r"\*\*PharmCAT:\*\*\s*([^\s|]+)", md.read_text(encoding="utf-8"))
        if m:
            return m.group(1)
    return ""


def write_report(label: str, out_dir: Path, old: dict | None, new: dict, rows: list[dict]) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    md_path = out_dir / f"{REPORT_STEM}.md"
    with open(md_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(render_changes_report(label, old, new, rows) + "\n")
    cols = ["section", "key", "gene", "old", "new", "direction", "note"]
    pd.DataFrame(rows, columns=cols).to_csv(out_dir / f"{REPORT_STEM}.tsv", sep="\t", index=False, lineterminator="\n")
    return md_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Annual snapshot and 'what changed since last run' report")
    parser.add_argument("--sample", required=True, help="Sample label from config/samples.tsv")
    parser.add_argument("--force", action="store_true", help="Overwrite today's snapshot if it exists")
    parser.add_argument("--old", type=Path, help="Diff only: earlier snapshot JSON")
    parser.add_argument("--new", type=Path, help="Diff only: later snapshot JSON")
    parser.add_argument("--out-dir", type=Path, help="Diff only: where to write the report (default: the sample's reports)")
    args = parser.parse_args()
    label = args.sample
    rep_dir = REPO_ROOT / "results" / "reports" / label

    if args.old or args.new:
        if not (args.old and args.new):
            parser.error("--old and --new go together")
        old, new = load_snapshot(args.old), load_snapshot(args.new)
        rows = diff_snapshots(old, new)
        path = write_report(label, args.out_dir or rep_dir, old, new, rows)
        print(f"[{label}] {len(rows)} change(s): {path}")
        return

    proc_dir = REPO_ROOT / "data" / "processed" / label
    ref_dir = REPO_ROOT / "data" / "reference"
    var_vcf = proc_dir / "snp_indel.variants.norm.vcf.gz"
    clinvar_vcf = ref_dir / "clinvar_GRCh38.norm.vcf.gz"
    missing = [p.name for p in (var_vcf, clinvar_vcf) if not p.exists()]
    if missing:
        print(f"[{label}] Step 13 skipped: missing {', '.join(missing)} (run setup_references.sh and Steps 1-7 first).")
        return

    run_date = datetime.date.today().isoformat()
    snap_dir = REPO_ROOT / "data" / "snapshots" / label
    snap_path = snap_dir / f"{run_date}.json"
    if snap_path.exists() and not args.force:
        print(f"[{label}] Snapshot {snap_path.name} already exists (use --force to overwrite); re-writing the report from it.")
        new = load_snapshot(snap_path)
    else:
        print(f"[{label}] Looking up carried variants in the full ClinVar (P/LP, conflicting, VUS)...")
        with tempfile.TemporaryDirectory(dir=proc_dir) as tmp:
            variants = extract_tracked_variants(var_vcf, clinvar_vcf, Path(tmp))
        qc_rows = read_tsv(rep_dir / "carrier_screening_plp_enriched.tsv") + read_tsv(rep_dir / "snp_indel_conflicting_plp.tsv")
        for r in qc_rows:
            v = variants.get(r.get("variant_id", ""))
            if v is not None:
                v["qc_tier"] = r.get("qc_tier", "") or v["qc_tier"]
                v["artifact_flag"] = r.get("artifact_flag", "") or v["artifact_flag"]
        enriched = rep_dir / "carrier_screen_enriched.json"
        lof = json.loads(enriched.read_text(encoding="utf-8")).get("novel_acmg_lof", []) if enriched.exists() else []
        refs = get_reference_metadata(ref_dir)
        identity = load_identity_qc(rep_dir) or {}
        meta = {
            "clinvar_date": refs.get("clinvar_date", ""),
            "gencc_date": refs.get("gencc_date", ""),
            "acmg_sf_version": refs.get("acmg_sf_version", ""),
            "pharmcat_version": pharmcat_version(rep_dir),
            "git_commit": git_commit(),
            "identity_qc": identity.get("status", ""),
        }
        new = build_snapshot(
            label, run_date, meta, variants, lof,
            card_rows=read_tsv(rep_dir / "pgx_medication_card.tsv"),
            hla_rows=read_tsv(rep_dir / "hla_disease_associations.tsv"),
            pharmcat_rows=read_tsv(rep_dir / "pharmcat_crosscheck.tsv"),
        )
        write_snapshot(snap_path, new)
        print(f"[{label}] Snapshot written: {snap_path} ({len(variants)} tracked variants)")

    prev_path = find_previous_snapshot(snap_dir, run_date)
    old = load_snapshot(prev_path) if prev_path else None
    rows = diff_snapshots(old, new, make_lookup(var_vcf, clinvar_vcf)) if old else []
    path = write_report(label, rep_dir, old, new, rows)
    against = f"against {prev_path.name}" if prev_path else "first snapshot, nothing to compare"
    print(f"[{label}] What-changed report ({against}): {path} ({len(rows)} change(s))")


if __name__ == "__main__":
    main()
