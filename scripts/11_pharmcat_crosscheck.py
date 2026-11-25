#!/usr/bin/env python3
"""
Step 11: PharmCAT cross-check of the medication card.
  1. Builds a PharmCAT input VCF: the sample's variants near PharmCAT's defining positions, plus explicit
     0/0 calls where the gVCF shows a well-covered reference block (absent positions stay no-calls).
  2. Passes Step 10's Cyrius CYP2D6 and untied T1K HLA-A / HLA-B calls as PharmCAT outside calls.
  3. Runs pharmcat_pipeline (installed by scripts/setup_pharmcat.sh).
  4. Compares PharmCAT gene by gene with pgx_medication_card.tsv and writes
     results/reports/<label>/pharmcat_crosscheck.{md,tsv}, with PharmCAT's own reports in pharmcat/.
Skips cleanly when PharmCAT is not installed or Step 7 has not written the card yet.
"""

from __future__ import annotations

import argparse
import datetime
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.genomics import read_sample_sex
from lib.identity import bam_results_usable, load_bam_results
from lib.pharmcat import (
    REF_BLOCK_MIN_DP,
    build_outside_calls,
    compare_card_with_pharmcat,
    describe_input_at,
    parse_calls_only_tsv,
    parse_phenotype_json,
    render_crosscheck_report,
    write_outside_calls,
    write_pharmcat_input,
)

REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_TOOLS_DIR = Path.home() / ".local" / "opt" / "wgs_tools"


def find_pharmcat() -> tuple[str | None, Path | None]:
    """Return (pharmcat_pipeline path, positions VCF) or Nones when setup_pharmcat.sh has not run."""
    exe = shutil.which("pharmcat_pipeline") or str(Path.home() / ".local" / "bin" / "pharmcat_pipeline")
    if not Path(exe).exists():
        return None, None
    pc_dir = Path(os.environ.get("WGS_TOOLS_DIR", str(DEFAULT_TOOLS_DIR))) / "pharmcat"
    positions = sorted(pc_dir.glob("pharmcat_positions_*.vcf.bgz"))
    positions = [p for p in positions if "uniallelic" not in p.name]
    return exe, (positions[-1] if positions else None)


def main() -> None:
    parser = argparse.ArgumentParser(description="Cross-check the medication card against PharmCAT")
    parser.add_argument("--sample", required=True, help="Sample label from config/samples.tsv")
    parser.add_argument("--min-ref-dp", type=int, default=REF_BLOCK_MIN_DP, help="Minimum gVCF reference-block depth for a 0/0 call")
    args = parser.parse_args()
    label = args.sample

    proc_dir = REPO_ROOT / "data" / "processed" / label
    rep_dir = REPO_ROOT / "results" / "reports" / label
    card_tsv = rep_dir / "pgx_medication_card.tsv"
    gvcf = proc_dir / "snp_indel.genome.vcf.gz"
    variants = proc_dir / "snp_indel.variants.norm.vcf.gz"

    exe, positions = find_pharmcat()
    if not exe or not positions:
        print(f"[{label}] PharmCAT not installed (run 'bash scripts/setup_pharmcat.sh' to enable Step 11).")
        return
    missing = [p.name for p in (card_tsv, gvcf, variants) if not p.exists()]
    if missing:
        print(f"[{label}] Step 11 skipped: missing {', '.join(missing)} (run Steps 1-7 first).")
        return

    work = proc_dir / "pharmcat"
    if work.exists():
        shutil.rmtree(work)
    in_dir = work / "in"
    in_dir.mkdir(parents=True)
    input_vcf = in_dir / f"{label}.vcf.bgz"
    print(f"[{label}] Building PharmCAT input from the variants VCF and gVCF reference blocks...")
    male = read_sample_sex(REPO_ROOT / "config" / "samples.tsv", label) == "M"
    stats = write_pharmcat_input(variants, gvcf, positions, input_vcf, min_dp=args.min_ref_dp, male=male)

    bam = load_bam_results(rep_dir)
    if not bam_results_usable(rep_dir):
        print(f"[{label}] Sample identity QC failed: running PharmCAT without the BAM outside calls.")
    outside = build_outside_calls(bam)
    if outside:
        # pharmcat_pipeline copies *.outside.tsv only when its input is a directory, and pairs
        # <name>.outside.tsv with the VCF of the same base name
        write_outside_calls(in_dir / f"{label}.outside.tsv", outside)

    print(f"[{label}] Running PharmCAT (outside calls: {', '.join(sorted(outside)) or 'none'})...")
    out_dir = work / "out"
    cmd = [exe, str(in_dir), "-o", str(out_dir), "-bf", label, "-reporterJson", "-reporterHtml", "-reporterCallsOnlyTsv"]
    # PharmCAT's bundled logback.xml writes pharmcat.log (and daily rollovers) to a relative path, i.e. into the
    # working directory, so run it from the work dir: the log stays in data/processed/<label>/pharmcat/ instead
    # of the repo root. All paths in cmd are absolute, and pharmcat.jar is found next to the pipeline script.
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=work)
    tsvs = sorted(p for p in out_dir.glob("*.tsv") if not p.name.endswith(".outside.tsv")) if out_dir.exists() else []
    calls_tsv = tsvs[0] if tsvs else None
    if proc.returncode != 0 or calls_tsv is None:
        lines = [ln for ln in (proc.stdout + "\n" + proc.stderr).splitlines() if "bcf_hdr_check_sanity" not in ln]
        tail = "\n".join(lines)[-2500:].strip() or f"exit {proc.returncode}"
        print(f"[{label}] PharmCAT failed:\n{tail}")
        sys.exit(1)

    version = ""
    for line in (proc.stdout + proc.stderr).splitlines():
        if "PharmCAT" in line and any(ch.isdigit() for ch in line):
            version = line.strip().removeprefix("PharmCAT version:").strip()
            break

    card_rows = pd.read_csv(card_tsv, sep="\t", dtype=str).fillna("").to_dict("records")
    tsv_rows = parse_calls_only_tsv(calls_tsv)
    pheno_json = out_dir / f"{label}.phenotype.json"
    pc_rows = parse_phenotype_json(pheno_json, tsv_rows) if pheno_json.exists() else tsv_rows
    rows = compare_card_with_pharmcat(
        card_rows, pc_rows, outside,
        ugt1a1_ta=(bam or {}).get("ugt1a1_ta_repeat"),
        ugt1a1_pharmcat_input=describe_input_at(input_vcf),
    )
    meta = {
        "run_date": datetime.date.today().isoformat(),
        "pharmcat_version": version or positions.name.removeprefix("pharmcat_positions_").removesuffix(".vcf.bgz"),
        "min_dp": args.min_ref_dp,
        **stats,
    }
    pd.DataFrame(rows).to_csv(rep_dir / "pharmcat_crosscheck.tsv", sep="\t", index=False)
    (rep_dir / "pharmcat_crosscheck.md").write_text(render_crosscheck_report(label, rows, meta), encoding="utf-8")

    keep = rep_dir / "pharmcat"
    if keep.exists():
        shutil.rmtree(keep)
    keep.mkdir()
    for f in out_dir.iterdir():
        if f.suffix in (".html", ".json", ".tsv") and not f.name.endswith(".outside.tsv"):
            shutil.copy2(f, keep / f.name)
    disagree = sum(1 for r in rows if r["status"] == "DISAGREE")
    print(f"[{label}] PharmCAT cross-check written to {rep_dir / 'pharmcat_crosscheck.md'} ({disagree} disagreement(s)).")


if __name__ == "__main__":
    main()
