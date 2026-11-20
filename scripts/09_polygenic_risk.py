#!/usr/bin/env python3
"""
Step 9: Gold-Standard Local Polygenic Risk Scores (PRS) with 1000 Genomes GRCh38 Calibration
  - --build-ref:
      1. Parses curated GRCh38 harmonized PGS Catalog scoring files.
      2. Intersects against 1000 Genomes Phase 3 GRCh38 biallelic ACGT SNPs.
      3. Orients weights to GRCh38 ALT (preserving exact raw score via constant offset C)
         so PLINK 2.0 scores all traits in a single pass across unrelated 1000G samples.
      4. Computes empirical score distributions per superpopulation (EUR/AFR/EAS/SAS/AMR/ALL)
         and superpopulation allele frequencies for ~1,500 ancestry-informative markers.
  - --sample <label>:
      1. Estimates closest 1000G superpopulation via HWE genotype log-likelihood.
      2. Computes exact 30x WGS raw PRS on the identical SNP set as 1000G (handling 0/0 REF
         sites exactly without missing-genotype imputation bias).
      3. Computes empirical percentile, z-score, and relative risk tier against the matched
         1000G superpopulation (and saves all 5 superpopulations in TSV/JSON).
      4. Writes results/reports/<sample>/polygenic_risk_report.md, .tsv, and .json.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import math
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.genomics import get_reference_metadata

REPO_ROOT = SCRIPT_DIR.parent
REFERENCE_DIR = REPO_ROOT / "data" / "reference"
PRS_DIR = REFERENCE_DIR / "prs"
KG_DIR = PRS_DIR / "1000g_hg38"
SCORES_DIR = PRS_DIR / "scores"
PREPARED_DIR = PRS_DIR / "prepared"

AUTOSOMES = {str(i) for i in range(1, 23)}
ACGT = {"A", "C", "G", "T"}
SUPERPOPS = ("EUR", "AFR", "EAS", "SAS", "AMR")

# Curated GRCh38 harmonized PGS Catalog clinical panel
# Tuple: (pgs_id, trait_label, category, sex_scope, citation_note)
PGS_PANEL: list[tuple[str, str, str, str, str]] = [
    ("PGS000018", "Coronary Artery Disease", "Cardiometabolic", "Both", "Inouye et al. metaGRS (~1.7M SNPs)"),
    ("PGS000014", "Type 2 Diabetes", "Cardiometabolic", "Both", "Khera et al. genome-wide (~6.9M SNPs)"),
    ("PGS000016", "Atrial Fibrillation", "Cardiometabolic", "Both", "Khera et al. genome-wide (~6.7M SNPs)"),
    ("PGS000688", "LDL Cholesterol", "Cardiometabolic", "Both", "Sinnott-Armstrong et al. (~16k SNPs)"),
    ("PGS002257", "Systolic Blood Pressure", "Cardiometabolic", "Both", "PGS002257 (~884 SNPs)"),
    ("PGS000004", "Breast Cancer", "Oncology", "Female", "Mavaddat et al. 313-SNP clinical benchmark"),
    ("PGS000662", "Prostate Cancer", "Oncology", "Male", "Conti et al. multi-ancestry 269-SNP benchmark"),
    ("PGS000734", "Colorectal Cancer", "Oncology", "Both", "Huyghe / Thomas et al. (95 SNPs)"),
    ("PGS000334", "Late-Onset Alzheimer's Disease", "Neurological", "Both", "22-SNP PRS (interpret alongside APOE diplotype)"),
]


def _find_plink2() -> str:
    local = REPO_ROOT / ".venv-wsl" / "bin" / "plink2"
    if local.exists():
        return str(local)
    found = shutil.which("plink2")
    if found:
        return found
    raise FileNotFoundError("plink2 not found. Run `bash scripts/setup_prs_references.sh` first.")


def _stream_pgs_to_chrom_buckets(
    score_gz: Path,
    pgs_idx: int,
    chrom_files: dict[str, Any],
) -> int:
    """
    Stream a PGS Catalog GRCh38 harmonized scoring file directly into per-chromosome
    bucket files so peak RAM stays tiny even for 7M-SNP genome-wide scores.
    """
    count = 0
    with gzip.open(score_gz, "rt", encoding="utf-8") as f:
        header: list[str] = []
        for line in f:
            if line.startswith("#"):
                continue
            header = line.rstrip("\r\n").split("\t")
            break
        if not header:
            return 0

        col = {name: i for i, name in enumerate(header)}
        chr_idx = col.get("hm_chr", col.get("chr_name"))
        pos_idx = col.get("hm_pos", col.get("chr_position"))
        eff_idx = col.get("effect_allele")
        oth_idx = col.get("hm_inferOtherAllele", col.get("other_allele"))
        alt_oth_idx = col.get("other_allele")
        wt_idx = col.get("effect_weight")
        if None in (chr_idx, pos_idx, eff_idx, wt_idx):
            raise ValueError(f"Missing required PGS columns in {score_gz.name}: {header}")

        min_len = max(chr_idx, pos_idx, eff_idx, wt_idx)
        for line in f:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) <= min_len:
                continue
            chrom = parts[chr_idx].replace("chr", "")
            if chrom not in AUTOSOMES:
                continue
            pos_s = parts[pos_idx]
            if not pos_s or not pos_s.isdigit():
                continue
            eff = parts[eff_idx].upper()
            if eff not in ACGT:
                continue
            oth = ""
            if alt_oth_idx is not None and alt_oth_idx < len(parts):
                oth = parts[alt_oth_idx].upper()
            if oth not in ACGT and oth_idx is not None and oth_idx < len(parts):
                oth = parts[oth_idx].upper()
            if oth and oth not in ACGT:
                continue
            try:
                wt = float(parts[wt_idx])
            except ValueError:
                continue
            if wt == 0.0 or math.isnan(wt) or math.isinf(wt):
                continue
            chrom_files[chrom].write(f"{pos_s}\t{pgs_idx}\t{eff}\t{oth}\t{wt:.8g}\n")
            count += 1
    return count


def orient_weight_to_alt(
    effect_allele: str,
    other_allele: str,
    weight: float,
    ref: str,
    alt: str,
) -> tuple[float, float] | None:
    """
    Orient a GRCh38 harmonized SNP weight to the ALT allele.
    Returns (alt_weight, ref_offset) such that:
        dosage(effect_allele) * weight == ref_offset + dosage(ALT) * alt_weight
    Returns None if alleles are incompatible with (ref, alt).
    """
    if effect_allele == alt and (not other_allele or other_allele == ref):
        return (weight, 0.0)
    if effect_allele == ref and (not other_allele or other_allele == alt):
        return (-weight, 2.0 * weight)
    return None


def _load_psam_superpops(psam_path: Path) -> dict[str, str]:
    """Read 1000G hg38_corrected.psam -> {IID: SuperPop}."""
    pop_map: dict[str, str] = {}
    with open(psam_path, encoding="utf-8") as f:
        header = f.readline().lstrip("#").rstrip("\r\n").split("\t")
        col = {h.upper(): i for i, h in enumerate(header)}
        iid_idx = col.get("IID", 0)
        sp_idx = col.get("SUPERPOP", col.get("SUPER_POP"))
        if sp_idx is None:
            raise ValueError(f"SuperPop column not found in {psam_path}: {header}")
        for line in f:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) > max(iid_idx, sp_idx):
                pop_map[parts[iid_idx]] = parts[sp_idx]
    return pop_map


def build_reference_distributions() -> None:
    """Build matched SNP score tables, 1000G empirical distributions, and ancestry AFs."""
    plink2 = _find_plink2()
    pgen = KG_DIR / "all_hg38.pgen"
    pvar_zst = KG_DIR / "all_hg38_noannot.pvar.zst"
    psam = KG_DIR / "hg38_corrected.psam"
    king_out = KG_DIR / "deg2_hg38.king.cutoff.out.id"

    for req in (pgen, pvar_zst, psam):
        if not req.exists():
            raise FileNotFoundError(f"Missing 1000G reference file: {req}. Run setup_prs_references.sh first.")

    PREPARED_DIR.mkdir(parents=True, exist_ok=True)
    pop_map = _load_psam_superpops(psam)
    print(f"[PRS Build] Loaded {len(pop_map):,} 1000 Genomes samples from {psam.name}.")

    combined_score_tsv = PREPARED_DIR / "combined_pgs_alt_scores.tsv"
    extract_ids_txt = PREPARED_DIR / "combined_pgs_var_ids.txt"
    ancestry_ids_txt = PREPARED_DIR / "ancestry_aim_ids.txt"
    ancestry_bed = PREPARED_DIR / "ancestry_aim.range"
    build_meta_path = PREPARED_DIR / "build_meta.json"
    pgs_ids = [p[0] for p in PGS_PANEL]

    if build_meta_path.exists() and combined_score_tsv.exists() and extract_ids_txt.exists():
        print("[PRS Build] Reusing prepared score table & metadata from disk...")
        meta = json.loads(build_meta_path.read_text(encoding="utf-8"))
        pgs_totals = meta["pgs_totals"]
        matched_counts = meta["matched_counts"]
        ref_offsets = meta["ref_offsets"]
        aim_records = [tuple(x) for x in meta["aim_records"]]
    else:
        # Step 1: Stream all PGS files into 22 per-chromosome bucket files (<20 MB peak RAM)
        pgs_totals = {}
        tmp_chrom_dir = PREPARED_DIR / "_tmp_chrom_buckets"
        tmp_chrom_dir.mkdir(parents=True, exist_ok=True)
        chrom_paths = {c: tmp_chrom_dir / f"chr{c}.tsv" for c in AUTOSOMES}
        chrom_handles = {c: open(p, "w", encoding="utf-8") for c, p in chrom_paths.items()}
        try:
            for idx, (pid, trait, *_rest) in enumerate(PGS_PANEL):
                sfile = SCORES_DIR / f"{pid}_hmPOS_GRCh38.txt.gz"
                if not sfile.exists():
                    raise FileNotFoundError(f"Missing PGS scoring file: {sfile}")
                n_snps = _stream_pgs_to_chrom_buckets(sfile, idx, chrom_handles)
                pgs_totals[pid] = n_snps
                print(f"  Streamed {pid} ({trait}): {n_snps:,} autosomal ACGT SNPs")
        finally:
            for fh in chrom_handles.values():
                fh.close()

        def _load_chrom_bucket(chrom: str) -> dict[int, list[tuple[int, str, str, float]]]:
            bucket: dict[int, list[tuple[int, str, str, float]]] = defaultdict(list)
            seen_per_score: set[tuple[int, int]] = set()
            with open(chrom_paths[chrom], encoding="utf-8") as fh:
                for line in fh:
                    pos_s, pidx_s, eff, oth, wt_s = line.rstrip("\n").split("\t")
                    pos = int(pos_s)
                    pidx = int(pidx_s)
                    if (pos, pidx) in seen_per_score:
                        continue
                    seen_per_score.add((pos, pidx))
                    bucket[pos].append((pidx, eff, oth, float(wt_s)))
            return bucket

        # Step 2: Stream 1000G pvar.zst once, loading only ONE chromosome bucket into RAM at a time
        print("[PRS Build] Streaming 1000G GRCh38 pvar.zst (1 chromosome in RAM at a time)...")
        ref_offsets = [0.0] * len(PGS_PANEL)
        matched_counts = [0] * len(PGS_PANEL)
        total_unique_matched = 0
        aim_records = []  # (vid, chrom, pos, ref, alt)
        last_aim_pos: dict[str, int] = defaultdict(lambda: -10_000_000)

        proc = subprocess.Popen(
            [plink2, "--zst-decompress", str(pvar_zst)],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            bufsize=1024 * 1024,
        )
        assert proc.stdout is not None

        current_chrom: str | None = None
        chrom_dict: dict[int, list[tuple[int, str, str, float]]] = {}
        seen_chrom_vids: set[str] = set()

        try:
            with open(combined_score_tsv, "w", encoding="utf-8") as fout, \
                 open(extract_ids_txt, "w", encoding="utf-8") as fext:
                fout.write("ID\tALT\t" + "\t".join(pgs_ids) + "\n")
                for line in proc.stdout:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\r\n").split("\t")
                    if len(parts) < 5:
                        continue
                    chrom = parts[0].replace("chr", "")
                    if chrom not in AUTOSOMES:
                        continue
                    if chrom != current_chrom:
                        current_chrom = chrom
                        chrom_dict = _load_chrom_bucket(chrom)
                        seen_chrom_vids = set()
                        print(f"  Matching chr{chrom} ({len(chrom_dict):,} unique PGS positions)...")

                    pos = int(parts[1])
                    hits = chrom_dict.get(pos)
                    if not hits:
                        continue
                    ref = parts[3].upper()
                    alt = parts[4].upper()
                    if ref not in ACGT or alt not in ACGT:
                        continue

                    vid = f"{chrom}:{pos}:{ref}:{alt}"
                    if vid in seen_chrom_vids:
                        continue

                    row_weights = [0.0] * len(PGS_PANEL)
                    any_match = False
                    for pidx, eff, oth, wt in hits:
                        oriented = orient_weight_to_alt(eff, oth, wt, ref, alt)
                        if oriented is None:
                            continue
                        alt_w, offset = oriented
                        row_weights[pidx] = alt_w
                        ref_offsets[pidx] += offset
                        matched_counts[pidx] += 1
                        any_match = True

                    if any_match:
                        seen_chrom_vids.add(vid)
                        total_unique_matched += 1
                        w_strs = [f"{w:.8g}" if w != 0.0 else "0" for w in row_weights]
                        fout.write(f"{vid}\t{alt}\t" + "\t".join(w_strs) + "\n")
                        fext.write(vid + "\n")
                        if pos - last_aim_pos[chrom] >= 1_800_000:
                            last_aim_pos[chrom] = pos
                            aim_records.append((vid, chrom, pos, ref, alt))

            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError(f"plink2 --zst-decompress failed with exit code {proc.returncode}")
        finally:
            chrom_dict.clear()
            seen_chrom_vids.clear()
            shutil.rmtree(tmp_chrom_dir, ignore_errors=True)

        print(f"[PRS Build] Matched {total_unique_matched:,} unique 1000G GRCh38 SNPs across {len(PGS_PANEL)} PGS scores.")
        for idx, (pid, trait, *_rest) in enumerate(PGS_PANEL):
            pct = 100.0 * matched_counts[idx] / max(1, pgs_totals[pid])
            print(f"  {pid} ({trait}): {matched_counts[idx]:,}/{pgs_totals[pid]:,} ({pct:.1f}%) matched")

        with open(ancestry_ids_txt, "w", encoding="utf-8") as f, \
             open(ancestry_bed, "w", encoding="utf-8") as fbed:
            for vid, chrom, pos, _ref, _alt in aim_records:
                f.write(vid + "\n")
                fbed.write(f"{chrom}\t{pos}\t{pos}\t{vid}\n")

        build_meta_path.write_text(
            json.dumps({
                "pgs_totals": pgs_totals,
                "matched_counts": matched_counts,
                "ref_offsets": ref_offsets,
                "aim_records": aim_records,
            }),
            encoding="utf-8",
        )

    dist_path = PRS_DIR / "ref_distributions.json"
    if dist_path.exists():
        print(f"[PRS Build] Reusing existing 1000G empirical distributions: {dist_path}")
    else:
        # Step 3: Score all unrelated 1000G individuals across all 9 traits in ONE PLINK 2.0 pass
        print("[PRS Build] Scoring 1000 Genomes reference panel with PLINK 2.0 (single pass)...")
        sscore_prefix = PREPARED_DIR / "1000g_ref_scores"
        cmd = [
            plink2,
            "--pgen", str(pgen),
            "--pvar", str(pvar_zst),
            "--psam", str(psam),
            "--autosome",
            "--snps-only", "just-acgt",
            "--max-alleles", "2",
            "--set-all-var-ids", "@:#:$r:$a",
            "--rm-dup", "exclude-all",
            "--extract", str(extract_ids_txt),
            "--score", str(combined_score_tsv), "1", "2", "header-read", "no-mean-imputation", "cols=+scoresums",
            "--score-col-nums", f"3-{2 + len(PGS_PANEL)}",
            "--memory", "8000", "require",
            "--threads", "4",
            "--out", str(sscore_prefix),
        ]
        if king_out.exists():
            cmd.extend(["--remove", str(king_out)])
        subprocess.run(cmd, check=True)

        sscore_df = pd.read_csv(f"{sscore_prefix}.sscore", sep="\t")
        iid_col = "#IID" if "#IID" in sscore_df.columns else "IID"
        sscore_df["SuperPop"] = sscore_df[iid_col].map(pop_map).fillna("UNKNOWN")

        distributions: dict[str, Any] = {
            "n_ref_samples_total": int(len(sscore_df)),
            "superpop_counts": {
                sp: int((sscore_df["SuperPop"] == sp).sum()) for sp in SUPERPOPS
            },
            "scores": {},
        }

        for idx, (pid, trait, category, sex_scope, citation) in enumerate(PGS_PANEL):
            sum_col = f"{pid}_SUM"
            if sum_col not in sscore_df.columns:
                candidates = [c for c in sscore_df.columns if c.startswith(pid) and c.endswith("_SUM")]
                if not candidates:
                    raise KeyError(f"Expected score sum column for {pid} in {list(sscore_df.columns)}")
                sum_col = candidates[0]

            raw_vals_all = sscore_df[sum_col].astype(float).to_numpy() + ref_offsets[idx]
            pop_stats: dict[str, Any] = {}
            for sp in (*SUPERPOPS, "ALL"):
                vals = (
                    raw_vals_all
                    if sp == "ALL"
                    else (sscore_df.loc[sscore_df["SuperPop"] == sp, sum_col].astype(float).to_numpy() + ref_offsets[idx])
                )
                sorted_vals = np.sort(vals)
                pop_stats[sp] = {
                    "n": int(len(sorted_vals)),
                    "mean": float(np.mean(sorted_vals)),
                    "sd": float(np.std(sorted_vals, ddof=1)) if len(sorted_vals) > 1 else 1.0,
                    "sorted_scores": [round(float(x), 6) for x in sorted_vals],
                }

            distributions["scores"][pid] = {
                "pgs_id": pid,
                "trait": trait,
                "category": category,
                "sex_scope": sex_scope,
                "citation": citation,
                "n_total_pgs": pgs_totals[pid],
                "n_matched_1000g": matched_counts[idx],
                "ref_offset": ref_offsets[idx],
                "populations": pop_stats,
            }

        dist_path.write_text(json.dumps(distributions), encoding="utf-8")
        print(f"[PRS Build] Saved 1000G empirical distributions to {dist_path}")
        del sscore_df, distributions

    # Step 4: Compute superpopulation AFs for the ~1,500 ancestry-informative markers
    # First extract the ~1,500 AIM sites into a tiny pfile via `--extract range`
    aim_tsv = PRS_DIR / "ancestry_afs.tsv"
    if aim_tsv.exists():
        print(f"[PRS Build] Reusing existing ancestry AFs: {aim_tsv}")
        return

    print(f"[PRS Build] Extracting {len(aim_records)} ancestry markers & computing 1000G superpop AFs...")
    aim_pfile = PREPARED_DIR / "aim_subset"
    cmd_sub = [
        plink2,
        "--pgen", str(pgen),
        "--pvar", str(pvar_zst),
        "--psam", str(psam),
        "--autosome",
        "--extract", "range", str(ancestry_bed),
        "--snps-only", "just-acgt",
        "--max-alleles", "2",
        "--set-all-var-ids", "@:#:$r:$a",
        "--rm-dup", "exclude-all",
        "--memory", "5500", "require",
        "--threads", "2",
        "--make-pgen",
        "--out", str(aim_pfile),
    ]
    if king_out.exists():
        cmd_sub.extend(["--remove", str(king_out)])
    subprocess.run(cmd_sub, check=True, stdout=subprocess.DEVNULL)

    aim_af_by_pop: dict[str, dict[str, float]] = {}
    for sp in SUPERPOPS:
        sp_keep = PREPARED_DIR / f"keep_{sp}.txt"
        sp_iids = [iid for iid, p in pop_map.items() if p == sp]
        sp_keep.write_text("\n".join(sp_iids) + "\n", encoding="utf-8")
        af_prefix = PREPARED_DIR / f"aim_af_{sp}"
        cmd_af = [
            plink2,
            "--pfile", str(aim_pfile),
            "--keep", str(sp_keep),
            "--extract", str(ancestry_ids_txt),
            "--freq",
            "--memory", "1000",
            "--threads", "4",
            "--out", str(af_prefix),
        ]
        subprocess.run(cmd_af, check=True, stdout=subprocess.DEVNULL)
        af_df = pd.read_csv(f"{af_prefix}.afreq", sep="\t")
        aim_af_by_pop[sp] = dict(zip(af_df["ID"], af_df["ALT_FREQS"].astype(float)))
        sp_keep.unlink(missing_ok=True)
        for ext in (".afreq", ".log"):
            Path(f"{af_prefix}{ext}").unlink(missing_ok=True)
    for ext in (".pgen", ".pvar", ".psam", ".log"):
        Path(f"{aim_pfile}{ext}").unlink(missing_ok=True)

    aim_rows = []
    for vid, chrom, pos, ref, alt in aim_records:
        afs = {sp: aim_af_by_pop[sp].get(vid) for sp in SUPERPOPS}
        if any(v is None for v in afs.values()):
            continue
        # Only keep informative markers where AF differs across superpopulations
        vals = list(afs.values())
        if max(vals) - min(vals) >= 0.10:
            aim_rows.append({"id": vid, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, **afs})
    aim_tsv = PRS_DIR / "ancestry_afs.tsv"
    pd.DataFrame(aim_rows).to_csv(aim_tsv, sep="\t", index=False)
    print(f"[PRS Build] Saved {len(aim_rows)} ancestry-informative markers to {aim_tsv}")


def estimate_superpopulation(
    personal_vcf: Path,
    ancestry_afs_tsv: Path,
) -> dict[str, Any]:
    """
    Estimate closest 1000 Genomes superpopulation (EUR/AFR/EAS/SAS/AMR) via HWE
    genotype log-likelihood across ~1,000 ancestry-informative SNPs.
    """
    df = pd.read_csv(ancestry_afs_tsv, sep="\t", dtype={"chrom": str})
    if df.empty:
        return {"inferred_superpop": "EUR", "n_markers": 0, "log_likelihoods": {}}

    # Lookup personal ALT dosages for AIMs via targeted tabix fetch
    loglik = {sp: 0.0 for sp in SUPERPOPS}
    eps = 0.005
    n_used = 0
    with pysam.VariantFile(str(personal_vcf)) as vf:
        for _, r in df.iterrows():
            chrom, pos, ref, alt = str(r["chrom"]), int(r["pos"]), str(r["ref"]), str(r["alt"])
            dosage = 0
            try:
                for rec in vf.fetch(chrom, pos - 1, pos):
                    if rec.pos == pos and rec.ref == ref and rec.alts and alt in rec.alts:
                        flt = ";".join(rec.filter.keys()) if len(rec.filter) > 0 else "PASS"
                        if flt != "PASS":
                            continue
                        sample = next(iter(rec.samples.values()))
                        gt = sample.get("GT") or ()
                        alt_idx = rec.alts.index(alt) + 1
                        dosage = sum(1 for a in gt if a == alt_idx)
                        break
            except ValueError:
                continue

            n_used += 1
            for sp in SUPERPOPS:
                p = min(1.0 - eps, max(eps, float(r[sp])))
                q = 1.0 - p
                if dosage == 0:
                    prob = q * q
                elif dosage == 1:
                    prob = 2.0 * p * q
                else:
                    prob = p * p
                loglik[sp] += math.log(prob)

    ranked = sorted(loglik.items(), key=lambda x: x[1], reverse=True)
    best_pop = ranked[0][0]
    return {
        "inferred_superpop": best_pop,
        "n_markers": n_used,
        "log_likelihoods": {k: round(v, 2) for k, v in ranked},
    }


def compute_empirical_percentile(value: float, sorted_ref: list[float]) -> float:
    """Return empirical percentile (0-100) using mid-rank interpolation."""
    if not sorted_ref:
        return 50.0
    lo = bisect.bisect_left(sorted_ref, value)
    hi = bisect.bisect_right(sorted_ref, value)
    mid = 0.5 * (lo + hi)
    return round(100.0 * mid / len(sorted_ref), 1)


def risk_tier_from_percentile(pct: float) -> str:
    if pct >= 95.0:
        return "High (>=95th)"
    if pct >= 80.0:
        return "Moderately Elevated (80th-94th)"
    if pct <= 20.0:
        return "Low (<=20th)"
    return "Average (21st-79th)"


def score_personal_wgs(
    personal_vcf: Path,
    combined_score_tsv: Path,
    distributions: dict[str, Any],
) -> tuple[dict[str, float], dict[str, int]]:
    """
    Compute exact 30x WGS raw PRS across all traits in combined_score_tsv.
    Uses a memory-light chromosome-by-chromosome streaming pass over combined_score_tsv
    paired with targeted tabix reads of personal_vcf so peak RAM stays <500 MB.
    """
    with open(combined_score_tsv, encoding="utf-8") as f:
        header = f.readline().rstrip("\r\n").split("\t")
    pgs_ids = header[2:]

    alt_sums = {pid: 0.0 for pid in pgs_ids}
    nonref_hits = {pid: 0 for pid in pgs_ids}

    with pysam.VariantFile(str(personal_vcf)) as pvf:
        current_chrom = None
        personal_dosages: dict[tuple[int, str, str], int] = {}

        def _load_chrom(chrom: str) -> dict[tuple[int, str, str], int]:
            dos: dict[tuple[int, str, str], int] = {}
            try:
                for rec in pvf.fetch(chrom):
                    if not rec.alts or len(rec.ref) != 1:
                        continue
                    if len(rec.filter) > 0 and "PASS" not in rec.filter:
                        continue
                    sample = next(iter(rec.samples.values()))
                    gt = sample.get("GT")
                    if not gt:
                        continue
                    for alt_i, alt in enumerate(rec.alts, 1):
                        if len(alt) != 1:
                            continue
                        cnt = sum(1 for a in gt if a == alt_i)
                        if cnt > 0:
                            dos[(rec.pos, rec.ref, alt)] = cnt
            except ValueError:
                pass
            return dos

        with open(combined_score_tsv, encoding="utf-8") as f:
            next(f)  # skip header
            for line in f:
                parts = line.rstrip("\r\n").split("\t")
                vid = parts[0]
                c, pos_s, ref, alt = vid.split(":")
                if c != current_chrom:
                    current_chrom = c
                    personal_dosages = _load_chrom(c)
                dosage = personal_dosages.get((int(pos_s), ref, alt))
                if not dosage:
                    continue
                for i, pid in enumerate(pgs_ids):
                    w_str = parts[2 + i]
                    if w_str != "0":
                        alt_sums[pid] += dosage * float(w_str)
                        nonref_hits[pid] += 1

    raw_scores: dict[str, float] = {}
    for pid in pgs_ids:
        offset = float(distributions["scores"][pid]["ref_offset"])
        raw_scores[pid] = alt_sums[pid] + offset
    return raw_scores, nonref_hits


def render_prs_report(
    sample_label: str,
    ref_meta: dict[str, str],
    ancestry: dict[str, Any],
    rows: list[dict[str, Any]],
) -> str:
    pop = ancestry["inferred_superpop"]
    md = [
        f"# Local Polygenic Risk Score (PRS) Report — {sample_label.title()}",
        "",
        f"**Run Date:** {ref_meta['run_date']} | **Reference Panel:** 1000 Genomes Phase 3 GRCh38 (unrelated) | "
        f"**Matched Superpopulation:** `{pop}` ({ancestry['n_markers']} ancestry markers)",
        "",
        "> **Clinical & Scientific Disclaimer:** Polygenic Risk Scores aggregate thousands to millions of common variants "
        "to estimate relative genetic predisposition compared to a reference population. They are **not** diagnostic, do **not** "
        "account for lifestyle, family history, or monogenic variants (see `carrier_screening_report.md`), and predictive accuracy "
        "varies by ancestry. Discuss elevated scores with a physician rather than making medical changes independently.",
        "",
        f"## 1. Polygenic Risk Summary (Calibrated vs. 1000 Genomes `{pop}`)",
        "",
        "| Category | Trait | PGS ID | Matched SNPs | Z-Score | Percentile | Risk Tier | Notes |",
        "|---|---|---|---:|---:|---:|---|---|",
    ]
    for r in rows:
        bold = "**" if r["percentile"] >= 80.0 else ""
        md.append(
            f"| {r['category']} | {bold}{r['trait']}{bold} | `{r['pgs_id']}` | "
            f"{r['n_matched_1000g']:,} ({r['match_pct']:.1f}%) | {r['z_score']:+.2f} | "
            f"{bold}{r['percentile']:.1f}%{bold} | {bold}{r['risk_tier']}{bold} | {r['note']} |"
        )

    md.extend([
        "",
        "## 2. How to Interpret These Percentiles",
        f"- **Reference Distribution:** Your exact raw score was compared against unrelated individuals in the 1000 Genomes GRCh38 **`{pop}`** superpopulation scored on the **exact same biallelic SNP set** using PLINK 2.0.",
        "- **Average (21st–79th percentile):** Typical population-level polygenic background.",
        "- **Moderately Elevated (80th–94th percentile):** Top quintile of polygenic risk; worth factoring into standard preventive care (e.g., lipid/BP monitoring, age-appropriate screening).",
        "- **High (≥95th percentile):** Top 5% of polygenic risk; commonly associated with ~1.5x–3x relative risk compared to population median depending on the trait.",
        "- **Sex-Specific Scores:** Breast cancer (`PGS000004`) and Prostate cancer (`PGS000662`) are computed for both samples for completeness but are clinically relevant primarily to the respective biological sex.",
        "- **Alzheimer's (`PGS000334`):** Interpret alongside your *APOE* diplotype in `deep_dive_report.md`.",
        "",
        "## 3. Ancestry Calibration QC",
        f"- **Best-matching 1000G superpopulation:** **`{pop}`** (evaluated across {ancestry['n_markers']} autosomal ancestry-informative SNPs).",
        "- **Superpopulation log-likelihoods:** " + ", ".join(f"`{k}`={v:,.1f}" for k, v in ancestry.get("log_likelihoods", {}).items()),
        "- Full multi-ancestry percentiles (`EUR`, `AFR`, `EAS`, `SAS`, `AMR`) are saved in `polygenic_risk_scores.tsv`.",
    ])
    return "\n".join(md) + "\n"


def run_sample_prs(sample_label: str, override_pop: str | None = None) -> None:
    dist_path = PRS_DIR / "ref_distributions.json"
    combined_score_tsv = PREPARED_DIR / "combined_pgs_alt_scores.tsv"
    ancestry_afs_tsv = PRS_DIR / "ancestry_afs.tsv"

    if not dist_path.exists() or not combined_score_tsv.exists():
        raise FileNotFoundError(
            "PRS reference distributions not built yet. Run `bash scripts/setup_prs_references.sh` first."
        )

    processed_dir = REPO_ROOT / "data" / "processed" / sample_label
    reports_dir = REPO_ROOT / "results" / "reports" / sample_label
    reports_dir.mkdir(parents=True, exist_ok=True)
    personal_vcf = processed_dir / "snp_indel.variants.norm.vcf.gz"

    ref_meta = get_reference_metadata(REFERENCE_DIR)
    distributions = json.loads(dist_path.read_text(encoding="utf-8"))

    print(f"[PRS: {sample_label}] Estimating 1000G reference superpopulation...")
    ancestry = estimate_superpopulation(personal_vcf, ancestry_afs_tsv)
    pop = override_pop.upper() if override_pop else ancestry["inferred_superpop"]
    if override_pop:
        ancestry["inferred_superpop"] = pop
    print(f"[PRS: {sample_label}] Matched superpopulation: {pop} ({ancestry['n_markers']} AIMs)")

    print(f"[PRS: {sample_label}] Scoring {len(PGS_PANEL)} PGS traits across normalized 30x WGS...")
    raw_scores, nonref_hits = score_personal_wgs(personal_vcf, combined_score_tsv, distributions)

    rows: list[dict[str, Any]] = []
    for pid, trait, category, sex_scope, citation in PGS_PANEL:
        smeta = distributions["scores"][pid]
        raw = raw_scores[pid]
        pop_dist = smeta["populations"][pop]
        z = round((raw - pop_dist["mean"]) / max(1e-12, pop_dist["sd"]), 3)
        pct = compute_empirical_percentile(raw, pop_dist["sorted_scores"])
        all_pcts = {
            f"pct_{sp}": compute_empirical_percentile(raw, smeta["populations"][sp]["sorted_scores"])
            for sp in SUPERPOPS
        }
        note = citation
        if sex_scope != "Both":
            note = f"[{sex_scope}-relevant] {citation}"
        rows.append({
            "pgs_id": pid,
            "trait": trait,
            "category": category,
            "sex_scope": sex_scope,
            "n_total_pgs": smeta["n_total_pgs"],
            "n_matched_1000g": smeta["n_matched_1000g"],
            "match_pct": round(100.0 * smeta["n_matched_1000g"] / max(1, smeta["n_total_pgs"]), 1),
            "personal_nonref_snps": nonref_hits[pid],
            "raw_score": round(raw, 6),
            "ref_superpop": pop,
            "z_score": z,
            "percentile": pct,
            "risk_tier": risk_tier_from_percentile(pct),
            **all_pcts,
            "note": note,
        })

    pd.DataFrame(rows).to_csv(reports_dir / "polygenic_risk_scores.tsv", sep="\t", index=False)
    (reports_dir / "polygenic_risk_scores.json").write_text(
        json.dumps({"sample": sample_label, "reference_metadata": ref_meta, "ancestry": ancestry, "scores": rows}, indent=2),
        encoding="utf-8",
    )
    report_md = render_prs_report(sample_label, ref_meta, ancestry, rows)
    out_md = reports_dir / "polygenic_risk_report.md"
    out_md.write_text(report_md, encoding="utf-8")
    print(f"[PRS: {sample_label}] Report written to {out_md}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Gold-standard local PRS calibrated against 1000 Genomes GRCh38.")
    parser.add_argument("--build-ref", action="store_true", help="Build 1000G empirical reference distributions via PLINK 2.0")
    parser.add_argument("--sample", help="Sample label to score (e.g. partner_a, partner_b)")
    parser.add_argument("--superpop", help="Optional override for 1000G reference superpopulation (EUR/AFR/EAS/SAS/AMR)")
    args = parser.parse_args()

    if args.build_ref:
        build_reference_distributions()
    if args.sample:
        run_sample_prs(args.sample, override_pop=args.superpop)
    if not args.build_ref and not args.sample:
        parser.error("Specify --build-ref and/or --sample <label>")


if __name__ == "__main__":
    main()
