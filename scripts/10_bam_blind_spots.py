#!/usr/bin/env python3
"""
Step 10: BAM/CRAM-Based Blind-Spot Callers & Read-Level Verification
Runs specialized short-read alignment callers on data/alignments/<label>/*.{bam,cram}:
  1. ExpansionHunter (v5.x): Repeat expansions (FMR1 Fragile X, HTT, C9ORF72, DMPK, FXN, ATXN*, etc.)
  2. SMNCopyNumberCaller: SMN1/SMN2 exon 7 copy number & g.27134T>G (2+0 silent carrier) status
  3. Cyrius: Full CYP2D6 star-allele, CNV (*5 del, duplications), and CYP2D7 hybrid calling
  4. BAM Read-Level Verification of Exonic SVs (inside vs. flanking depth ratio, MAPQ, split/discordant reads)
  5. BAM Read-Backed Phasing (cis vs. trans via shared reads / read pairs) for candidate compound hets
  6. T1K: HLA class I/II typing (two-field) for CPIC drug-hypersensitivity alleles
  7. HLA disease associations from the same typing: celiac (HLA-DQ2/DQ8) and HLA-B27
  8. UGT1A1 promoter TA repeat (*28 / *36 / *37), counted from reads that span it

Runs only after Step 3 (scripts/03_sample_identity_qc.py) has matched this alignment to the sample's VCF.

Also upgrades results/reports/<label>/pgx_medication_card.{md,tsv} with Cyrius CYP2D6, T1K HLA and UGT1A1 TA-repeat calls.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.genomics import get_reference_metadata, normalize_gt_unphased, read_sample_sex
from lib.hla_traits import interpret_hla_associations
from lib.identity import bam_verified, identity_line, load_identity_qc
from lib.pgx import (
    HLA_PGX_GENES,
    interpret_cyp2d6_cyrius,
    interpret_hla_typing,
    interpret_pgx_profile,
    render_pgx_medication_card,
)
from lib.alignment import caller_cram_problem, is_cram, open_alignment, t1k_input_bam
from lib.ugt1a1 import analyze_ta_repeat, counts_label

REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_TOOLS_DIR = Path.home() / ".local" / "opt" / "wgs_tools"
DEFAULT_BIN_DIR = Path.home() / ".local" / "bin"

# Clinical repeat-size thresholds (repeat units) for key ExpansionHunter GRCh38 loci.
# Sources: ACMG / GeneReviews / Illumina ExpansionHunter clinical reference thresholds.
REPEAT_THRESHOLDS: dict[str, dict[str, Any]] = {
    "FMR1": {"condition": "Fragile X Syndrome / FXTAS / FXPOI", "moi": "XL", "normal_max": 44, "intermediate_max": 54, "premutation_max": 200, "pathogenic_min": 201},
    "HTT": {"condition": "Huntington Disease", "moi": "AD", "normal_max": 26, "intermediate_max": 35, "pathogenic_min": 36},
    "C9ORF72": {"condition": "ALS / Frontotemporal Dementia (C9orf72)", "moi": "AD", "normal_max": 24, "intermediate_max": 29, "pathogenic_min": 30},
    "DMPK": {"condition": "Myotonic Dystrophy Type 1 (DM1)", "moi": "AD", "normal_max": 34, "intermediate_max": 49, "pathogenic_min": 50},
    "FXN": {"condition": "Friedreich Ataxia", "moi": "AR", "normal_max": 33, "intermediate_max": 65, "pathogenic_min": 66},
    "ATXN1": {"condition": "Spinocerebellar Ataxia Type 1 (SCA1)", "moi": "AD", "normal_max": 35, "intermediate_max": 38, "pathogenic_min": 39},
    "ATXN2": {"condition": "Spinocerebellar Ataxia Type 2 (SCA2)", "moi": "AD", "normal_max": 31, "intermediate_max": 32, "pathogenic_min": 33},
    "ATXN3": {"condition": "Spinocerebellar Ataxia Type 3 (SCA3)", "moi": "AD", "normal_max": 44, "intermediate_max": 59, "pathogenic_min": 60},
    "ATXN7": {"condition": "Spinocerebellar Ataxia Type 7 (SCA7)", "moi": "AD", "normal_max": 19, "intermediate_max": 36, "pathogenic_min": 37},
    "ATXN8OS": {"condition": "Spinocerebellar Ataxia Type 8 (SCA8)", "moi": "AD", "normal_max": 50, "intermediate_max": 70, "pathogenic_min": 71},
    "ATXN10": {"condition": "Spinocerebellar Ataxia Type 10 (SCA10)", "moi": "AD", "normal_max": 32, "intermediate_max": 279, "pathogenic_min": 280},
    "CACNA1A": {"condition": "Spinocerebellar Ataxia Type 6 (SCA6)", "moi": "AD", "normal_max": 18, "intermediate_max": 19, "pathogenic_min": 20},
    "AR": {"condition": "Spinal & Bulbar Muscular Atrophy (Kennedy Disease)", "moi": "XL", "normal_max": 34, "intermediate_max": 37, "pathogenic_min": 38},
    "ATN1": {"condition": "Dentatorubral-Pallidoluysian Atrophy (DRPLA)", "moi": "AD", "normal_max": 35, "intermediate_max": 47, "pathogenic_min": 48},
    "CNBP": {"condition": "Myotonic Dystrophy Type 2 (DM2)", "moi": "AD", "normal_max": 26, "intermediate_max": 74, "pathogenic_min": 75},
    "CSTB": {"condition": "Unverricht-Lundborg Disease (EPM1)", "moi": "AR", "normal_max": 3, "intermediate_max": 29, "pathogenic_min": 30},
    "JPH3": {"condition": "Huntington Disease-Like 2 (HDL2)", "moi": "AD", "normal_max": 28, "intermediate_max": 39, "pathogenic_min": 40},
    "PABPN1": {"condition": "Oculopharyngeal Muscular Dystrophy (OPMD)", "moi": "AD/AR", "normal_max": 10, "intermediate_max": 10, "pathogenic_min": 11},
    "PHOX2B": {"condition": "Congenital Central Hypoventilation Syndrome", "moi": "AD", "normal_max": 20, "intermediate_max": 23, "pathogenic_min": 24},
    "PPP2R2B": {"condition": "Spinocerebellar Ataxia Type 12 (SCA12)", "moi": "AD", "normal_max": 32, "intermediate_max": 42, "pathogenic_min": 43},
    "TBP": {"condition": "Spinocerebellar Ataxia Type 17 (SCA17)", "moi": "AD", "normal_max": 40, "intermediate_max": 42, "pathogenic_min": 43},
    "TCF4": {"condition": "Fuchs Endothelial Corneal Dystrophy (FECD3)", "moi": "AD/Risk", "normal_max": 40, "intermediate_max": 49, "pathogenic_min": 50},
    "AFF2": {"condition": "FRAXE Intellectual Disability", "moi": "XL", "normal_max": 30, "intermediate_max": 199, "pathogenic_min": 200},
    "RFC1": {"condition": "CANVAS / Late-Onset Ataxia (note: motif-dependent)", "moi": "AR", "normal_max": 100, "intermediate_max": 399, "pathogenic_min": 400},
    "NOP56": {"condition": "Spinocerebellar Ataxia Type 36 (SCA36)", "moi": "AD", "normal_max": 14, "intermediate_max": 649, "pathogenic_min": 650},
    "GLS": {"condition": "Global Developmental Delay / Ataxia (GLS)", "moi": "AR", "normal_max": 90, "intermediate_max": 399, "pathogenic_min": 400},
    "DIP2B": {"condition": "FRA12A Intellectual Disability", "moi": "AD", "normal_max": 23, "intermediate_max": 249, "pathogenic_min": 250},
    "CBL": {"condition": "Jacobsen Syndrome (FRA11B)", "moi": "AD", "normal_max": 35, "intermediate_max": 99, "pathogenic_min": 100},
    "NIPA1": {"condition": "ALS Susceptibility Modifier (NIPA1)", "moi": "Risk", "normal_max": 10, "intermediate_max": 10, "pathogenic_min": 11},
    "GIPC1": {"condition": "Oculopharyngodistal Myopathy 2 (OPDM2)", "moi": "AD", "normal_max": 32, "intermediate_max": 72, "pathogenic_min": 73},
    "NOTCH2NL": {"condition": "Neuronal Intranuclear Inclusion Disease (NIID)", "moi": "AD", "normal_max": 41, "intermediate_max": 64, "pathogenic_min": 65},
}

PRIORITY_REPEAT_LOCI = ("FMR1", "HTT", "C9ORF72", "DMPK", "FXN", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A", "AR", "RFC1")


def classify_repeat_size(locus: str, genotype: str) -> tuple[str, int | None]:
    """Classify an ExpansionHunter genotype string (e.g. '29/31' or '30') against clinical thresholds."""
    if not genotype or genotype in (".", "./."):
        return "NO_CALL", None
    alleles: list[int] = []
    for part in str(genotype).split("/"):
        part = part.strip()
        if part.isdigit():
            alleles.append(int(part))
    if not alleles:
        return "NO_CALL", None
    max_rep = max(alleles)
    thresh = REPEAT_THRESHOLDS.get(locus)
    if not thresh:
        return "ASSESSED", max_rep
    if locus == "FMR1":
        if max_rep <= thresh["normal_max"]:
            return "NORMAL", max_rep
        if max_rep <= thresh["intermediate_max"]:
            return "INTERMEDIATE_GRAY_ZONE", max_rep
        if max_rep <= thresh["premutation_max"]:
            return "PREMUTATION", max_rep
        return "FULL_MUTATION", max_rep
    if max_rep <= thresh["normal_max"]:
        return "NORMAL", max_rep
    if max_rep < thresh["pathogenic_min"]:
        return "INTERMEDIATE", max_rep
    return "EXPANDED_PATHOGENIC_RANGE", max_rep


def find_alignment_file(sample_label: str) -> Path | None:
    """Locate the first BAM/CRAM file (regular file or symlink) in data/alignments/<label>/."""
    align_dir = REPO_ROOT / "data" / "alignments" / sample_label
    if not align_dir.exists():
        return None
    for ext in ("*.bam", "*.cram"):
        matches = sorted(align_dir.glob(ext))
        if matches:
            return matches[0]
    return None


def get_sample_sex(sample_label: str) -> str:
    """Look up sample sex ('M' or 'F') from config/samples.tsv."""
    return read_sample_sex(REPO_ROOT / "config" / "samples.tsv", sample_label)


def detect_bam_has_chr(bam_path: Path) -> bool:
    with open_alignment(bam_path) as af:
        refs = set(af.references[:30])
    return "chr1" in refs


def find_tool(name: str) -> str | None:
    found = shutil.which(name)
    if found:
        return found
    candidate = DEFAULT_BIN_DIR / name
    if candidate.exists() and os.access(candidate, os.X_OK):
        return str(candidate)
    return None


def run_expansion_hunter(
    bam_path: Path,
    fasta_path: Path,
    sex: str,
    out_dir: Path,
    sample_label: str,
    threads: int,
    has_chr: bool,
) -> dict[str, Any]:
    """Run ExpansionHunter v5.x and parse all repeat loci."""
    eh_bin = find_tool("ExpansionHunter")
    if not eh_bin:
        return {"status": "TOOL_NOT_FOUND", "error": "ExpansionHunter not found in PATH or ~/.local/bin"}
    if not fasta_path.exists():
        return {"status": "MISSING_REFERENCE", "error": f"Reference FASTA not found: {fasta_path}"}
    cram_problem = caller_cram_problem(bam_path)
    if cram_problem:
        return {"status": "UNSUPPORTED_CRAM", "error": cram_problem}

    tools_dir = Path(os.environ.get("WGS_TOOLS_DIR", str(DEFAULT_TOOLS_DIR)))
    eh_cat_dir = tools_dir / "ExpansionHunter-v5.0.0-linux_x86_64" / "variant_catalog" / "grch38"
    catalog = eh_cat_dir / ("variant_catalog.json" if has_chr else "variant_catalog.nochr.json")
    if not catalog.exists():
        return {"status": "MISSING_CATALOG", "error": f"Variant catalog not found: {catalog}"}

    eh_dir = out_dir / "expansion_hunter"
    eh_dir.mkdir(parents=True, exist_ok=True)
    prefix = eh_dir / sample_label
    out_json = Path(f"{prefix}.json")

    cmd = [
        eh_bin,
        "--reads", str(bam_path.resolve()),
        "--reference", str(fasta_path.resolve()),
        "--variant-catalog", str(catalog),
        "--output-prefix", str(prefix),
        "--threads", str(threads),
    ]
    if sex == "M":
        cmd.extend(["--sex", "male"])
    elif sex == "F":
        cmd.extend(["--sex", "female"])

    if not out_json.exists():
        print(f"[{sample_label}] Running ExpansionHunter ({len(REPEAT_THRESHOLDS)} GRCh38 repeat loci)...")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            return {
                "status": "ERROR",
                "error": (proc.stderr or proc.stdout or f"exit {proc.returncode}")[-500:],
            }
    else:
        print(f"[{sample_label}] Reusing existing ExpansionHunter output: {out_json}")

    with open(out_json, encoding="utf-8") as f:
        raw = json.load(f)

    loci_rows: list[dict[str, Any]] = []
    flagged: list[dict[str, Any]] = []
    for locus_id, ldata in sorted(raw.get("LocusResults", {}).items()):
        variants = ldata.get("Variants", {})
        # Primary repeat variant (usually key == locus_id)
        vdata = variants.get(locus_id) or next(iter(variants.values()), {})
        gt = str(vdata.get("Genotype") or "./.")
        ci = str(vdata.get("GenotypeConfidenceInterval") or "")
        unit = str(vdata.get("RepeatUnit") or "")
        ref_region = str(vdata.get("ReferenceRegion") or "")
        cls, max_rep = classify_repeat_size(locus_id, gt)
        thresh = REPEAT_THRESHOLDS.get(locus_id, {})
        row = {
            "locus": locus_id,
            "condition": thresh.get("condition", ""),
            "moi": thresh.get("moi", ""),
            "repeat_unit": unit,
            "region": ref_region,
            "genotype": gt,
            "ci": ci,
            "max_repeat": max_rep,
            "normal_max": thresh.get("normal_max"),
            "pathogenic_min": thresh.get("pathogenic_min"),
            "classification": cls,
            "coverage": round(float(ldata.get("Coverage") or 0.0), 1),
        }
        loci_rows.append(row)
        if cls not in ("NORMAL", "ASSESSED"):
            flagged.append(row)

    return {
        "status": "OK",
        "tool": "ExpansionHunter v5.0.0",
        "sex_used": raw.get("SampleParameters", {}).get("Sex", sex),
        "n_loci": len(loci_rows),
        "n_flagged": len(flagged),
        "loci": loci_rows,
        "flagged": flagged,
        "fmr1": next((r for r in loci_rows if r["locus"] == "FMR1"), None),
        "output_json": str(out_json),
    }


def parse_smn_results(tsv_path: Path, json_path: Path) -> dict[str, Any]:
    """Parse SMNCopyNumberCaller TSV + JSON outputs into a structured summary."""
    row: dict[str, Any] = {}
    with open(tsv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        row = next(reader, {})

    raw_json: dict[str, Any] = {}
    if json_path.exists():
        with open(json_path, encoding="utf-8") as f:
            j = json.load(f)
            if isinstance(j, dict) and j:
                raw_json = next(iter(j.values())) if not ("SMN1_CN" in j) else j

    smn1_cn = row.get("SMN1_CN")
    smn2_cn = row.get("SMN2_CN")
    smn2_del = row.get("SMN2delta7-8_CN")
    is_sma = str(row.get("isSMA", "")).strip()
    is_carrier = str(row.get("isCarrier", "")).strip()
    g27134_cn = row.get("g.27134T>G_CN")
    if g27134_cn is None:
        g27134_cn = raw_json.get("g27134TG_CN")
    ref_c = raw_json.get("g27134TG_REF_count")
    alt_c = raw_json.get("g27134TG_ALT_count")
    if g27134_cn is not None and str(g27134_cn) != "":
        detail = f" (REF={ref_c[0]}, ALT={alt_c[0]})" if (ref_c and alt_c) else ""
        g27134 = f"CN={g27134_cn}{detail}"
    else:
        g27134 = "Not reported"

    # Build plain-language clinical interpretation
    if is_sma.lower() == "true":
        interp = "AFFECTED_RANGE (SMN1 CN=0; clinical confirmation urgently required)"
    elif is_carrier.lower() == "true":
        interp = f"SMA_CARRIER (SMN1 CN={smn1_cn}; confirm clinically)"
    elif str(smn1_cn) in ("2", "3", "4", "5") and is_carrier.lower() == "false":
        interp = f"NON_CARRIER_LOW_RISK (SMN1 CN={smn1_cn}, SMN2 CN={smn2_cn}; silent 2+0 SNP negative; low residual risk)"
    else:
        interp = f"INDETERMINATE (SMN1 CN={smn1_cn}, isCarrier={is_carrier}; confirm via clinical MLPA/qPCR)"

    return {
        "status": "OK",
        "tool": "Illumina SMNCopyNumberCaller",
        "is_sma": is_sma,
        "is_carrier": is_carrier,
        "smn1_cn": smn1_cn,
        "smn2_cn": smn2_cn,
        "smn2_delta7_8_cn": smn2_del,
        "g27134tg": str(g27134),
        "interpretation": interp,
        "raw_row": row,
        "output_tsv": str(tsv_path),
    }


def run_smn_caller(
    bam_path: Path, out_dir: Path, sample_label: str, threads: int, fasta_path: Path | None = None
) -> dict[str, Any]:
    smn_bin = find_tool("smn_caller.py")
    if not smn_bin:
        return {"status": "TOOL_NOT_FOUND", "error": "smn_caller.py not found in PATH or ~/.local/bin"}

    smn_dir = out_dir / "smn_caller"
    smn_dir.mkdir(parents=True, exist_ok=True)
    manifest = smn_dir / "manifest.txt"
    manifest.write_text(f"{bam_path.resolve()}\n", encoding="utf-8")

    tsv_path = smn_dir / f"{sample_label}.tsv"
    json_path = smn_dir / f"{sample_label}.json"

    if not tsv_path.exists():
        print(f"[{sample_label}] Running SMNCopyNumberCaller (GRCh38)...")
        cmd = [
            smn_bin,
            "--manifest", str(manifest),
            "--genome", "38",
            "--prefix", sample_label,
            "--outDir", str(smn_dir),
            "--threads", str(threads),
        ]
        if is_cram(bam_path) and fasta_path is not None and fasta_path.exists():
            cmd += ["--reference", str(fasta_path.resolve())]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0 or not tsv_path.exists():
            return {
                "status": "ERROR",
                "error": (proc.stderr or proc.stdout or f"exit {proc.returncode}")[-500:],
            }
    else:
        print(f"[{sample_label}] Reusing existing SMNCopyNumberCaller output: {tsv_path}")

    return parse_smn_results(tsv_path, json_path)


def run_cyrius(
    bam_path: Path, out_dir: Path, sample_label: str, threads: int, fasta_path: Path | None = None
) -> dict[str, Any]:
    cyrius_bin = find_tool("star_caller.py")
    if not cyrius_bin:
        return {"status": "TOOL_NOT_FOUND", "error": "star_caller.py not found in PATH or ~/.local/bin"}

    cyp_dir = out_dir / "cyrius"
    cyp_dir.mkdir(parents=True, exist_ok=True)
    manifest = cyp_dir / "manifest.txt"
    manifest.write_text(f"{bam_path.resolve()}\n", encoding="utf-8")

    tsv_path = cyp_dir / f"{sample_label}.tsv"
    json_path = cyp_dir / f"{sample_label}.json"

    if not tsv_path.exists():
        print(f"[{sample_label}] Running Cyrius (CYP2D6 star-allele / CNV / hybrid caller, GRCh38)...")
        cmd = [
            cyrius_bin,
            "--manifest", str(manifest),
            "--genome", "38",
            "--prefix", sample_label,
            "--outDir", str(cyp_dir),
            "--threads", str(threads),
        ]
        if is_cram(bam_path) and fasta_path is not None and fasta_path.exists():
            cmd += ["--reference", str(fasta_path.resolve())]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0 or not tsv_path.exists():
            return {
                "status": "ERROR",
                "error": (proc.stderr or proc.stdout or f"exit {proc.returncode}")[-500:],
            }
    else:
        print(f"[{sample_label}] Reusing existing Cyrius output: {tsv_path}")

    row: dict[str, Any] = {}
    with open(tsv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        row = next(reader, {})

    raw_json: dict[str, Any] = {}
    if json_path.exists():
        with open(json_path, encoding="utf-8") as f:
            j = json.load(f)
            if isinstance(j, dict) and j:
                raw_json = next(iter(j.values())) if "Genotype" not in j else j

    gt = row.get("Genotype") or raw_json.get("Genotype")
    flt = row.get("Filter") or raw_json.get("Filter")
    cpic = interpret_cyp2d6_cyrius(gt, flt)

    return {
        "status": "OK",
        "tool": "Illumina Cyrius",
        "genotype": gt,
        "filter": flt,
        "cnv_group": raw_json.get("CNV_group"),
        "total_cn": raw_json.get("Total_CN"),
        "total_cn_raw": raw_json.get("Total_CN_raw"),
        "spacer_cn": raw_json.get("Spacer_CN"),
        "exon9_cn": raw_json.get("Exon9_CN"),
        "cpic_phenotype": cpic["phenotype"],
        "cpic_tier": cpic["tier"],
        "cpic_activity_score": cpic.get("activity_score"),
        "cpic_guidance": cpic["guidance"],
        "output_tsv": str(tsv_path),
    }


HLA_REPORT_GENES: tuple[str, ...] = ("HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPB1")


def parse_t1k_genotype(tsv_path: Path) -> dict[str, Any]:
    """
    Parse a T1K *_genotype.tsv (no header):
    gene_name num_diff_alleles allele_1 abundance_1 quality_1 allele_2 abundance_2 quality_2 secondary_alleles
    Alleles with quality <= 0 are dropped (T1K README recommendation); '.' marks a missing/homozygous slot.
    """
    genotypes: dict[str, list[str]] = {}
    qualities: dict[str, list[float]] = {}
    rows: list[dict[str, Any]] = []
    with open(tsv_path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            gene = cols[0]
            called: list[str] = []
            called_q: list[float] = []
            for allele, qual in ((cols[2], cols[4]), (cols[5], cols[7])):
                try:
                    q = float(qual)
                except ValueError:
                    q = -1.0
                if allele and allele != "." and q > 0:
                    called.append(allele)  # may be a comma-separated list of tied candidates
                    called_q.append(q)
            genotypes[gene] = called
            qualities[gene] = called_q
            rows.append({
                "gene": gene,
                "num_diff_alleles": cols[1],
                "allele_1": cols[2],
                "quality_1": cols[4],
                "allele_2": cols[5],
                "quality_2": cols[7],
                "called": "/".join(called) if called else "no confident call",
            })
    return {"genotypes": genotypes, "qualities": qualities, "rows": rows}


def _space_free_alignment(bam_path: Path, link_dir: Path) -> Path:
    """Return a path to the BAM/CRAM with no spaces, symlinking it (and its index) into link_dir if needed.

    run-t1k hands its arguments to a shell unquoted, so a path such as
    "/mnt/d/My Data/x.bam" splits at the space and bam-extractor crashes.
    """
    target = bam_path.resolve()
    if " " not in str(target):
        return target
    link = link_dir / f"input{target.suffix}"
    for index_suffix in ("", ".bai", ".crai", ".csi"):
        src = Path(str(target) + index_suffix)
        dst = Path(str(link) + index_suffix)
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        if src.exists():
            dst.symlink_to(src)
    # Also cover "x.bai" (index named without the .bam part)
    alt_index = target.with_suffix(".bai")
    if alt_index.exists() and not Path(str(link) + ".bai").exists():
        Path(str(link) + ".bai").symlink_to(alt_index)
    return link


def run_t1k_hla(
    bam_path: Path, out_dir: Path, sample_label: str, threads: int, has_chr: bool, fasta_path: Path | None = None
) -> dict[str, Any]:
    """Run T1K (--preset hla-wgs) on the BAM for two-field HLA typing."""
    t1k_bin = find_tool("run-t1k")
    if not t1k_bin:
        return {"status": "TOOL_NOT_FOUND", "error": "run-t1k not found in PATH or ~/.local/bin (run scripts/setup_bam_tools.sh)"}

    tools_dir = Path(os.environ.get("WGS_TOOLS_DIR", str(DEFAULT_TOOLS_DIR)))
    idx_dir = tools_dir / "t1k_hlaidx"
    seq_fa = idx_dir / "hlaidx_dna_seq.fa"
    coord_fa = idx_dir / ("hlaidx_dna_coord.fa" if has_chr else "hlaidx_dna_coord.nochr.fa")
    if not seq_fa.exists() or not coord_fa.exists():
        return {"status": "MISSING_INDEX", "error": f"T1K HLA index not found in {idx_dir} (run scripts/setup_bam_tools.sh)"}

    t1k_dir = out_dir / "t1k"
    t1k_dir.mkdir(parents=True, exist_ok=True)
    geno_tsv = t1k_dir / f"{sample_label}_genotype.tsv"

    if not geno_tsv.exists():
        print(f"[{sample_label}] Running T1K HLA typing (hla-wgs preset)...")
        t1k_bam = bam_path
        if is_cram(bam_path):
            # T1K's bam-extractor cannot read CRAM: give it a BAM of the reads its whole-file scan would pick
            print(f"[{sample_label}] Extracting the T1K input reads from the CRAM into a temporary BAM...")
            try:
                t1k_bam = t1k_input_bam(bam_path, coord_fa, t1k_dir, fasta_path, threads)
            except RuntimeError as e:
                return {"status": "ERROR", "error": str(e)}
        cmd = [
            t1k_bin,
            "-b", str(_space_free_alignment(t1k_bam, t1k_dir)),
            "-f", str(seq_fa),
            "-c", str(coord_fa),
            "--preset", "hla-wgs",
            "-t", str(threads),
            "-o", sample_label,
            "--od", str(t1k_dir),
            "--alleleDigitUnits", "2",
            "--alleleDelimiter", ":",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0 or not geno_tsv.exists():
            return {
                "status": "ERROR",
                "error": (proc.stderr or proc.stdout or f"exit {proc.returncode}")[-500:],
            }
    else:
        print(f"[{sample_label}] Reusing existing T1K output: {geno_tsv}")

    parsed = parse_t1k_genotype(geno_tsv)
    return {
        "status": "OK",
        "tool": "T1K (--preset hla-wgs)",
        "genotypes": parsed["genotypes"],
        "qualities": parsed["qualities"],
        "rows": parsed["rows"],
        "output_tsv": str(geno_tsv),
    }


def _mean_depth_and_mapq(
    af: pysam.AlignmentFile, chrom: str, start: int, end: int
) -> tuple[float, float]:
    """Compute mean read depth and mean MAPQ over [start, end) using sampled windows for large spans."""
    length = max(1, end - start)
    windows: list[tuple[int, int]] = []
    if length <= 10_000:
        windows = [(start, end)]
    else:
        # Sample 10 evenly spaced 1 kb windows inside large SVs for fast O(ms) depth estimation
        step = max(1, (length - 1000) // 9)
        for i in range(10):
            ws = start + i * step
            we = min(end, ws + 1000)
            if we > ws:
                windows.append((ws, we))

    total_bases = 0
    total_cov = 0
    mapq_sum = 0
    mapq_n = 0
    for ws, we in windows:
        wlen = we - ws
        total_bases += wlen
        # Fast approximate depth from aligned read overlap length (avoids slow per-base pileup)
        for r in af.fetch(chrom, ws, we):
            if r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary:
                continue
            ov = r.get_overlap(ws, we)
            total_cov += ov
            mapq_sum += r.mapping_quality
            mapq_n += 1
    mean_dp = total_cov / total_bases if total_bases else 0.0
    mean_mq = mapq_sum / mapq_n if mapq_n else 0.0
    return round(mean_dp, 1), round(mean_mq, 1)


def _breakpoint_evidence(
    af: pysam.AlignmentFile, chrom: str, bp: int
) -> tuple[int, int]:
    """Count soft/hard-clipped reads within ±50 bp and discordant read pairs within ±500 bp of bp."""
    clipped = 0
    discordant = 0
    lo = max(0, bp - 500)
    hi = bp + 500
    for r in af.fetch(chrom, lo, hi):
        if r.is_unmapped or r.is_duplicate or r.is_secondary:
            continue
        if r.is_paired and not r.is_proper_pair:
            discordant += 1
        if abs(r.reference_start - bp) <= 50 or (r.reference_end and abs(r.reference_end - bp) <= 50):
            if r.cigartuples and any(op in (4, 5) and ln >= 10 for op, ln in r.cigartuples):
                clipped += 1
    return clipped, discordant


def classify_sv_bam_support(
    svtype: str,
    depth_ratio: float,
    mapq_inside: float,
    clipped_total: int,
    discordant_total: int,
) -> str:
    """Classify read-level BAM support for a short-read CNV/SV call."""
    st = svtype.upper()
    if st in ("DEL", "LOSS"):
        if depth_ratio <= 0.15:
            return "CONFIRMED_HOM_DEL"
        if depth_ratio <= 0.72:
            return "CONFIRMED_HET_DEL"
        if mapq_inside < 20.0:
            return "LOW_MAPQ_COMPLEX_REGION"
        return "NO_DEPTH_DROP_LIKELY_ARTIFACT_OR_REPEAT"
    if st in ("DUP", "GAIN"):
        if depth_ratio >= 1.30:
            return "CONFIRMED_DUP_GAIN"
        if mapq_inside < 20.0:
            return "LOW_MAPQ_COMPLEX_REGION"
        return "NO_DEPTH_GAIN_LIKELY_ARTIFACT_OR_REPEAT"
    if st == "BND":
        if clipped_total >= 5 or discordant_total >= 10:
            return "BREAKPOINT_READ_SUPPORT"
        return "WEAK_BREAKPOINT_SUPPORT"
    return "INSPECTED"


def verify_exonic_svs_in_bam(
    bam_path: Path, sv_tsv: Path, has_chr: bool, reference: Path | None = None
) -> list[dict[str, Any]]:
    """Inspect BAM depth ratio, MAPQ, clipped reads, and discordant pairs for exonic SVs."""
    if not sv_tsv.exists():
        return []
    df = pd.read_csv(sv_tsv, sep="\t")
    if df.empty:
        return []

    rows: list[dict[str, Any]] = []
    with open_alignment(bam_path, reference) as af:
        ref_lens = dict(zip(af.references, af.lengths))
        for _, rec in df.iterrows():
            raw_chrom = str(rec["chrom"])
            chrom = f"chr{raw_chrom}" if (has_chr and not raw_chrom.startswith("chr")) else raw_chrom.removeprefix("chr")
            if chrom not in ref_lens:
                continue
            clen = ref_lens[chrom]
            start = max(0, int(rec["start"]))
            end = min(clen, int(rec["end"]))
            if end <= start:
                end = min(clen, start + 1)
            svtype = str(rec["svtype"])

            in_dp, in_mq = _mean_depth_and_mapq(af, chrom, start, end)
            l_start, l_end = max(0, start - 6000), max(0, start - 1000)
            r_start, r_end = min(clen, end + 1000), min(clen, end + 6000)
            l_dp, _ = _mean_depth_and_mapq(af, chrom, l_start, l_end) if l_end > l_start else (0.0, 0.0)
            r_dp, _ = _mean_depth_and_mapq(af, chrom, r_start, r_end) if r_end > r_start else (0.0, 0.0)
            flanks = [d for d in (l_dp, r_dp) if d > 0]
            flank_dp = round(sum(flanks) / len(flanks), 1) if flanks else 0.0
            ratio = round(in_dp / flank_dp, 3) if flank_dp > 0 else 0.0

            c1, d1 = _breakpoint_evidence(af, chrom, start)
            c2, d2 = _breakpoint_evidence(af, chrom, end) if end - start > 200 else (0, 0)
            clipped_total = c1 + c2
            disc_total = d1 + d2

            verdict = classify_sv_bam_support(svtype, ratio, in_mq, clipped_total, disc_total)
            rows.append({
                "source": str(rec["source"]),
                "chrom": raw_chrom,
                "start": start,
                "end": end,
                "length_bp": int(rec["length_bp"]),
                "svtype": svtype,
                "gt_or_cn": str(rec["gt_or_cn"]),
                "cds_genes": str(rec.get("cds_genes") or ""),
                "inside_dp": in_dp,
                "flank_dp": flank_dp,
                "depth_ratio": ratio,
                "mapq_inside": in_mq,
                "clipped_reads": clipped_total,
                "discordant_pairs": disc_total,
                "bam_verdict": verdict,
            })
    return rows


def check_pair_phasing_in_bam(
    af: pysam.AlignmentFile,
    chrom: str,
    pos1: int,
    ref1: str,
    alt1: str,
    pos2: int,
    ref2: str,
    alt2: str,
) -> dict[str, Any]:
    """
    Check cis vs. trans phasing between two 1-based SNV sites using shared reads / read pairs.
    """
    p1, p2 = pos1 - 1, pos2 - 1
    dist = abs(pos2 - pos1)

    def _base_at(read: pysam.AlignedSegment, target_0: int) -> str | None:
        if read.reference_start is None or read.reference_end is None:
            return None
        if not (read.reference_start <= target_0 < read.reference_end):
            return None
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            if rpos == target_0 and qpos is not None and read.query_sequence:
                return read.query_sequence[qpos].upper()
        return None

    q1: dict[str, str] = {}
    for r in af.fetch(chrom, p1, p1 + 1):
        if r.is_unmapped or r.is_duplicate or r.is_secondary or r.mapping_quality < 10:
            continue
        b = _base_at(r, p1)
        if b in (ref1.upper(), alt1.upper()):
            q1[r.query_name] = b

    q2: dict[str, str] = {}
    for r in af.fetch(chrom, p2, p2 + 1):
        if r.is_unmapped or r.is_duplicate or r.is_secondary or r.mapping_quality < 10:
            continue
        b = _base_at(r, p2)
        if b in (ref2.upper(), alt2.upper()):
            q2[r.query_name] = b

    shared = set(q1) & set(q2)
    cis_n = 0
    trans_n = 0
    for qn in shared:
        a1_is_alt = q1[qn] == alt1.upper()
        a2_is_alt = q2[qn] == alt2.upper()
        if a1_is_alt == a2_is_alt:
            cis_n += 1
        else:
            trans_n += 1

    if not shared:
        verdict = "UNPHASED_EXCEEDS_INSERT_SIZE" if dist > 700 else "NO_SPANNING_READS"
    elif cis_n >= 2 and trans_n == 0:
        verdict = "CIS"
    elif trans_n >= 2 and cis_n == 0:
        verdict = "TRANS"
    elif cis_n > trans_n * 3:
        verdict = "LIKELY_CIS"
    elif trans_n > cis_n * 3:
        verdict = "LIKELY_TRANS"
    else:
        verdict = "AMBIGUOUS"

    return {
        "distance_bp": dist,
        "shared_read_pairs": len(shared),
        "cis_supporting": cis_n,
        "trans_supporting": trans_n,
        "phase_verdict": verdict,
    }


def check_candidate_phasing(
    bam_path: Path, rep_dir: Path, has_chr: bool, reference: Path | None = None
) -> list[dict[str, Any]]:
    """Check read-backed phasing for any candidate compound hets + double-het MTHFR."""
    pairs_to_check: list[dict[str, Any]] = []

    # 1. Any same-gene compound het candidates from Step 4/5
    comphet_tsv = rep_dir / "snp_indel_hom_hemi_comphet.tsv"
    if comphet_tsv.exists() and comphet_tsv.stat().st_size > 10:
        try:
            cdf = pd.read_csv(comphet_tsv, sep="\t")
            if not cdf.empty and "gene" in cdf.columns:
                for gene, grp in cdf.groupby("gene"):
                    hets = grp[grp["zygosity"] == "HET"]
                    if len(hets) >= 2:
                        recs = hets.to_dict("records")
                        for i in range(len(recs)):
                            for j in range(i + 1, len(recs)):
                                pairs_to_check.append({
                                    "category": "CLINICAL_COMPHET",
                                    "gene": str(gene),
                                    "chrom": str(recs[i]["chrom"]),
                                    "var1": f"{recs[i]['variant_id']} ({recs[i].get('rsid', '.')})",
                                    "pos1": int(recs[i]["pos"]),
                                    "ref1": str(recs[i]["ref"]),
                                    "alt1": str(recs[i]["alt"]),
                                    "var2": f"{recs[j]['variant_id']} ({recs[j].get('rsid', '.')})",
                                    "pos2": int(recs[j]["pos"]),
                                    "ref2": str(recs[j]["ref"]),
                                    "alt2": str(recs[j]["alt"]),
                                })
        except pd.errors.EmptyDataError:
            pass

    # 2. Check MTHFR C677T (rs1801133) + A1298C (rs1801131) if both are HET
    tm_tsv = rep_dir / "targeted_markers.tsv"
    if tm_tsv.exists():
        tdf = pd.read_csv(tm_tsv, sep="\t")
        by_rs = {str(r["rsid"]): r for _, r in tdf.iterrows()}
        r1, r2 = by_rs.get("rs1801133"), by_rs.get("rs1801131")
        if r1 is not None and r2 is not None:
            if normalize_gt_unphased(str(r1["gt"])) == "0/1" and normalize_gt_unphased(str(r2["gt"])) == "0/1":
                # locus format: 1:11796321 G>A
                def _parse_locus(loc: str) -> tuple[str, int, str, str]:
                    coord, alleles = loc.split()
                    c, p = coord.split(":")
                    ref, alt = alleles.split(">")
                    return c, int(p), ref, alt

                c1, p1, ref1, alt1 = _parse_locus(str(r1["locus"]))
                _, p2, ref2, alt2 = _parse_locus(str(r2["locus"]))
                pairs_to_check.append({
                    "category": "COMMON_TRAIT_DOUBLE_HET",
                    "gene": "MTHFR",
                    "chrom": c1,
                    "var1": f"rs1801133 ({r1['locus']})",
                    "pos1": p1, "ref1": ref1, "alt1": alt1,
                    "var2": f"rs1801131 ({r2['locus']})",
                    "pos2": p2, "ref2": ref2, "alt2": alt2,
                })

    results: list[dict[str, Any]] = []
    if not pairs_to_check:
        return results

    with open_alignment(bam_path, reference) as af:
        for pair in pairs_to_check:
            raw_c = pair["chrom"]
            chrom = f"chr{raw_c}" if (has_chr and not raw_c.startswith("chr")) else raw_c.removeprefix("chr")
            ph = check_pair_phasing_in_bam(
                af, chrom,
                pair["pos1"], pair["ref1"], pair["alt1"],
                pair["pos2"], pair["ref2"], pair["alt2"],
            )
            results.append({**pair, **ph})
    return results


def render_bam_report(
    sample_label: str,
    bam_path: Path,
    eh_res: dict[str, Any],
    smn_res: dict[str, Any],
    cyrius_res: dict[str, Any],
    sv_rows: list[dict[str, Any]],
    phasing_rows: list[dict[str, Any]],
    hla_res: dict[str, Any] | None = None,
    ugt_res: dict[str, Any] | None = None,
    identity_qc: dict[str, Any] | None = None,
) -> str:
    today = datetime.date.today().isoformat()
    md = [
        f"# BAM-Based Blind-Spot & Read-Level Verification Report — {sample_label.title()}",
        "",
        f"**Run Date:** {today} | **Alignment:** `{bam_path.name}`",
        "",
        identity_line(identity_qc) + " (details: `sample_identity_qc.md`)",
        "",
        "> **Clinical Disclaimer:** Research-grade BAM bioinformatics analysis (ExpansionHunter, SMNCopyNumberCaller, Cyrius, T1K). "
        "Confirm any actionable or reproductive finding via a CLIA-certified clinical laboratory before medical decisions.",
        "",
        "## 1. Spinal Muscular Atrophy (*SMN1* / *SMN2* Exon 7 & 2+0 Silent Carrier)",
    ]
    if smn_res.get("status") == "OK":
        md.extend([
            "| Metric | Result |",
            "|---|---|",
            f"| **Caller** | {smn_res['tool']} |",
            f"| **SMN1 Copy Number** | `{smn_res['smn1_cn']}` |",
            f"| **SMN2 Copy Number** | `{smn_res['smn2_cn']}` |",
            f"| **SMN2Δ7-8 Copy Number** | `{smn_res['smn2_delta7_8_cn']}` |",
            f"| **isSMA** | `{smn_res['is_sma']}` |",
            f"| **isCarrier** | `{smn_res['is_carrier']}` |",
            f"| **g.27134T>G (2+0 Silent Carrier SNP)** | `{smn_res['g27134tg']}` |",
            f"| **Interpretation** | **{smn_res['interpretation']}** |",
        ])
    else:
        md.append(f"_Unavailable: {smn_res.get('status')} ({smn_res.get('error', '')})_")

    md.extend([
        "",
        "## 2. Full *CYP2D6* Star-Allele, CNV & *CYP2D7* Hybrid Calling (Cyrius)",
    ])
    if cyrius_res.get("status") == "OK":
        md.extend([
            "| Metric | Result |",
            "|---|---|",
            f"| **Caller** | {cyrius_res['tool']} |",
            f"| **Diplotype** | `{cyrius_res['genotype']}` |",
            f"| **Filter** | `{cyrius_res['filter']}` |",
            f"| **CNV Group / Total CYP2D6+2D7 CN / Spacer CN** | `{cyrius_res.get('cnv_group')}` / `{cyrius_res.get('total_cn')}` / `{cyrius_res.get('spacer_cn')}` |",
            f"| **CPIC Phenotype** | **{cyrius_res['cpic_phenotype']}** (`{cyrius_res['cpic_tier']}`) |",
            f"| **CPIC Guidance** | {cyrius_res['cpic_guidance']} |",
        ])
    else:
        md.append(f"_Unavailable: {cyrius_res.get('status')} ({cyrius_res.get('error', '')})_")

    md.extend([
        "",
        "## 3. Repeat Expansions (*FMR1* Fragile X & Neurogenetic Loci — ExpansionHunter)",
    ])
    if eh_res.get("status") == "OK":
        md.append(
            f"Screened **{eh_res['n_loci']}** GRCh38 repeat loci (sex=`{eh_res.get('sex_used')}`); "
            f"**{eh_res['n_flagged']}** locus/loci outside normal range."
        )
        md.extend([
            "",
            "### Key Clinical Repeat Loci",
            "| Locus | Condition | Repeat Unit | Genotype (Units) | 95% CI | Normal Max | Classification | Depth |",
            "|---|---|---|---|---|---:|---|---:|",
        ])
        by_locus = {r["locus"]: r for r in eh_res["loci"]}
        shown: set[str] = set()
        for loc in list(PRIORITY_REPEAT_LOCI) + [r["locus"] for r in eh_res["flagged"]]:
            if loc in shown or loc not in by_locus:
                continue
            shown.add(loc)
            r = by_locus[loc]
            md.append(
                f"| **{r['locus']}** | {r['condition']} | `{r['repeat_unit']}` | `{r['genotype']}` | "
                f"`{r['ci']}` | {r['normal_max'] if r['normal_max'] is not None else '-'} | "
                f"**{r['classification']}** | {r['coverage']}x |"
            )
    else:
        md.append(f"_Unavailable: {eh_res.get('status')} ({eh_res.get('error', '')})_")

    md.extend([
        "",
        "## 4. BAM Read-Level Verification of Coding-Exon CNV/SV Calls",
    ])
    if sv_rows:
        md.extend([
            "| Gene(s) | Call | Length | Inside DP | Flank DP | Depth Ratio | Mean MAPQ | Clipped / Disc. | BAM Verdict |",
            "|---|---|---:|---:|---:|---:|---:|---:|---|",
        ])
        for s in sv_rows:
            md.append(
                f"| **{s['cds_genes']}** | {s['svtype']} `{s['chrom']}:{s['start']:,}-{s['end']:,}` ({s['gt_or_cn']}) | "
                f"{s['length_bp']:,} bp | {s['inside_dp']}x | {s['flank_dp']}x | **{s['depth_ratio']:.2f}** | "
                f"{s['mapq_inside']} | {s['clipped_reads']} / {s['discordant_pairs']} | **{s['bam_verdict']}** |"
            )
    else:
        md.append("_No coding-exon CNV/SV calls to verify._")

    md.extend([
        "",
        "## 5. BAM Read-Backed Phasing (Candidate Compound Heterozygotes)",
    ])
    if phasing_rows:
        md.extend([
            "| Category | Gene | Variant 1 | Variant 2 | Distance | Shared Pairs | Cis / Trans | Verdict |",
            "|---|---|---|---|---:|---:|---:|---|",
        ])
        for p in phasing_rows:
            md.append(
                f"| {p['category']} | **{p['gene']}** | `{p['var1']}` | `{p['var2']}` | "
                f"{p['distance_bp']:,} bp | {p['shared_read_pairs']} | {p['cis_supporting']} / {p['trans_supporting']} | **{p['phase_verdict']}** |"
            )
    else:
        md.append("_No candidate compound heterozygotes identified in Step 4/5._")

    md.extend([
        "",
        "## 6. HLA Typing & CPIC Drug-Hypersensitivity Alleles (T1K)",
    ])
    hla_res = hla_res or {"status": "NOT_RUN"}
    if hla_res.get("status") == "OK":
        md.extend([
            "| Gene | Allele 1 (quality) | Allele 2 (quality) | Called (quality > 0) |",
            "|---|---|---|---|",
        ])
        by_gene = {r["gene"]: r for r in hla_res.get("rows", [])}
        for gene in HLA_REPORT_GENES:
            r = by_gene.get(gene)
            if not r:
                md.append(f"| **{gene}** | - | - | not reported |")
                continue
            md.append(
                f"| **{gene}** | `{r['allele_1']}` ({r['quality_1']}) | `{r['allele_2']}` ({r['quality_2']}) | `{r['called']}` |"
            )
        md.extend(["", "**CPIC risk-allele check (HLA-A / HLA-B):**", ""])
        for e in interpret_hla_typing(hla_res):
            md.append(f"- **[{e['tier']}] `{e['gene']}` — {e['phenotype']}:** {e['guidance']}")
        md.append("")
        md.append(
            f"_Class II genes ({', '.join(g for g in HLA_REPORT_GENES if g not in HLA_PGX_GENES)}) are shown for information only; "
            "HLA typing from short reads is research-grade, so confirm any positive with clinical HLA typing._"
        )
    else:
        md.append(f"_Unavailable: {hla_res.get('status')} ({hla_res.get('error', '')})_")

    md.extend(["", "## 7. HLA Disease Associations (Celiac HLA-DQ2/DQ8, HLA-B27)", ""])
    associations = interpret_hla_associations(hla_res)
    if associations:
        md.extend([
            "| Trait | Typed alleles | Result | What it means | Next step |",
            "|---|---|---|---|---|",
        ])
        for a in associations:
            md.append(
                f"| **{a['trait']}** | `{a['genotype']}` | **[{a['tier']}] {a['result']}** | {a['meaning']} | {a['next_step']} |"
            )
        md.extend([
            "",
            "_Risk associations, not diagnoses: both are most informative when absent. "
            "DQA1 and DQB1 are not phased; a DQ2/DQ8 molecule formed in trans carries the same risk. "
            "Confirm with clinical HLA typing before relying on a result._",
        ])
    else:
        md.append("_Unavailable: needs T1K HLA typing (section 6)._")

    md.extend(["", "## 8. *UGT1A1* Promoter TA Repeat (\\*28, rs3064744)", ""])
    md.extend(render_ugt1a1_section(ugt_res))
    return "\n".join(md) + "\n"


def render_ugt1a1_section(ugt: dict[str, Any] | None) -> list[str]:
    ugt = ugt or {"status": "NOT_RUN"}
    if "counts" not in ugt:
        return [f"_Unavailable: {ugt.get('status')} ({ugt.get('error', '')})_"]
    if ugt["status"] == "OK":
        verdict = f"**{ugt['diplotype']}**"
    elif ugt["status"] == "LOW_DEPTH":
        verdict = f"**No call:** fewer than {ugt['min_spanning']} spanning read pairs"
    else:
        verdict = "**Not resolved:** " + " or ".join(f"`{c}`" for c in ugt.get("candidate_labels") or [])
    lod = "—" if ugt.get("lod") is None else f"{ugt['lod']}"
    return [
        "| Metric | Result |",
        "|---|---|",
        f"| **Repeat lengths (read pairs spanning both anchors)** | {counts_label(ugt['counts']) or 'none'} |",
        f"| **Spanning read pairs** | {ugt['spanning_fragments']} (minimum {ugt['min_spanning']}) |",
        f"| **Reads lacking an anchor / non-(TA)n / mates disagreeing** | {ugt['partial_reads']} / "
        f"{ugt['noncanonical_fragments']} / {ugt['mate_conflict_fragments']} |",
        f"| **Diplotype** | {verdict} |",
        f"| **Confidence (log10 likelihood ratio, best vs next)** | {lod} (a call needs >= {ugt['min_lod']}) |",
        "",
        "_TA5 = \\*36, TA6 = \\*1 (reference), TA7 = \\*28, TA8 = \\*37. \\*28 and \\*37 lower UGT1A1 activity "
        "(irinotecan toxicity, atazanavir jaundice, Gilbert syndrome). Research-grade read count: a 30x genome gives "
        "only about 10-25 spanning read pairs, so an unresolved call needs a clinical UGT1A1 test. The medication card "
        "combines this call with the VCF's \\*6 call._",
    ]


def refresh_pgx_card_with_bam(
    sample_label: str,
    rep_dir: Path,
    cyrius_res: dict[str, Any],
    hla_res: dict[str, Any],
    sex: str,
    ugt_res: dict[str, Any] | None = None,
) -> None:
    """Re-render pgx_medication_card.{md,tsv} with BAM-resolved Cyrius CYP2D6, T1K HLA and UGT1A1 TA-repeat calls."""
    tm_tsv = rep_dir / "targeted_markers.tsv"
    bam_calls = cyrius_res.get("status") == "OK" or hla_res.get("status") == "OK" or "counts" in (ugt_res or {})
    if not tm_tsv.exists() or not bam_calls:
        return
    df = pd.read_csv(tm_tsv, sep="\t", dtype={"gt": str})
    marker_rows = df.to_dict("records")
    pgx_entries = interpret_pgx_profile(marker_rows, cyrius_cyp2d6=cyrius_res, sex=sex, hla=hla_res, ugt1a1_ta=ugt_res)
    pd.DataFrame(pgx_entries).to_csv(rep_dir / "pgx_medication_card.tsv", sep="\t", index=False)
    ref_meta = get_reference_metadata(REPO_ROOT / "data" / "reference")
    card_md = render_pgx_medication_card(sample_label, ref_meta, pgx_entries)
    (rep_dir / "pgx_medication_card.md").write_text(card_md, encoding="utf-8")
    print(f"[{sample_label}] Updated PGx Medication Card with BAM calls (Cyrius {cyrius_res.get('status')}, "
          f"T1K {hla_res.get('status')}, UGT1A1 TA repeat {(ugt_res or {}).get('status')}).")


def reset_caller_cache_if_alignment_changed(work_dir: Path, bam_path: Path, sample_label: str) -> None:
    """The callers reuse their outputs in work_dir; drop them when the alignment link now points at another file."""
    stamp = work_dir / "alignment.json"
    current = {"alignment_target": os.path.realpath(bam_path), "size": bam_path.stat().st_size}
    if work_dir.exists() and stamp.exists():
        try:
            previous = json.loads(stamp.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = {}
        if previous != current:
            print(f"[{sample_label}] Alignment changed since the cached caller outputs were made; re-running all callers.")
            shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    stamp.write_text(json.dumps(current) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run BAM-based blind-spot callers and read-level verification")
    parser.add_argument("--sample", required=True, help="Sample label from config/samples.tsv")
    parser.add_argument("--threads", type=int, default=4, help="Thread count for callers")
    parser.add_argument(
        "--reference",
        default=str(REPO_ROOT / "data" / "reference" / "GRCh38.primary_assembly.fa.bgz"),
        help="Path to GRCh38 reference FASTA (.fa or .fa.bgz)",
    )
    args = parser.parse_args()

    sample_label = args.sample
    bam_path = find_alignment_file(sample_label)
    if not bam_path:
        print(f"[{sample_label}] No BAM/CRAM found in data/alignments/{sample_label}/ — skipping Step 10.")
        return
    if not bam_path.exists():
        # Broken symlink (e.g. external drive holding the BAM is not mounted): keep the previous
        # bam_blind_spots.json/report untouched rather than failing the whole annual run.
        print(
            f"[{sample_label}] Alignment link {bam_path.name} points to a missing file "
            f"({os.path.realpath(bam_path)}); is its drive mounted? Skipping Step 10 and keeping previous BAM results."
        )
        return

    rep_dir = REPO_ROOT / "results" / "reports" / sample_label
    work_dir = REPO_ROOT / "data" / "processed" / sample_label / "bam_callers"
    verified, why = bam_verified(rep_dir, bam_path)
    if not verified:
        print(
            f"[{sample_label}] ERROR: Step 10 refuses to run: {why}. Run "
            f"'.venv-wsl/bin/python scripts/03_sample_identity_qc.py --sample {sample_label}' and fix any mismatch first.",
            file=sys.stderr,
        )
        sys.exit(1)
    identity_qc = load_identity_qc(rep_dir)
    rep_dir.mkdir(parents=True, exist_ok=True)
    reset_caller_cache_if_alignment_changed(work_dir, bam_path, sample_label)

    sex = get_sample_sex(sample_label)
    has_chr = detect_bam_has_chr(bam_path)
    fasta_path = Path(args.reference)

    eh_res = run_expansion_hunter(bam_path, fasta_path, sex, work_dir, sample_label, args.threads, has_chr)
    smn_res = run_smn_caller(bam_path, work_dir, sample_label, args.threads, fasta_path)
    cyrius_res = run_cyrius(bam_path, work_dir, sample_label, args.threads, fasta_path)
    hla_res = run_t1k_hla(bam_path, work_dir, sample_label, args.threads, has_chr, fasta_path)
    print(f"[{sample_label}] Counting the UGT1A1 promoter TA repeat from spanning reads...")
    ugt_res = analyze_ta_repeat(bam_path, has_chr, reference=fasta_path)

    print(f"[{sample_label}] Verifying coding-exon SVs in BAM...")
    sv_rows = verify_exonic_svs_in_bam(bam_path, rep_dir / "cnv_sv_exonic_clinical.tsv", has_chr, fasta_path)
    print(f"[{sample_label}] Checking read-backed phasing for candidate compound hets...")
    phasing_rows = check_candidate_phasing(bam_path, rep_dir, has_chr, fasta_path)

    payload = {
        "sample": sample_label,
        "run_date": datetime.date.today().isoformat(),
        "alignment_file": bam_path.name,
        "sex": sex,
        "has_chr_contigs": has_chr,
        "expansion_hunter": eh_res,
        "smn": smn_res,
        "cyrius_cyp2d6": cyrius_res,
        "hla_t1k": hla_res,
        "hla_associations": interpret_hla_associations(hla_res),
        "ugt1a1_ta_repeat": ugt_res,
        "identity_qc": {k: (identity_qc or {}).get(k) for k in ("status", "summary", "run_date")},
        "sv_verification": sv_rows,
        "phasing_checks": phasing_rows,
    }
    out_json = rep_dir / "bam_blind_spots.json"
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    if eh_res.get("status") == "OK" and eh_res.get("loci"):
        pd.DataFrame(eh_res["loci"]).to_csv(rep_dir / "bam_repeat_expansions.tsv", sep="\t", index=False)
    if sv_rows:
        pd.DataFrame(sv_rows).to_csv(rep_dir / "bam_sv_verification.tsv", sep="\t", index=False)
    if payload["hla_associations"]:
        pd.DataFrame(payload["hla_associations"]).to_csv(rep_dir / "hla_disease_associations.tsv", sep="\t", index=False)

    report_md = render_bam_report(
        sample_label, bam_path, eh_res, smn_res, cyrius_res, sv_rows, phasing_rows, hla_res, ugt_res, identity_qc
    )
    out_md = rep_dir / "bam_blind_spots_report.md"
    out_md.write_text(report_md, encoding="utf-8")
    print(f"[{sample_label}] BAM blind-spot report written to {out_md}")

    refresh_pgx_card_with_bam(sample_label, rep_dir, cyrius_res, hla_res, sex, ugt_res)


if __name__ == "__main__":
    main()
