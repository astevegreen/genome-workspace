"""Sample identity QC: does the BAM come from the same person as the VCF, and does genetic sex match samples.tsv?

BAM vs VCF: about 3,000 well-called SNVs (one per ~1 Mb of autosome, heterozygous or homozygous ALT in
the VCF) are genotyped from the BAM with `bcftools mpileup | bcftools call` and compared allele by allele.
The comparison is counted here rather than with `bcftools gtcheck`: gtcheck matches records on POS+ALT,
so it drops every site where the BAM call is homozygous reference (ALT '.'), and those are exactly the
sites where another person's BAM disagrees.

Sex: chrX and chrY read depth relative to the autosomes (from the BAM index), and the chrX non-PAR
heterozygous / haploid call pattern of the VCF, are each compared with the sex in config/samples.tsv.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Any

AUTOSOMES = tuple(str(i) for i in range(1, 23))
GRCH38_X_PARS = ((10_001, 2_781_479), (155_701_383, 156_030_895))

SITE_SPACING = 1_000_000
SITE_WINDOW = 200_000
SITE_MIN_DP = 15
SITE_MIN_GQ = 30
BAM_MIN_DP = 8
MIN_SITES = 500
MATCH_MAX_DISCORDANCE = 0.05  # same person: typically < 1%
MISMATCH_MIN_DISCORDANCE = 0.20  # unrelated people: ~30-50% at one person's variant sites

# Relative depth (per bp, versus the autosomes) by sex: XY ~0.5 / ~0.4, XX ~1.0 / <0.05.
MALE_MAX_X_RATIO = 0.70
FEMALE_MIN_X_RATIO = 0.80
MALE_MIN_Y_RATIO = 0.15
FEMALE_MAX_Y_RATIO = 0.10
# VCF chrX non-PAR: share of variant calls that are heterozygous (XX ~0.5+, XY ~0) or haploid (XY callers)
FEMALE_MIN_X_HET = 0.25
MALE_MAX_X_HET = 0.15
MALE_MIN_X_HAPLOID = 0.90

QC_JSON = "sample_identity_qc.json"
QC_MD = "sample_identity_qc.md"


# --------------------------------------------------------------------------- BAM vs VCF genotypes

def _gt_alleles(rec: Any) -> tuple[str, ...] | None:
    s = rec.samples[0]
    gt = s.get("GT")
    if not gt or any(a is None for a in gt):
        return None
    alleles = rec.alleles or ()
    if any(a >= len(alleles) for a in gt):
        return None
    return tuple(sorted(alleles[a] for a in gt))


def _site_ok(rec: Any) -> bool:
    if any(f != "PASS" for f in rec.filter.keys()):
        return False
    if len(rec.ref or "") != 1 or not rec.alts or len(rec.alts) != 1 or rec.alts[0] not in ("A", "C", "G", "T"):
        return False
    s = rec.samples[0]
    gt = s.get("GT")
    if not gt or len(gt) != 2 or any(a is None for a in gt) or not any(a > 0 for a in gt):
        return False
    try:
        dp, gq = s.get("DP"), s.get("GQ")
    except KeyError:
        return False
    return dp is not None and dp >= SITE_MIN_DP and (gq is None or gq >= SITE_MIN_GQ)


def select_identity_sites(vcf_path: Path, spacing: int = SITE_SPACING, window: int = SITE_WINDOW) -> list[dict[str, Any]]:
    """One well-called biallelic SNV (het or hom ALT) per `spacing` bp of each autosome."""
    import pysam

    sites: list[dict[str, Any]] = []
    with pysam.VariantFile(str(vcf_path)) as vf:
        for chrom in AUTOSOMES:
            contig = chrom if chrom in vf.header.contigs else f"chr{chrom}"
            if contig not in vf.header.contigs:
                continue
            length = vf.header.contigs[contig].length or 0
            for start in range(spacing // 2, max(0, length - window), spacing):
                recs = list(vf.fetch(contig, start, start + window))
                by_pos: dict[int, list[Any]] = {}
                for r in recs:
                    by_pos.setdefault(r.pos, []).append(r)
                for pos in sorted(by_pos):
                    if len(by_pos[pos]) == 1 and _site_ok(by_pos[pos][0]):
                        r = by_pos[pos][0]
                        sites.append({"chrom": contig.removeprefix("chr"), "pos": pos, "ref": r.ref,
                                      "alt": r.alts[0], "gt": _gt_alleles(r)})
                        break
    return sites


def genotype_bam_at_sites(
    bam_path: Path, fasta: Path, sites: list[dict[str, Any]], work_dir: Path, has_chr: bool, bcftools: str = "bcftools"
) -> dict[tuple[str, int], dict[str, Any]]:
    """Genotype the BAM at the sites with bcftools mpileup/call; {(chrom, pos): {'gt': alleles, 'dp': n}}."""
    import pysam

    work_dir.mkdir(parents=True, exist_ok=True)
    regions = work_dir / "identity_sites.tsv"
    with open(regions, "w", encoding="utf-8", newline="\n") as fh:
        for s in sites:
            fh.write(f"{'chr' if has_chr else ''}{s['chrom']}\t{s['pos']}\n")
    out_vcf = work_dir / "identity_bam_calls.vcf"
    cmd = (
        f"{bcftools} mpileup -f '{fasta}' -R '{regions}' -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 -Ou '{bam_path}' "
        f"| {bcftools} call -m -Ov -o '{out_vcf}'"
    )
    proc = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"bcftools mpileup/call failed: {proc.stderr.strip()[-800:]}")
    calls: dict[tuple[str, int], dict[str, Any]] = {}
    with pysam.VariantFile(str(out_vcf)) as vf:
        for r in vf:
            s = r.samples[0]
            calls[(r.chrom.removeprefix("chr"), r.pos)] = {"gt": _gt_alleles(r), "dp": s.get("DP")}
    return calls


def compare_genotypes(
    sites: list[dict[str, Any]], bam_calls: dict[tuple[str, int], dict[str, Any]], min_dp: int = BAM_MIN_DP
) -> dict[str, Any]:
    """Allele-level comparison of VCF and BAM genotypes at sites the BAM covers at >= min_dp."""
    compared = matching = low_depth = 0
    for s in sites:
        c = bam_calls.get((s["chrom"], s["pos"]))
        if not c or c.get("gt") is None or (c.get("dp") or 0) < min_dp:
            low_depth += 1
            continue
        compared += 1
        if tuple(sorted(c["gt"])) == tuple(sorted(s["gt"])):
            matching += 1
    disc = round(1 - matching / compared, 4) if compared else None
    return {"sites": len(sites), "compared": compared, "matching": matching, "low_depth": low_depth, "discordance": disc}


def classify_identity(discordance: float | None, compared: int) -> str:
    if discordance is None or compared < MIN_SITES:
        return "INCONCLUSIVE"
    if discordance <= MATCH_MAX_DISCORDANCE:
        return "MATCH"
    if discordance >= MISMATCH_MIN_DISCORDANCE:
        return "MISMATCH"
    return "INCONCLUSIVE"


# --------------------------------------------------------------------------- sex

def coverage_ratios(index_stats: list[tuple[str, int, int]]) -> dict[str, float | None]:
    """chrX / chrY reads per bp relative to the autosomes, from (contig, length, mapped) rows."""
    auto_len = auto_reads = 0
    x = y = None
    for contig, length, mapped in index_stats:
        c = contig.removeprefix("chr")
        if c in AUTOSOMES:
            auto_len += length
            auto_reads += mapped
        elif c == "X":
            x = (length, mapped)
        elif c == "Y":
            y = (length, mapped)
    if not auto_len or not auto_reads:
        return {"x_ratio": None, "y_ratio": None}
    base = auto_reads / auto_len
    return {
        "x_ratio": round((x[1] / x[0]) / base, 3) if x and x[0] else None,
        "y_ratio": round((y[1] / y[0]) / base, 3) if y and y[0] else 0.0,
    }


def sex_from_coverage(x_ratio: float | None, y_ratio: float | None) -> str:
    if x_ratio is None or y_ratio is None:
        return "UNCLEAR"
    if x_ratio <= MALE_MAX_X_RATIO and y_ratio >= MALE_MIN_Y_RATIO:
        return "M"
    if x_ratio >= FEMALE_MIN_X_RATIO and y_ratio <= FEMALE_MAX_Y_RATIO:
        return "F"
    return "UNCLEAR"


def vcf_chrx_metrics(vcf_path: Path) -> dict[str, Any]:
    """Counts of chrX non-PAR variant calls: heterozygous, homozygous/hemizygous ALT, and haploid; chrY ALT calls."""
    import pysam

    het = hom = haploid = y_alt = 0
    with pysam.VariantFile(str(vcf_path)) as vf:
        contigs = set(vf.header.contigs)
        xc = "X" if "X" in contigs else ("chrX" if "chrX" in contigs else None)
        yc = "Y" if "Y" in contigs else ("chrY" if "chrY" in contigs else None)
        if xc:
            for rec in vf.fetch(xc):
                if any(a <= rec.pos <= b for a, b in GRCH38_X_PARS):
                    continue
                gt = rec.samples[0].get("GT")
                if not gt or any(a is None for a in gt):
                    continue
                if len(gt) == 1:
                    haploid += gt[0] > 0
                elif gt[0] != gt[1]:
                    het += 1
                elif gt[0] > 0:
                    hom += 1
        if yc:
            for rec in vf.fetch(yc):
                gt = rec.samples[0].get("GT")
                if gt and any(a is not None and a > 0 for a in gt):
                    y_alt += 1
    called = het + hom + haploid
    return {
        "x_het": het, "x_hom_alt": hom, "x_haploid": haploid,
        "x_het_ratio": round(het / called, 4) if called else None,
        "x_haploid_ratio": round(haploid / called, 4) if called else None,
        "y_alt_calls": y_alt,
    }


def sex_from_vcf(m: dict[str, Any]) -> str:
    if m.get("x_het_ratio") is None:
        return "UNCLEAR"
    if (m.get("x_haploid_ratio") or 0) >= MALE_MIN_X_HAPLOID:
        return "M"
    if m["x_het_ratio"] >= FEMALE_MIN_X_HET:
        return "F"
    if m["x_het_ratio"] <= MALE_MAX_X_HET:
        return "M"
    return "UNCLEAR"


def sex_status(declared: str, inferred: str | None) -> str:
    if inferred is None:
        return "NOT_CHECKED"
    if inferred == "UNCLEAR" or declared not in ("M", "F"):
        return "UNCLEAR"
    return "MATCH" if inferred == declared else "MISMATCH"


# --------------------------------------------------------------------------- verdict and reuse

def overall_status(bam_status: str, sex_statuses: list[str]) -> str:
    """FAIL stops the BAM steps; WARN is reported but does not."""
    if bam_status in ("MISMATCH", "INCONCLUSIVE") or "MISMATCH" in sex_statuses:
        return "FAIL"
    if "UNCLEAR" in sex_statuses:
        return "WARN"
    return "PASS"


def alignment_fingerprint(bam_path: Path) -> dict[str, Any]:
    """What the BAM check was run on: the link target, size and mtime (a swapped link changes all three)."""
    st = bam_path.stat()
    return {"alignment": bam_path.name, "alignment_target": os.path.realpath(bam_path),
            "size": st.st_size, "mtime": int(st.st_mtime)}


def load_identity_qc(rep_dir: Path) -> dict[str, Any] | None:
    p = Path(rep_dir) / QC_JSON
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def bam_results_usable(rep_dir: Path) -> bool:
    """False when the identity QC failed, so BAM-derived results must not reach the reports."""
    qc = load_identity_qc(rep_dir)
    if qc is None:
        return True  # results from before this check existed
    return qc.get("status") != "FAIL"


def load_bam_results(rep_dir: Path) -> dict[str, Any] | None:
    """Step 10's bam_blind_spots.json, or None if missing or the identity QC failed."""
    p = Path(rep_dir) / "bam_blind_spots.json"
    if not p.exists() or not bam_results_usable(rep_dir):
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def bam_verified(rep_dir: Path, bam_path: Path) -> tuple[bool, str]:
    """Step 10 guard: the identity QC passed and its BAM check was run on this exact alignment."""
    qc = load_identity_qc(rep_dir)
    if qc is None:
        return False, "no sample identity QC result"
    if qc.get("status") == "FAIL":
        return False, f"sample identity QC failed ({qc.get('summary', '')})"
    bam = qc.get("bam_check") or {}
    if bam.get("status") != "MATCH":
        return False, f"BAM identity status is {bam.get('status', 'NOT_CHECKED')}"
    if bam.get("alignment_target") != os.path.realpath(bam_path):
        return False, "the identity check was run on a different alignment file"
    return True, "ok"


def identity_line(qc: dict[str, Any] | None) -> str:
    """One-line summary for other reports."""
    if not qc:
        return "Sample identity QC: not run."
    return f"Sample identity QC: **{qc.get('status')}** — {qc.get('summary', '')}"


def _pct(x: float | None) -> str:
    return "—" if x is None else f"{100 * x:.1f}%"


def render_identity_report(label: str, qc: dict[str, Any]) -> str:
    bam = qc.get("bam_check") or {}
    sex = qc.get("sex_check") or {}
    md = [
        f"# Sample Identity QC — {label.title()}",
        "",
        f"**Run date:** {qc.get('run_date', '')} | **Result:** **{qc.get('status')}** — {qc.get('summary', '')}",
        "",
        "Checks that the BAM and the VCF come from the same person, and that genetic sex matches "
        "`config/samples.tsv`. A failure stops the BAM-based steps (Step 10) and keeps their results off the reports.",
        "",
        "| Check | Result | Detail |",
        "|---|---|---|",
    ]
    if bam.get("status") == "NOT_CHECKED":
        md.append(f"| BAM vs VCF genotypes | NOT CHECKED | {bam.get('reason', 'no BAM')} |")
    else:
        carried = " (carried over: the BAM was not readable this run)" if bam.get("carried_over") else ""
        md.append(
            f"| BAM vs VCF genotypes | **{bam.get('status')}** | {bam.get('matching', 0):,} of {bam.get('compared', 0):,} "
            f"sites agree ({_pct(bam.get('discordance'))} discordant; same person <= {_pct(MATCH_MAX_DISCORDANCE)}, "
            f"different person >= {_pct(MISMATCH_MIN_DISCORDANCE)}){carried} |"
        )
        for o in bam.get("other_samples") or []:
            md.append(
                f"| BAM vs {o['label']}'s VCF | **{o['status']}** | {_pct(o.get('discordance'))} discordant over "
                f"{o.get('compared', 0):,} sites |"
            )
        if bam.get("read_group_samples"):
            md.append(
                f"| BAM read-group sample | {'matches' if bam.get('read_group_matches') else 'differs'} | "
                f"`{', '.join(bam['read_group_samples'])}` vs sample_id `{qc.get('sample_id', '')}` (information only) |"
            )
    cov = sex.get("bam") or {}
    if cov:
        md.append(
            f"| Genetic sex, BAM depth | **{cov.get('inferred')}** ({cov.get('status')}) | chrX {cov.get('x_ratio')}x, "
            f"chrY {cov.get('y_ratio')}x of autosomal depth (XY ~0.5 / ~0.4; XX ~1.0 / <0.1) |"
        )
    v = sex.get("vcf") or {}
    if v:
        md.append(
            f"| Genetic sex, VCF chrX | **{v.get('inferred')}** ({v.get('status')}) | non-PAR chrX calls: "
            f"{v.get('x_het', 0):,} het, {v.get('x_hom_alt', 0):,} hom ALT, {v.get('x_haploid', 0):,} haploid "
            f"(het share {_pct(v.get('x_het_ratio'))}) |"
        )
    md.append(f"| Declared sex (`samples.tsv`) | {sex.get('declared') or '—'} | |")
    md.extend([
        "",
        f"_Method: {bam.get('sites', 0):,} autosomal SNVs called in the VCF (PASS, DP >= {SITE_MIN_DP}, GQ >= {SITE_MIN_GQ}, "
        f"one per {SITE_SPACING // 1_000_000} Mb) re-genotyped from the BAM with bcftools mpileup/call and compared where the "
        f"BAM has >= {BAM_MIN_DP} reads._",
    ])
    return "\n".join(md) + "\n"
