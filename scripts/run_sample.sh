#!/usr/bin/env bash
set -euo pipefail

# Reproducible end-to-end single-sample runner using config/samples.tsv
# Usage: bash scripts/run_sample.sh <sample_label> [--force-norm]

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <sample_label> [--force-norm]" >&2
  exit 1
fi

LABEL="$1"
FORCE_NORM=0
if [[ "${2:-}" == "--force-norm" ]]; then
  FORCE_NORM=1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

SAMPLES_TSV="config/samples.tsv"
if [[ ! -f "${SAMPLES_TSV}" ]]; then
  echo "ERROR: ${SAMPLES_TSV} not found. Copy config/samples.tsv.example first." >&2
  exit 1
fi

IFS=$'\t' read -r _lbl SAMPLE_ID SEX RAW_DIR SNP_VCF CNV_VCF SV_VCF < <(
  awk -F'\t' -v target="${LABEL}" 'NR>1 && $1==target {print; exit}' "${SAMPLES_TSV}"
)

if [[ -z "${SAMPLE_ID:-}" ]]; then
  echo "ERROR: Sample '${LABEL}' not found in ${SAMPLES_TSV}" >&2
  exit 1
fi

PROCESSED_DIR="data/processed/${LABEL}"
QC_DIR="results/qc/${LABEL}"
REFERENCE_DIR="data/reference"
FASTA_BGZ="${REFERENCE_DIR}/GRCh38.primary_assembly.fa.bgz"
CHR_MAP="${REFERENCE_DIR}/chr_to_nochr.map"
PRIMARY_CHROMS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT"
THREADS="${THREADS:-4}"

mkdir -p "${PROCESSED_DIR}" "${QC_DIR}" "results/reports/${LABEL}" "data/alignments/${LABEL}"

echo "============================================================"
echo "Running pipeline for sample: ${LABEL} (${SAMPLE_ID}, sex=${SEX})"
echo "Read-only source dir: ${RAW_DIR}"

index_raw_vcf() {
  local tag="$1"
  local fname="$2"
  local raw_path="${RAW_DIR}/${fname}"
  local out_gz="${PROCESSED_DIR}/${tag}.vcf.gz"
  local header_file="${QC_DIR}/${tag}.header.txt"

  if [[ -f "${out_gz}" && -f "${out_gz}.tbi" ]]; then
    echo "[${LABEL}:${tag}] Already indexed: ${out_gz}"
  else
    if [[ ! -r "${raw_path}" ]]; then
      echo "ERROR: Cannot read ${raw_path}" >&2
      exit 1
    fi
    if [[ ! -f "${header_file}" ]]; then
      awk '/^#/ {print; next} {exit}' "${raw_path}" > "${header_file}"
    fi
    local has_chr
    has_chr="$(awk '/^##contig=<ID=chr/ {print "yes"; exit} /^[^#]/ {if ($1 ~ /^chr/) print "yes"; exit}' "${raw_path}")"
    if [[ "${has_chr}" == "yes" ]]; then
      echo "[${LABEL}:${tag}] Harmonizing 'chr' contigs + bgzip (-@ ${THREADS}) + tabix -> ${out_gz}"
      bcftools annotate --threads "${THREADS}" --rename-chrs "${CHR_MAP}" -Oz -o "${out_gz}" "${raw_path}"
    else
      echo "[${LABEL}:${tag}] bgzip (-@ ${THREADS}) + tabix -> ${out_gz}"
      bgzip -@ "${THREADS}" -c "${raw_path}" > "${out_gz}"
    fi
    tabix -f -p vcf "${out_gz}"
  fi
  bcftools index -s "${out_gz}" > "${QC_DIR}/${tag}.contigs.tsv"
}

index_raw_vcf "snp_indel.genome" "${SNP_VCF}"
index_raw_vcf "cnv" "${CNV_VCF}"
index_raw_vcf "sv" "${SV_VCF}"

VAR_NORM="${PROCESSED_DIR}/snp_indel.variants.norm.vcf.gz"
GVCF_IN="${PROCESSED_DIR}/snp_indel.genome.vcf.gz"
if [[ "${FORCE_NORM}" -eq 1 || ! -f "${VAR_NORM}" || ! -f "${VAR_NORM}.tbi" ]]; then
  echo "[${LABEL}] Extracting variant-only VCF & running bcftools norm -> ${VAR_NORM}"
  HAS_CHR="$(bcftools index -s "${GVCF_IN}" | awk '$1 ~ /^chr/ {print "yes"; exit}')"
  if [[ "${HAS_CHR}" == "yes" ]]; then
    echo "[${LABEL}] Detected 'chr' contig prefix; harmonizing to GRCh38 (1..22, X, Y, MT)..."
    bcftools annotate --threads "${THREADS}" --rename-chrs "${CHR_MAP}" -Ou "${GVCF_IN}" \
      | bcftools view --threads "${THREADS}" -e 'ALT="."' -Ou \
      | bcftools norm --threads "${THREADS}" -m -any -f "${FASTA_BGZ}" -c w -Oz -o "${VAR_NORM}"
  else
    bcftools view --threads "${THREADS}" -r "${PRIMARY_CHROMS}" -e 'ALT="."' -Ou "${GVCF_IN}" \
      | bcftools norm --threads "${THREADS}" -m -any -f "${FASTA_BGZ}" -c w -Oz -o "${VAR_NORM}"
  fi
  tabix -f -p vcf "${VAR_NORM}"
  echo "[${LABEL}] Normalized variant records: $(bcftools index -n "${VAR_NORM}")"
else
  echo "[${LABEL}] Already normalized: ${VAR_NORM} ($(bcftools index -n "${VAR_NORM}") records)"
fi

echo "[${LABEL}] Running Step 4 (ClinVar + GenCC + ACMG SF v3.3 screening)..."
.venv-wsl/bin/python scripts/04_carrier_screen.py --sample "${LABEL}"

echo "[${LABEL}] Running Step 5 (cached gnomAD v4 + Ensembl + novel ACMG LoF rescue)..."
.venv-wsl/bin/python scripts/05_annotate_candidates.py --sample "${LABEL}"

echo "[${LABEL}] Running Step 7 (exonic CNV/SV, strict PASS gVCF coverage, PGx & CPIC Medication Card)..."
.venv-wsl/bin/python scripts/07_deep_analysis.py --sample "${LABEL}"

if [[ -f "data/reference/prs/ref_distributions.json" ]]; then
  echo "[${LABEL}] Running Step 9 (1000G-calibrated Polygenic Risk Scores)..."
  .venv-wsl/bin/python scripts/09_polygenic_risk.py --sample "${LABEL}"
else
  echo "[${LABEL}] PRS reference distributions not built yet (run 'bash scripts/setup_prs_references.sh' to enable Step 9)."
fi

ALIGN_FOUND="$(find "data/alignments/${LABEL}" -maxdepth 1 \( -type f -o -type l \) \( -name '*.bam' -o -name '*.cram' \) 2>/dev/null | head -n 1 || true)"
if [[ -n "${ALIGN_FOUND}" ]]; then
  echo "[${LABEL}] Alignment file detected: ${ALIGN_FOUND}"
  if [[ ! -f "${QC_DIR}/bam.idxstats.tsv" ]]; then
    echo "[${LABEL}] Running samtools idxstats -> ${QC_DIR}/bam.idxstats.tsv"
    samtools idxstats "${ALIGN_FOUND}" > "${QC_DIR}/bam.idxstats.tsv"
  fi
  echo "[${LABEL}] Running Step 10 (BAM blind-spot callers: ExpansionHunter, SMNCopyNumberCaller, Cyrius, SV/phasing check)..."
  .venv-wsl/bin/python scripts/10_bam_blind_spots.py --sample "${LABEL}" --threads "${THREADS}"
else
  echo "[${LABEL}] No BAM/CRAM alignment file found in data/alignments/${LABEL}/ (skipping alignment-based blind-spot callers)."
fi

echo "============================================================"
echo "Pipeline complete for ${LABEL}!"
