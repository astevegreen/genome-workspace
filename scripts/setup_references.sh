#!/usr/bin/env bash
set -euo pipefail

# Unified Reference Setup & Annual Refresh Script
# Downloads/indexes static GRCh38 references (FASTA, GENCODE v46 CDS BED) once,
# and downloads/normalizes dynamic clinical databases (NCBI ClinVar GRCh38, GenCC).
#
# Usage:
#   bash scripts/setup_references.sh            # Download only missing files
#   bash scripts/setup_references.sh --refresh  # Annual refresh of ClinVar & GenCC

REFRESH=0
for arg in "$@"; do
  case "${arg}" in
    --refresh) REFRESH=1 ;;
    *) echo "Unknown option: ${arg} (use --refresh)" >&2; exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

REFERENCE_DIR="data/reference"
mkdir -p "${REFERENCE_DIR}"
THREADS="${THREADS:-4}"
PRIMARY_CHROMS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT"

# Optional custom clinical VCF from config/paths.env if present
if [[ -f "config/paths.env" ]]; then
  # shellcheck disable=SC1091
  source "config/paths.env"
fi

# 1. Build reusable contig harmonization map (chr1 -> 1, chrM -> MT)
CHR_MAP="${REFERENCE_DIR}/chr_to_nochr.map"
cat > "${CHR_MAP}" <<'EOF'
chr1 1
chr2 2
chr3 3
chr4 4
chr5 5
chr6 6
chr7 7
chr8 8
chr9 9
chr10 10
chr11 11
chr12 12
chr13 13
chr14 14
chr15 15
chr16 16
chr17 17
chr18 18
chr19 19
chr20 20
chr21 21
chr22 22
chrX X
chrY Y
chrM MT
chrMT MT
EOF

# 2. Static: Ensembl GRCh38 primary assembly FASTA
FASTA_URL="https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
FASTA_BGZ="${REFERENCE_DIR}/GRCh38.primary_assembly.fa.bgz"
if [[ ! -f "${FASTA_BGZ}" || ! -f "${FASTA_BGZ}.fai" ]]; then
  RAW_FA_GZ="${REFERENCE_DIR}/GRCh38.primary_assembly.raw.fa.gz"
  echo "[1/4] Downloading GRCh38 primary assembly FASTA (~900 MB)..."
  curl -fL --retry 3 -o "${RAW_FA_GZ}" "${FASTA_URL}"
  gzip -dc "${RAW_FA_GZ}" | bgzip -@ "${THREADS}" > "${FASTA_BGZ}"
  samtools faidx "${FASTA_BGZ}"
  rm -f "${RAW_FA_GZ}"
else
  echo "[1/4] GRCh38 FASTA ready: ${FASTA_BGZ}"
fi

# 3. Static: GENCODE v46 CDS BED
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz"
GTF_GZ="${REFERENCE_DIR}/gencode.v46.basic.annotation.gtf.gz"
CDS_BED_GZ="${REFERENCE_DIR}/gencode.v46.cds.nochr.bed.gz"
if [[ ! -f "${GTF_GZ}" ]]; then
  echo "[2/4] Downloading GENCODE v46 basic annotation GTF (~30 MB)..."
  curl -fL --retry 3 -o "${GTF_GZ}" "${GTF_URL}"
fi
if [[ ! -f "${CDS_BED_GZ}" || ! -f "${CDS_BED_GZ}.tbi" ]]; then
  echo "[2/4] Extracting harmonized GENCODE v46 CDS BED..."
  gzip -dc "${GTF_GZ}" \
    | awk -F'\t' '
      $0 !~ /^#/ && $3 == "CDS" {
        chr = $1
        sub(/^chr/, "", chr)
        if (chr == "M") chr = "MT"
        if (chr ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/) {
          gene = ""
          if (match($9, /gene_name "([^"]+)"/, m)) gene = m[1]
          if (gene != "") {
            start = $4 - 1
            if (start < 0) start = 0
            print chr "\t" start "\t" $5 "\t" gene
          }
        }
      }
    ' \
    | sort -k1,1 -k2,2n -k3,3n -u \
    | bgzip -c > "${CDS_BED_GZ}"
  tabix -f -p bed "${CDS_BED_GZ}"
else
  echo "[2/4] GENCODE v46 CDS BED ready: ${CDS_BED_GZ}"
fi

# 4. Dynamic: NCBI ClinVar GRCh38 (re-downloaded on --refresh)
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_RAW="${REFERENCE_DIR}/clinvar_GRCh38.vcf.gz"
CLINVAR_NORM="${REFERENCE_DIR}/clinvar_GRCh38.norm.vcf.gz"
CLINVAR_CLINICAL="${REFERENCE_DIR}/clinvar_clinical.norm.vcf.gz"
if [[ "${REFRESH}" -eq 1 || ! -f "${CLINVAR_RAW}" || ! -f "${CLINVAR_RAW}.tbi" ]]; then
  echo "[3/4] Fetching latest NCBI ClinVar GRCh38 VCF..."
  curl -fL --retry 3 -o "${CLINVAR_RAW}" "${CLINVAR_URL}"
  curl -fL --retry 3 -o "${CLINVAR_RAW}.tbi" "${CLINVAR_URL}.tbi"
  curl -fL --retry 3 -o "${CLINVAR_RAW}.md5" "${CLINVAR_URL}.md5"
  EXPECTED_MD5="$(awk '{print $1}' "${CLINVAR_RAW}.md5")"
  ACTUAL_MD5="$(md5sum "${CLINVAR_RAW}" | awk '{print $1}')"
  if [[ "${EXPECTED_MD5}" != "${ACTUAL_MD5}" ]]; then
    echo "ERROR: ClinVar MD5 mismatch! Expected ${EXPECTED_MD5}, got ${ACTUAL_MD5}" >&2
    exit 1
  fi
  echo "[3/4] ClinVar MD5 verified: ${ACTUAL_MD5}"
  rm -f "${CLINVAR_NORM}" "${CLINVAR_NORM}.tbi" "${CLINVAR_CLINICAL}" "${CLINVAR_CLINICAL}.tbi"
  if [[ "${REFRESH}" -eq 1 && -f "${REFERENCE_DIR}/api_cache.json" ]]; then
    echo "[3/4] Refreshing local gnomAD/Ensembl API cache for annual checkup..."
    rm -f "${REFERENCE_DIR}/api_cache.json"
  fi
else
  echo "[3/4] ClinVar raw VCF ready: ${CLINVAR_RAW}"
fi

if [[ ! -f "${CLINVAR_NORM}" || ! -f "${CLINVAR_NORM}.tbi" ]]; then
  echo "[3/4] Normalizing ClinVar GRCh38 (split multiallelics + left-align indels)..."
  bcftools view --threads "${THREADS}" -r "${PRIMARY_CHROMS}" -Ou "${CLINVAR_RAW}" \
    | bcftools norm --threads "${THREADS}" -m -any -f "${FASTA_BGZ}" -c w -Oz -o "${CLINVAR_NORM}"
  tabix -f -p vcf "${CLINVAR_NORM}"
  echo "[3/4] Normalized ClinVar records: $(bcftools index -n "${CLINVAR_NORM}")"
else
  echo "[3/4] ClinVar normalized VCF ready: ${CLINVAR_NORM} ($(bcftools index -n "${CLINVAR_NORM}") records)"
fi

if [[ ! -f "${CLINVAR_CLINICAL}" || ! -f "${CLINVAR_CLINICAL}.tbi" ]]; then
  echo "[3/4] Building compact ClinVar clinical subset (P/LP, Conflicting, Risk/Drug)..."
  gzip -dc "${CLINVAR_NORM}" \
    | awk '/^#/ || $8 ~ /(Pathogenic|pathogenic|Conflicting|drug_response|risk_factor|protective|association)/' \
    | bgzip -@ "${THREADS}" -c > "${CLINVAR_CLINICAL}"
  tabix -f -p vcf "${CLINVAR_CLINICAL}"
  echo "[3/4] Compact ClinVar clinical records: $(bcftools index -n "${CLINVAR_CLINICAL}")"
else
  echo "[3/4] Compact ClinVar clinical VCF ready: ${CLINVAR_CLINICAL} ($(bcftools index -n "${CLINVAR_CLINICAL}") records)"
fi

# 5. Dynamic: GenCC Mode of Inheritance table (re-downloaded on --refresh)
GENCC_URL="https://search.thegencc.org/download/action/submissions-export-tsv"
GENCC_TSV="${REFERENCE_DIR}/gencc_submissions.tsv"
if [[ "${REFRESH}" -eq 1 || ! -s "${GENCC_TSV}" ]]; then
  echo "[4/4] Fetching latest GenCC Gene-Disease & Mode of Inheritance table..."
  curl -fL --retry 3 -o "${GENCC_TSV}" "${GENCC_URL}"
  echo "[4/4] GenCC ready: $(wc -l < "${GENCC_TSV}") lines"
else
  echo "[4/4] GenCC ready: ${GENCC_TSV} ($(wc -l < "${GENCC_TSV}") lines)"
fi

# Optional custom clinical VCF normalization
if [[ -n "${CUSTOM_CLINICAL_VCF:-}" && -f "${CUSTOM_CLINICAL_VCF}" ]]; then
  CUSTOM_NORM="${REFERENCE_DIR}/custom_clinical.norm.vcf.gz"
  echo "[Custom] Normalizing ${CUSTOM_CLINICAL_VCF} -> ${CUSTOM_NORM}"
  HAS_CHR="$(bcftools index -s "${CUSTOM_CLINICAL_VCF}" | awk '$1 ~ /^chr/ {print "yes"; exit}')"
  if [[ "${HAS_CHR}" == "yes" ]]; then
    bcftools annotate --threads "${THREADS}" --rename-chrs "${CHR_MAP}" -Ou "${CUSTOM_CLINICAL_VCF}" \
      | bcftools view --threads "${THREADS}" -Ou \
      | bcftools norm --threads "${THREADS}" -m -any -f "${FASTA_BGZ}" -c w -Oz -o "${CUSTOM_NORM}"
  else
    bcftools view --threads "${THREADS}" -r "${PRIMARY_CHROMS}" -Ou "${CUSTOM_CLINICAL_VCF}" \
      | bcftools norm --threads "${THREADS}" -m -any -f "${FASTA_BGZ}" -c w -Oz -o "${CUSTOM_NORM}"
  fi
  tabix -f -p vcf "${CUSTOM_NORM}"
fi

echo "Reference setup complete."
