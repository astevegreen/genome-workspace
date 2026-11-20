#!/usr/bin/env bash
set -euo pipefail

# Gold-Standard PRS Reference Setup
# Downloads:
#   1. PLINK 2.0 binary (into .venv-wsl/bin/plink2 if missing)
#   2. 1000 Genomes Project Phase 3 GRCh38 high-coverage PLINK2 panel (~3.8 GB one-time)
#   3. Curated PGS Catalog GRCh38 harmonized scoring files (9 clinical traits)
# Then builds cached 1000G empirical reference distributions via scripts/09_polygenic_risk.py --build-ref.
#
# Usage:
#   bash scripts/setup_prs_references.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

PRS_DIR="data/reference/prs"
KG_DIR="${PRS_DIR}/1000g_hg38"
SCORES_DIR="${PRS_DIR}/scores"
mkdir -p "${KG_DIR}" "${SCORES_DIR}"

PLINK2=".venv-wsl/bin/plink2"
if [[ ! -x "${PLINK2}" ]]; then
  echo "[1/4] Installing PLINK 2.0 into ${PLINK2}..."
  curl -fL --retry 3 -o /tmp/plink2.zip "https://s3.amazonaws.com/plink2-assets/beta1/plink2_linux_x86_64_20261107.zip"
  python3 -m zipfile -e /tmp/plink2.zip /tmp/plink2_extract
  mv /tmp/plink2_extract/plink2 "${PLINK2}"
  chmod +x "${PLINK2}"
  rm -rf /tmp/plink2.zip /tmp/plink2_extract
else
  echo "[1/4] PLINK 2.0 ready: $("${PLINK2}" --version)"
fi

# 2. 1000 Genomes GRCh38 PLINK2 reference panel (official PLINK 2.0 resources)
PGEN_ZST="${KG_DIR}/all_hg38.pgen.zst"
PGEN="${KG_DIR}/all_hg38.pgen"
PVAR_ZST="${KG_DIR}/all_hg38_noannot.pvar.zst"
PSAM="${KG_DIR}/hg38_corrected.psam"
KING_OUT="${KG_DIR}/deg2_hg38.king.cutoff.out.id"

if [[ ! -f "${PGEN}" ]]; then
  if [[ ! -f "${PGEN_ZST}" ]]; then
    echo "[2/4] Downloading 1000 Genomes GRCh38 pgen (~3.2 GB, one-time)..."
    curl -fL --retry 5 -C - -o "${PGEN_ZST}" "https://www.dropbox.com/s/e5n8yr4n7y91fyp/all_hg38.pgen.zst?dl=1"
  fi
  echo "[2/4] Decompressing all_hg38.pgen.zst via plink2..."
  "${PLINK2}" --zst-decompress "${PGEN_ZST}" > "${PGEN}"
  rm -f "${PGEN_ZST}"
else
  echo "[2/4] 1000G GRCh38 pgen ready: ${PGEN}"
fi

if [[ ! -f "${PVAR_ZST}" ]]; then
  echo "[2/4] Downloading 1000 Genomes GRCh38 pvar (~349 MB, one-time)..."
  curl -fL --retry 5 -C - -o "${PVAR_ZST}" "https://www.dropbox.com/s/n6ko1td6ltd81bv/all_hg38_noannot.pvar.zst?dl=1"
else
  echo "[2/4] 1000G GRCh38 pvar ready: ${PVAR_ZST}"
fi

if [[ ! -f "${PSAM}" ]]; then
  echo "[2/4] Downloading 1000 Genomes GRCh38 psam..."
  curl -fL --retry 3 -o "${PSAM}" "https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1"
fi

if [[ ! -f "${KING_OUT}" ]]; then
  echo "[2/4] Downloading 1000 Genomes unrelated cutoff list..."
  curl -fL --retry 3 -o "${KING_OUT}" "https://www.dropbox.com/s/4zhmxpk5oclfplp/deg2_hg38.king.cutoff.out.id?dl=1"
fi

# 3. Curated PGS Catalog GRCh38 harmonized scoring files
PGS_IDS=(
  "PGS000018"  # Coronary artery disease (Inouye metaGRS, 1.7M SNPs)
  "PGS000014"  # Type 2 diabetes (Khera et al., 6.9M SNPs)
  "PGS000016"  # Atrial fibrillation (Khera et al., 6.7M SNPs)
  "PGS000688"  # LDL cholesterol (16k SNPs)
  "PGS002257"  # Systolic blood pressure (884 SNPs)
  "PGS000004"  # Breast cancer (Mavaddat 313 SNPs)
  "PGS000662"  # Prostate cancer (Conti multi-ancestry 269 SNPs)
  "PGS000734"  # Colorectal cancer (95 SNPs)
  "PGS000334"  # Late-onset Alzheimer's disease (22 SNPs)
)

echo "[3/4] Downloading curated PGS Catalog GRCh38 harmonized scoring files..."
for pid in "${PGS_IDS[@]}"; do
  out="${SCORES_DIR}/${pid}_hmPOS_GRCh38.txt.gz"
  if [[ ! -s "${out}" ]]; then
    echo "  Fetching ${pid}..."
    curl -fL --retry 3 -o "${out}" \
      "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/${pid}/ScoringFiles/Harmonized/${pid}_hmPOS_GRCh38.txt.gz"
  else
    echo "  Ready: ${out}"
  fi
done

# 4. Build cached 1000G empirical reference distributions & ancestry-informative markers
if [[ -f "scripts/09_polygenic_risk.py" ]]; then
  echo "[4/4] Building 1000 Genomes empirical PRS reference distributions..."
  .venv-wsl/bin/python scripts/09_polygenic_risk.py --build-ref
fi

echo "PRS reference setup complete."
