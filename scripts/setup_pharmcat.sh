#!/usr/bin/env bash
set -euo pipefail

# Install PharmCAT (CPIC star-allele caller) for Step 11's medication-card cross-check, in WSL2/Linux.
#   - Java, Python 3.12 and bcftools/htslib 1.24 come from conda-forge/Bioconda via micromamba
#     (no sudo or apt, same as T1K in setup_bam_tools.sh).
#   - The PharmCAT pipeline (jar + VCF preprocessor) and its positions VCF come from the GitHub release.
#
# Installs into ~/.local/opt/wgs_tools/pharmcat and creates ~/.local/bin/pharmcat_pipeline.

PHARMCAT_VER="${PHARMCAT_VER:-3.4.0}"
TOOLS_DIR="${WGS_TOOLS_DIR:-${HOME}/.local/opt/wgs_tools}"
BIN_DIR="${HOME}/.local/bin"
MAMBA_BIN="${TOOLS_DIR}/bin/micromamba"
PC_ENV="${TOOLS_DIR}/pharmcat-env"
PC_DIR="${TOOLS_DIR}/pharmcat"
RELEASE_URL="https://github.com/PharmGKB/PharmCAT/releases/download/v${PHARMCAT_VER}"
export MAMBA_ROOT_PREFIX="${TOOLS_DIR}/mamba"

mkdir -p "${TOOLS_DIR}/bin" "${BIN_DIR}" "${PC_DIR}"

# 1. micromamba (skipped when setup_bam_tools.sh already installed it)
if [[ ! -x "${MAMBA_BIN}" ]]; then
  echo "[1/4] Installing micromamba..."
  curl -fL --retry 3 -o "${TOOLS_DIR}/micromamba.tar.bz2" "https://micro.mamba.pm/api/micromamba/linux-64/latest"
  # Python's tarfile handles .bz2 itself, so the bzip2 binary (absent on minimal WSL) isn't needed
  python3 - "${TOOLS_DIR}/micromamba.tar.bz2" "${TOOLS_DIR}" <<'PY'
import sys, tarfile
with tarfile.open(sys.argv[1], "r:bz2") as tf:
    tf.extract(tf.getmember("bin/micromamba"), sys.argv[2])
PY
  chmod +x "${MAMBA_BIN}"
  rm -f "${TOOLS_DIR}/micromamba.tar.bz2"
fi

# 2. Java 21, Python 3.12, bcftools/htslib 1.24 (PharmCAT needs Java >= 17 and bcftools/bgzip >= 1.18)
if [[ ! -x "${PC_ENV}/bin/java" || ! -x "${PC_ENV}/bin/bcftools" ]]; then
  echo "[2/4] Creating the PharmCAT environment..."
  "${MAMBA_BIN}" create -y -p "${PC_ENV}" -c conda-forge -c bioconda \
    "openjdk=21" "python=3.12" "bcftools=1.24" "htslib=1.24"
fi

# 3. PharmCAT pipeline release + positions VCF
if [[ ! -s "${PC_DIR}/.installed-${PHARMCAT_VER}" ]]; then
  echo "[3/4] Downloading PharmCAT ${PHARMCAT_VER}..."
  rm -rf "${PC_DIR:?}"/*
  curl -fL --retry 3 -o "${PC_DIR}/pipeline.tar.gz" "${RELEASE_URL}/pharmcat-pipeline-${PHARMCAT_VER}.tar.gz"
  tar -xzf "${PC_DIR}/pipeline.tar.gz" -C "${PC_DIR}"
  rm -f "${PC_DIR}/pipeline.tar.gz"
  for f in "pharmcat_positions_${PHARMCAT_VER}.vcf.bgz" "pharmcat_positions_${PHARMCAT_VER}.vcf.bgz.csi"; do
    curl -fL --retry 3 -o "${PC_DIR}/${f}" "${RELEASE_URL}/${f}"
  done
  PIPELINE_SCRIPT="$(find "${PC_DIR}" -maxdepth 3 -type f -name pharmcat_pipeline | head -n 1)"
  [[ -n "${PIPELINE_SCRIPT}" ]] || { echo "ERROR: pharmcat_pipeline not found in the release tarball" >&2; exit 1; }
  REQS="$(find "${PC_DIR}" -maxdepth 3 -type f -name requirements.txt | head -n 1)"
  if [[ -n "${REQS}" ]]; then
    "${PC_ENV}/bin/python" -m pip install --quiet -r "${REQS}"
  fi
  echo "${PIPELINE_SCRIPT}" > "${PC_DIR}/.installed-${PHARMCAT_VER}"
fi
PIPELINE_SCRIPT="$(cat "${PC_DIR}/.installed-${PHARMCAT_VER}")"

# 4. Wrapper: run the pipeline with this environment's Java, Python and bcftools
echo "[4/4] Writing ${BIN_DIR}/pharmcat_pipeline..."
cat > "${BIN_DIR}/pharmcat_pipeline" <<EOF
#!/usr/bin/env bash
export JAVA_HOME="${PC_ENV}"
export BCFTOOLS_PATH="${PC_ENV}/bin/bcftools"
export BGZIP_PATH="${PC_ENV}/bin/bgzip"
export PATH="${PC_ENV}/bin:\${PATH}"
exec "${PC_ENV}/bin/python" "${PIPELINE_SCRIPT}" "\$@"
EOF
chmod +x "${BIN_DIR}/pharmcat_pipeline"

echo "=== Verification ==="
"${PC_ENV}/bin/java" -version 2>&1 | head -n 1
"${PC_ENV}/bin/bcftools" --version | head -n 1
"${BIN_DIR}/pharmcat_pipeline" --version 2>&1 | head -n 3 || true
ls -la "${PC_DIR}"/pharmcat_positions_*.vcf.bgz
echo "PharmCAT ${PHARMCAT_VER} installed; Step 11 (scripts/11_pharmcat_crosscheck.py) will now run."
