#!/usr/bin/env bash
set -euo pipefail

# Annual Checkup Runner:
#   1. Ensures references exist (or refreshes ClinVar & GenCC when --refresh-refs is passed)
#   2. Runs the per-sample pipeline for every sample listed in config/samples.tsv
#   3. Runs joint couple screening + auto-generates plain_language_summary.md
#
# Usage:
#   bash scripts/run_annual_checkup.sh                 # Re-run with current local references
#   bash scripts/run_annual_checkup.sh --refresh-refs  # Annual refresh of ClinVar & GenCC first

REFRESH_FLAG=""
for arg in "$@"; do
  case "${arg}" in
    --refresh-refs) REFRESH_FLAG="--refresh" ;;
    *) echo "Unknown option: ${arg} (use --refresh-refs)" >&2; exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

SAMPLES_TSV="config/samples.tsv"
if [[ ! -f "${SAMPLES_TSV}" ]]; then
  echo "ERROR: ${SAMPLES_TSV} not found." >&2
  exit 1
fi

echo "============================================================"
echo "Step A: Checking / Updating Reference Databases..."
bash scripts/setup_references.sh ${REFRESH_FLAG}

MALE_LABEL=""
FEMALE_LABEL=""
while IFS=$'\t' read -r label sample_id sex _rest; do
  [[ -z "${label}" || "${label}" == "label" ]] && continue
  bash scripts/run_sample.sh "${label}"
  if [[ "${sex}" == "M" && -z "${MALE_LABEL}" ]]; then
    MALE_LABEL="${label}"
  elif [[ "${sex}" == "F" && -z "${FEMALE_LABEL}" ]]; then
    FEMALE_LABEL="${label}"
  fi
done < "${SAMPLES_TSV}"

if [[ -n "${MALE_LABEL}" && -n "${FEMALE_LABEL}" ]]; then
  echo "============================================================"
  echo "Step B: Running Joint Couple Screening (${MALE_LABEL} + ${FEMALE_LABEL})..."
  .venv-wsl/bin/python scripts/08_couple_carrier_screen.py --male "${MALE_LABEL}" --female "${FEMALE_LABEL}"
fi

echo "============================================================"
echo "Annual checkup complete! Reports are in results/reports/."
