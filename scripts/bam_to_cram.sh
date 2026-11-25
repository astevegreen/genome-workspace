#!/usr/bin/env bash
set -euo pipefail

# Convert BAMs to verified CRAM 3.0 copies. Never modifies, moves or deletes a BAM, and never touches the
# data/alignments/ symlinks: it only reads the BAM and writes <name>.cram (+ .crai, .md5, .verify.txt).
#
# Usage:
#   bash scripts/bam_to_cram.sh [options] <bam> [<bam> ...]
#
# Options:
#   --out-dir DIR      Where to write the CRAMs (default: next to each BAM's real file, following symlinks)
#   --reference FASTA  Reference-based CRAM against FASTA (smaller). The FASTA must contain every @SQ contig of the
#                      BAM with the same length, or the script stops (see below). Default: no reference (no_ref=1).
#   --threads N        samtools threads (default: $THREADS or 4)
#   --force            Replace an existing CRAM of the same name (the BAM is never touched)
#
# Why no_ref by default: these BAMs were aligned to GRCh38 no-alt + hs38d1 decoys (2,581 contigs), and the local
# GRCh38.primary_assembly.fa.bgz lacks the decoy, EBV and UCSC-named unplaced contigs. With such a FASTA, samtools
# silently switches the whole file to embed_ref=2, and htslib 1.22 then writes a .crai whose slices all have span 0:
# quickcheck, flagstat, idxstats and checksums still pass, but region queries miss reads. A no_ref CRAM needs no
# reference to decode, indexes correctly, and is about 66% of the BAM's size (a matching reference gives ~59%).
# CRAM 3.0 rather than samtools' default 3.1: ExpansionHunter v5 cannot decode 3.1 codecs.
#
# Verification, all of which must pass before the .part files are renamed into place:
#   samtools quickcheck; flagstat and idxstats identical between BAM and CRAM; `samtools checksum -a` identical
#   (every record, in order: flags, sequence, qualities, positions, MAPQ/CIGAR, mates and tags); read counts identical in
#   a panel of regions (catches a broken index); then an md5 of the CRAM for checking backup copies.

OUT_DIR=""
REFERENCE=""
THREADS="${THREADS:-4}"
FORCE=0
BAMS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --reference) REFERENCE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) sed -n '4,30p' "$0"; exit 0 ;;
    -*) echo "Unknown option: $1" >&2; exit 2 ;;
    *) BAMS+=("$1"); shift ;;
  esac
done
if [[ ${#BAMS[@]} -eq 0 ]]; then
  echo "Usage: $0 [--out-dir DIR] [--reference FASTA] [--threads N] [--force] <bam> [<bam> ...]" >&2
  exit 2
fi
SAMTOOLS="${SAMTOOLS:-samtools}"
command -v "${SAMTOOLS}" >/dev/null || { echo "ERROR: samtools not found" >&2; exit 1; }
if ! "${SAMTOOLS}" help 2>&1 | grep -q -w checksum; then
  echo "ERROR: ${SAMTOOLS} has no 'checksum' command (needs samtools >= 1.20)" >&2
  exit 1
fi

# Regions whose read counts must agree between BAM and CRAM (GRCh38, no 'chr' prefix; 'chr' is added if the BAM uses it)
CHECK_REGIONS=(1:1000000-1100000 2:233759000-233762000 5:70900000-71000000 6:29900000-31500000 22:42120000-42140000 X:147900000-147960000 Y:2781480-2800000 MT)

log() { echo "[$(date +%H:%M:%S)] $*"; }

convert_one() {
  local bam_in="$1"
  local bam
  bam="$(readlink -f "${bam_in}")"
  [[ -f "${bam}" ]] || { echo "ERROR: ${bam_in} does not resolve to a file (drive mounted?)" >&2; return 1; }
  local stem dir
  stem="$(basename "${bam}" .bam)"
  dir="${OUT_DIR:-$(dirname "${bam}")}"
  mkdir -p "${dir}"
  dir="$(readlink -f "${dir}")"
  local cram="${dir}/${stem}.cram"
  local part="${dir}/${stem}.part.cram"
  local report="${cram}.verify.txt"

  [[ "${cram}" != "${bam}" ]] || { echo "ERROR: output would overwrite the input" >&2; return 1; }
  if [[ -e "${cram}" && "${FORCE}" -ne 1 ]]; then
    echo "ERROR: ${cram} exists (use --force to replace it; the BAM is never touched)" >&2
    return 1
  fi

  log "${stem}: preflight"
  "${SAMTOOLS}" quickcheck -v "${bam}"
  [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" || -f "${bam}.csi" ]] || { echo "ERROR: ${bam} has no index (.bai/.csi)" >&2; return 1; }
  local bam_bytes free_bytes
  bam_bytes="$(stat -c %s "${bam}")"
  free_bytes="$(df -B1 --output=avail "${dir}" | tail -1 | tr -d ' ')"
  if (( free_bytes < bam_bytes * 3 / 4 )); then
    echo "ERROR: ${dir} has $((free_bytes / 1000000000)) GB free; need about $((bam_bytes * 3 / 4 / 1000000000)) GB for the CRAM" >&2
    return 1
  fi

  local fmt_opts=(--output-fmt-option version=3.0)
  if [[ -n "${REFERENCE}" ]]; then
    [[ -f "${REFERENCE}.fai" ]] || { echo "ERROR: ${REFERENCE}.fai not found" >&2; return 1; }
    local missing
    missing="$("${SAMTOOLS}" view -H "${bam}" | awk -F'\t' -v fai="${REFERENCE}.fai" '
      BEGIN { while ((getline l < fai) > 0) { split(l, a, "\t"); len[a[1]] = a[2] } }
      /^@SQ/ { sn = ""; ln = ""; for (i = 2; i <= NF; i++) { if ($i ~ /^SN:/) sn = substr($i, 4); if ($i ~ /^LN:/) ln = substr($i, 4) }
               if (!(sn in len) || len[sn] != ln) n++ }
      END { print n + 0 }')"
    if [[ "${missing}" -gt 0 ]]; then
      echo "ERROR: ${missing} @SQ contig(s) of ${stem} are missing from ${REFERENCE} or have another length." >&2
      echo "       samtools would fall back to embed_ref=2 and write a CRAM index that breaks region queries." >&2
      echo "       Use the exact alignment reference, or drop --reference for a no_ref CRAM." >&2
      return 1
    fi
    fmt_opts+=(--reference "${REFERENCE}")
  else
    fmt_opts+=(--output-fmt-option no_ref=1)
  fi

  local work
  work="$(mktemp -d "${dir}/.${stem}.verify.XXXXXX")"
  log "${stem}: writing ${part}"
  rm -f "${part}" "${part}.crai"
  "${SAMTOOLS}" view -@ "${THREADS}" -C "${fmt_opts[@]}" -o "${part}" "${bam}"
  "${SAMTOOLS}" index -@ "${THREADS}" "${part}"

  local ref_opt=()
  [[ -n "${REFERENCE}" ]] && ref_opt=(--input-fmt-option "reference=${REFERENCE}")
  local failed=0
  log "${stem}: quickcheck"
  "${SAMTOOLS}" quickcheck -v "${part}" || failed=1

  log "${stem}: flagstat (BAM, CRAM)"
  "${SAMTOOLS}" flagstat -@ "${THREADS}" "${bam}" > "${work}/bam.flagstat"
  "${SAMTOOLS}" flagstat -@ "${THREADS}" "${ref_opt[@]}" "${part}" > "${work}/cram.flagstat"
  diff -q "${work}/bam.flagstat" "${work}/cram.flagstat" >/dev/null || { echo "FAIL: flagstat differs" >&2; failed=1; }

  log "${stem}: idxstats (BAM, CRAM)"
  "${SAMTOOLS}" idxstats "${bam}" > "${work}/bam.idxstats"
  "${SAMTOOLS}" idxstats "${ref_opt[@]}" "${part}" > "${work}/cram.idxstats"
  diff -q "${work}/bam.idxstats" "${work}/cram.idxstats" >/dev/null || { echo "FAIL: idxstats differs" >&2; failed=1; }

  log "${stem}: checksum of every record (BAM, CRAM)"
  "${SAMTOOLS}" checksum -a -@ "${THREADS}" "${bam}" | grep -v '^#' > "${work}/bam.checksum"
  "${SAMTOOLS}" checksum -a -@ "${THREADS}" "${ref_opt[@]}" "${part}" | grep -v '^#' > "${work}/cram.checksum"
  [[ -s "${work}/bam.checksum" ]] && diff -q "${work}/bam.checksum" "${work}/cram.checksum" >/dev/null \
    || { echo "FAIL: record checksums differ" >&2; failed=1; }

  log "${stem}: region read counts (index check)"
  local contigs has_chr=0
  contigs="$(cut -f1 "${work}/bam.idxstats")"
  grep -qx "chr1" <<< "${contigs}" && has_chr=1
  # plus one non-main contig (decoy / unplaced) that has reads, whole
  local extra
  extra="$(awk -F'\t' '$3 > 0 && $1 != "*" && $1 !~ /^(chr)?([0-9]+|X|Y|M|MT)$/ {print $1; exit}' "${work}/bam.idxstats")"
  : > "${work}/regions.tsv"
  local r q nb nc
  for r in "${CHECK_REGIONS[@]}" ${extra:+"${extra}"}; do
    q="${r}"
    if [[ "${has_chr}" -eq 1 && "${r}" != "${extra}" ]]; then q="chr${r/#MT/M}"; fi
    grep -qx "${q%%:*}" <<< "${contigs}" || continue  # contig not in this BAM
    nb="$("${SAMTOOLS}" view -c "${bam}" "${q}" 2>/dev/null || echo NA)"
    nc="$("${SAMTOOLS}" view -c "${ref_opt[@]}" "${part}" "${q}" 2>/dev/null || echo NA)"
    printf '%s\t%s\t%s\n' "${q}" "${nb}" "${nc}" >> "${work}/regions.tsv"
    [[ "${nb}" == "${nc}" && "${nb}" != "NA" ]] || { echo "FAIL: ${q}: BAM ${nb} reads, CRAM ${nc}" >&2; failed=1; }
  done

  if [[ "${failed}" -ne 0 ]]; then
    echo "VERIFICATION FAILED for ${stem}: kept ${part} and ${work}/ for inspection; nothing renamed, BAM untouched." >&2
    return 1
  fi

  mv -f "${part}" "${cram}"
  mv -f "${part}.crai" "${cram}.crai"
  log "${stem}: md5 of the CRAM"
  (cd "${dir}" && md5sum "$(basename "${cram}")" "$(basename "${cram}").crai") > "${cram}.md5"

  local cram_bytes
  cram_bytes="$(stat -c %s "${cram}")"
  {
    echo "CRAM verification report"
    echo "date:        $(date -Iseconds)"
    echo "bam:         ${bam} (${bam_bytes} bytes)"
    echo "cram:        ${cram} (${cram_bytes} bytes, $((cram_bytes * 100 / bam_bytes))% of the BAM)"
    echo "mode:        CRAM 3.0, ${REFERENCE:+reference ${REFERENCE}}${REFERENCE:-no_ref (decodes without a reference)}"
    echo "samtools:    $("${SAMTOOLS}" --version | head -1)"
    echo "quickcheck:  PASS"
    echo "flagstat:    identical"
    echo "idxstats:    identical"
    echo "checksum -a: identical ($(wc -l < "${work}/bam.checksum") lines)"
    echo "regions (region, BAM reads, CRAM reads):"
    sed 's/^/  /' "${work}/regions.tsv"
    echo "md5:"
    sed 's/^/  /' "${cram}.md5"
    echo
    echo "flagstat:"
    sed 's/^/  /' "${work}/cram.flagstat"
  } > "${report}"
  rm -rf "${work}"
  log "${stem}: VERIFIED -> ${cram} ($((cram_bytes / 1000000000)) GB, $((cram_bytes * 100 / bam_bytes))% of the BAM); report ${report}"
}

status=0
for b in "${BAMS[@]}"; do
  convert_one "${b}" || status=1
done
exit "${status}"
