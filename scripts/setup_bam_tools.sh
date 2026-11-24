#!/usr/bin/env bash
set -euo pipefail

# Install and configure BAM-based blind-spot callers in WSL2/Linux:
#   1. ExpansionHunter v5.0.0 (FMR1, HTT, C9orf72, DMPK, FXN, ATXN*, etc.)
#   2. SMNCopyNumberCaller (SMN1/SMN2 exon 7 copy number & 2+0 silent carrier)
#   3. Cyrius (full CYP2D6 star-allele, CNV, and CYP2D7 hybrid caller)
#
# Installs into ~/.local/opt/wgs_tools and creates wrappers in ~/.local/bin.

TOOLS_DIR="${WGS_TOOLS_DIR:-${HOME}/.local/opt/wgs_tools}"
BIN_DIR="${HOME}/.local/bin"
export PATH="${BIN_DIR}:${PATH}"

mkdir -p "${TOOLS_DIR}" "${BIN_DIR}"
cd "${TOOLS_DIR}"

# 1. ExpansionHunter v5.0.0
EH_VER="5.0.0"
EH_DIR="${TOOLS_DIR}/ExpansionHunter-v${EH_VER}-linux_x86_64"
if [[ ! -x "${EH_DIR}/bin/ExpansionHunter" ]]; then
  echo "[1/4] Downloading ExpansionHunter v${EH_VER}..."
  curl -fL --retry 3 -o "ExpansionHunter-v${EH_VER}-linux_x86_64.tar.gz" \
    "https://github.com/Illumina/ExpansionHunter/releases/download/v${EH_VER}/ExpansionHunter-v${EH_VER}-linux_x86_64.tar.gz"
  tar -xzf "ExpansionHunter-v${EH_VER}-linux_x86_64.tar.gz"
  rm -f "ExpansionHunter-v${EH_VER}-linux_x86_64.tar.gz"
fi
ln -sfn "${EH_DIR}/bin/ExpansionHunter" "${BIN_DIR}/ExpansionHunter"

# Generate no-chr GRCh38 variant catalog alongside the bundled chr-prefixed catalog
EH_CAT_CHR="${EH_DIR}/variant_catalog/grch38/variant_catalog.json"
EH_CAT_NOCHR="${EH_DIR}/variant_catalog/grch38/variant_catalog.nochr.json"
python3 - "${EH_CAT_CHR}" "${EH_CAT_NOCHR}" <<'PY'
import json, re, sys

src, dst = sys.argv[1], sys.argv[2]
with open(src, encoding="utf-8") as f:
    data = json.load(f)

def strip_chr(val):
    if isinstance(val, str):
        val = re.sub(r"^chrM:", "MT:", val)
        val = re.sub(r"^chr([0-9XY]+):", r"\1:", val)
        return val
    if isinstance(val, list):
        return [strip_chr(x) for x in val]
    return val

for entry in data:
    if "ReferenceRegion" in entry:
        entry["ReferenceRegion"] = strip_chr(entry["ReferenceRegion"])
    if "OfftargetRegions" in entry:
        entry["OfftargetRegions"] = strip_chr(entry["OfftargetRegions"])

with open(dst, "w", encoding="utf-8") as f:
    json.dump(data, f, indent=4)
print(f"Wrote {len(data)} repeat loci to {dst}")
PY

# 2. SMNCopyNumberCaller
if [[ ! -d "${TOOLS_DIR}/SMNCopyNumberCaller" ]]; then
  echo "[2/4] Cloning Illumina/SMNCopyNumberCaller..."
  git clone --depth 1 https://github.com/Illumina/SMNCopyNumberCaller.git "${TOOLS_DIR}/SMNCopyNumberCaller"
fi

# 3. Cyrius (CYP2D6)
if [[ ! -d "${TOOLS_DIR}/Cyrius" ]]; then
  echo "[3/4] Cloning Illumina/Cyrius..."
  git clone --depth 1 https://github.com/Illumina/Cyrius.git "${TOOLS_DIR}/Cyrius"
fi

# Patch open_alignment_file in SMNCopyNumberCaller & Cyrius to harmonize chr/no-chr contigs automatically
python3 - "${TOOLS_DIR}/SMNCopyNumberCaller/depth_calling/utilities.py" "${TOOLS_DIR}/Cyrius/depth_calling/utilities.py" <<'PY'
import pathlib, sys

WRAPPER = '''

class _ContigHarmonizingBam:
    """Proxy pysam.AlignmentFile that harmonizes chr/no-chr contig names on fetch/pileup."""
    def __init__(self, af):
        self._af = af
        self._refs = set(af.references)

    def _resolve_contig(self, contig):
        if contig is None or contig in self._refs:
            return contig
        if isinstance(contig, str):
            if contig.startswith("chr") and contig[3:] in self._refs:
                return contig[3:]
            if not contig.startswith("chr") and ("chr" + contig) in self._refs:
                return "chr" + contig
            if contig == "chrM" and "MT" in self._refs:
                return "MT"
            if contig == "MT" and "chrM" in self._refs:
                return "chrM"
        return contig

    def fetch(self, contig=None, *args, **kwargs):
        if "contig" in kwargs:
            kwargs["contig"] = self._resolve_contig(kwargs["contig"])
            return self._af.fetch(*args, **kwargs)
        if contig is not None:
            return self._af.fetch(self._resolve_contig(contig), *args, **kwargs)
        return self._af.fetch(*args, **kwargs)

    def pileup(self, contig=None, *args, **kwargs):
        if "contig" in kwargs:
            kwargs["contig"] = self._resolve_contig(kwargs["contig"])
            return self._af.pileup(*args, **kwargs)
        if contig is not None:
            return self._af.pileup(self._resolve_contig(contig), *args, **kwargs)
        return self._af.pileup(*args, **kwargs)

    def __getattr__(self, name):
        return getattr(self._af, name)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return self._af.__exit__(exc_type, exc, tb)
'''

for p_str in sys.argv[1:]:
    p = pathlib.Path(p_str)
    if not p.exists():
        continue
    txt = p.read_text(encoding="utf-8")
    if "_ContigHarmonizingBam" in txt:
        continue
    txt = txt + WRAPPER
    txt = txt.replace(
        "return pysam.AlignmentFile(",
        "return _ContigHarmonizingBam(pysam.AlignmentFile(",
    )
    # Close the extra parenthesis on each return _ContigHarmonizingBam(pysam.AlignmentFile(...) block
    lines = txt.splitlines()
    out_lines = []
    in_open = False
    depth = 0
    for line in lines:
        if "return _ContigHarmonizingBam(pysam.AlignmentFile(" in line:
            in_open = True
            depth = line.count("(") - line.count(")")
            if depth == 1:
                line = line + ")"
                in_open = False
                depth = 0
        elif in_open:
            depth += line.count("(") - line.count(")")
            if depth == 1:
                line = line + ")"
                in_open = False
                depth = 0
        out_lines.append(line)
    p.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    print(f"Patched contig harmonization in {p}")
PY

# 4. Shared Python venv for SMNCopyNumberCaller + Cyrius
echo "[4/4] Configuring Python environment for SMNCopyNumberCaller & Cyrius..."
if [[ ! -d "${TOOLS_DIR}/.venv" ]]; then
  if command -v uv >/dev/null 2>&1; then
    uv venv --python 3.12 "${TOOLS_DIR}/.venv"
  else
    python3 -m venv "${TOOLS_DIR}/.venv"
  fi
fi

if command -v uv >/dev/null 2>&1; then
  uv pip install --python "${TOOLS_DIR}/.venv/bin/python" numpy scipy pysam statsmodels reportlab
else
  "${TOOLS_DIR}/.venv/bin/pip" install --upgrade pip
  "${TOOLS_DIR}/.venv/bin/pip" install numpy scipy pysam statsmodels reportlab
fi

cat > "${BIN_DIR}/smn_caller.py" <<EOF
#!/usr/bin/env bash
exec "${TOOLS_DIR}/.venv/bin/python" "${TOOLS_DIR}/SMNCopyNumberCaller/smn_caller.py" "\$@"
EOF
chmod +x "${BIN_DIR}/smn_caller.py"

cat > "${BIN_DIR}/star_caller.py" <<EOF
#!/usr/bin/env bash
exec "${TOOLS_DIR}/.venv/bin/python" "${TOOLS_DIR}/Cyrius/star_caller.py" "\$@"
EOF
chmod +x "${BIN_DIR}/star_caller.py"

echo "=== Verification ==="
"${BIN_DIR}/ExpansionHunter" --version
"${BIN_DIR}/smn_caller.py" --help | head -n 10
"${BIN_DIR}/star_caller.py" --help | head -n 10
echo "All BAM blind-spot callers installed successfully."
