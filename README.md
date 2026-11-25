# Personal 30x WGS Bioinformatics Workspace

Reproducible local workspace for annual personal & couple carrier screening, ACMG SF v3.3 actionable variant analysis, and pharmacogenomics (CPIC genes including *G6PD*, *MT-RNR1*, and BAM-based *CYP2D6* and HLA typing), plus HLA disease associations (celiac HLA-DQ2/DQ8, HLA-B27) and a PharmCAT cross-check of the medication card, on 30x Whole Genome Sequencing (GRCh38) data, with BAM/CRAM alignment support.

## Quick Start: Annual Checkup (WSL2 Ubuntu)

```bash
# Re-run the full pipeline with the latest ClinVar & GenCC databases:
bash scripts/run_annual_checkup.sh --refresh-refs

# Or re-run with currently downloaded local references:
bash scripts/run_annual_checkup.sh
```

Reports are written to `results/reports/<sample>/` and `results/reports/couple/` (including the auto-generated `plain_language_summary.md`).

## Privacy & Read-Only Policy
- Raw source directories configured in `config/samples.tsv` are treated as **strictly read-only**.
- `.gitignore` excludes `data/`, `results/`, `config/samples.tsv`, `config/paths.env`, and all genomic/report file formats (`*.vcf*`, `*.bam*`, `*.cram*`, `*.fastq*`, `*.tsv`, `*.csv`, `*.pdf`, `*.html`).
- Never commit personal variant data or clinical reports.

## Directory Structure
- `config/` — Sample registry (`samples.tsv.example` committed; `samples.tsv` gitignored).
- `data/processed/<sample>/` — `bgzip`-compressed, `tabix`-indexed, and normalized personal VCFs.
- `data/alignments/<sample>/` — BAM or CRAM + index (symlinks to the external drive). `scripts/bam_to_cram.sh` converts a BAM to a verified CRAM; see [docs/storage_backup.md](docs/storage_backup.md) for the conversion and backup checklist.
- `data/reference/` — ClinVar GRCh38, GenCC MOI, GRCh38 FASTA, GENCODE v46 CDS BED, and local API cache (`api_cache.json`).
- `scripts/` — Reproducible pipeline entry points (`run_annual_checkup.sh`, `setup_references.sh`, `run_sample.sh`, `04`–`08`) and shared Python library (`scripts/lib/`).
- `tests/` — Automated unit tests (`pytest`).
- `results/qc/<sample>/` — Header snapshots and contig summaries.
- `results/reports/` — Per-sample and couple screening reports.
- `docs/` — Pipeline runbook ([docs/pipeline.md](docs/pipeline.md)) and short-read WGS blind-spot guide ([docs/wgs_blind_spots.md](docs/wgs_blind_spots.md)).

## Environment (WSL2 Ubuntu)
- Native bioinformatics tools: `bcftools`, `tabix` (`bgzip`/`tabix`), `samtools`, `bedtools`.
- Python environment: managed with `uv` (`.venv-wsl/`, Python 3.12 with `pysam`, `pandas`, `requests`, `pytest`).
