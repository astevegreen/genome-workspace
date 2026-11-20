# Reproducible 30x WGS Annual Checkup & Carrier Screening Runbook

All commands run from the workspace root in **WSL2 Ubuntu** (or Linux).
Raw source directories configured in `config/samples.tsv` are treated as **strictly read-only**.

## 1. Annual Checkup (Recommended One-Command Workflow)

Once a year (or whenever you want to re-screen against newly curated clinical variants in NCBI ClinVar and GenCC):

```bash
bash scripts/run_annual_checkup.sh --refresh-refs
```

What this does automatically:
1. **`scripts/setup_references.sh --refresh`** — Downloads fresh NCBI ClinVar GRCh38 (with MD5 verification) and GenCC Mode of Inheritance tables, re-normalizes ClinVar against the local GRCh38 FASTA, builds a compact clinical subset (`clinvar_clinical.norm.vcf.gz`) for fast C-level pre-filtering, and clears `api_cache.json` so gnomAD/Ensembl annotations refresh annually (reuses static GRCh38 FASTA and GENCODE v46 CDS BED).
2. **`scripts/run_sample.sh <label>`** — For every sample in `config/samples.tsv`:
   - Compresses/indexes raw VCFs (`snp_indel.genome`, `cnv`, `sv`) if not already present in `data/processed/<label>/`, harmonizing `chr1` → `1` prefixes across all VCFs, and normalizes variants.
   - **Step 4 (`04_carrier_screen.py`):** Single-pass ClinVar GRCh38 (+ optional `custom_clinical.norm.vcf.gz`) indexing + `bcftools view -T` C-level pre-filtering + multi-gene GenCC MOI resolution + ACMG SF v3.3 (84 genes).
   - **Step 5 (`05_annotate_candidates.py`):** Cached gnomAD v4 + Ensembl enrichment (with retry/backoff and 365-day TTL) for P/LP and Conflicting P/LP hits, automated BA1/BS1 artifact & polymorphism flagging, and **unreported rare `HIGH`-impact LoF rescue** across all 84 ACMG SF v3.3 genes (with call-QC and ClinVar Benign/LB exclusion).
   - **Step 7 (`07_deep_analysis.py`):** GENCODE v46 coding-exon CNV/SV filtering (reusing Step 4's ClinVar P/LP gene list and preserving full untruncated gene lists for couple matching), strict `FILTER == PASS` + conservative `MinDP`/`MIN_DP >= 10x` gVCF coverage audit, categorized Clinical / CPIC Tier-1 PGx / Common Trait genotyping (phased- and unphased-safe *APOE*), and **CPIC Medication Card** (`pgx_medication_card.md`).
   - **Step 9 (`09_polygenic_risk.py`, automatic once PRS references are built):** Gold-standard local Polygenic Risk Scores across 9 curated PGS Catalog GRCh38 traits, empirically calibrated against unrelated 1000 Genomes GRCh38 superpopulations via PLINK 2.0 (`polygenic_risk_report.md`).
   - Automatically checks `data/alignments/<label>/` for uploaded `.bam`/`.cram` files.
3. **`scripts/08_couple_carrier_screen.py`** — Runs sex-concordance QC, shared recessive SNV/indel and transparent SV+SV overlap checks (across full untruncated exonic SV gene lists), reciprocal partner coverage & Ensembl VEP novel-variant rescue, maternal X-linked/MT screening, AD/AD-AR & ACMG SF v3.3 screening, and **auto-generates `results/reports/couple/plain_language_summary.md`**.

## 2. Running Individual Steps

```bash
# Refresh references only
bash scripts/setup_references.sh --refresh

# Run a single sample (e.g., partner_a or partner_b from config/samples.tsv)
bash scripts/run_sample.sh partner_a

# Run couple comparison + plain-language summary only
.venv-wsl/bin/python scripts/08_couple_carrier_screen.py --male partner_a --female partner_b

# One-time setup for Gold-Standard 1000G-calibrated Polygenic Risk Scores (~3.8 GB download)
bash scripts/setup_prs_references.sh

# Run automated unit tests
.venv-wsl/bin/python -m pytest tests/
```

## 3. Quality Control & Automated Artifact Filtering
- **Call QC Tiers (`04_carrier_screen.py`):**
  - **`HIGH`**: `FILTER == PASS`, `DP >= 10` (`>= 6` for hemizygous X/Y), `GQ >= 20`, balanced `AB` (`0.25–0.75` for `HET`, `>= 0.85` for `HOM_ALT`/`HEMI_ALT`).
  - **`MODERATE`**: `FILTER == PASS`, slightly outside strict depth/GQ/AB thresholds.
  - **`LOW_CONFIDENCE`**: Failed caller `FILTER`, `DP < 6`, or severely skewed `AB < 0.20`.
- **Automated BA1/BS1 Artifact Flagging (`scripts/lib/genomics.py`):**
  - Variants with gnomAD AF ≥ 5% (BA1), Autosomal Dominant hits with AF ≥ 1% (BS1), benign blood-group traits, or 0★ legacy somatic entries are automatically demoted to a dedicated **Filtered Likely Artifacts & Common Polymorphisms** audit section rather than cluttering actionable tables.
- **Local API Cache (`data/reference/api_cache.json`):**
  - gnomAD v4 and Ensembl VEP responses are cached locally on disk so re-runs are fast and never confuse API timeouts with `"absent"`. Delete `data/reference/api_cache.json` if you ever want to force a fresh remote lookup.

## 4. Planned BAM / CRAM Alignment Workflow (Future Upload)
When raw alignment files are uploaded locally, place them in `data/alignments/<sample>/` (gitignored):
- `data/alignments/partner_a/<sample>.cram` + `.cram.crai` (or `.bam` + `.bam.bai`)
- `data/alignments/partner_b/<sample>.cram` + `.cram.crai`

Having BAM/CRAM alignments locally will unlock:
1. **Repeat Expansions (*FMR1* Fragile X, *HTT*, *C9orf72*, *DMPK*, *FXN*):** `ExpansionHunter`
2. **Spinal Muscular Atrophy (*SMN1*/*SMN2* exon 7 & 2+0 carriers):** `SMNCopyNumberCaller`
3. **Full *CYP2D6* Star-Allele & CNV/Hybrid Calling:** `Cyrius` or `Aldy`
4. **Visual Read-Level Verification & Cis/Trans Phasing:** `samtools tview` / IGV for any shared exonic SV calls (Section 2b of the couple report) or potential compound heterozygotes.
