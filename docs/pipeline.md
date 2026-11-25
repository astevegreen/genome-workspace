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
   - **Step 3 (`03_sample_identity_qc.py`):** Sample identity QC (`sample_identity_qc.md`). Re-genotypes about 2,700 of the VCF's own autosomal SNVs from the BAM (`bcftools mpileup | call`) and counts agreement: same person < 5% discordant, another person >= 20%; on a mismatch it also tries the BAM against the other samples' VCFs to name a swap. Checks genetic sex from chrX/chrY depth (BAM index) and the VCF's non-PAR chrX het/haploid calls against `config/samples.tsv`. On **FAIL**, Step 10 is skipped, earlier BAM results are kept off the Step 7/8/11 reports, and `run_sample.sh` exits 1 after the VCF-only steps. The BAM check is reused while the alignment link target and VCF are unchanged (`--force` re-runs it).
   - **Step 4 (`04_carrier_screen.py`):** Single-pass ClinVar GRCh38 (+ optional `custom_clinical.norm.vcf.gz`) indexing + `bcftools view -T` C-level pre-filtering + multi-gene GenCC MOI resolution + ACMG SF v3.3 (84 genes).
   - **Step 5 (`05_annotate_candidates.py`):** Cached gnomAD v4 + Ensembl enrichment (with retry/backoff and 365-day TTL) for P/LP and Conflicting P/LP hits, automated BA1/BS1 artifact & polymorphism flagging, and **unreported rare `HIGH`-impact LoF rescue** across all 84 ACMG SF v3.3 genes (with call-QC and ClinVar Benign/LB exclusion).
   - **Step 7 (`07_deep_analysis.py`):** GENCODE v46 coding-exon CNV/SV filtering (reusing Step 4's ClinVar P/LP gene list and preserving full untruncated gene lists for couple matching), strict `FILTER == PASS` + conservative `MinDP`/`MIN_DP >= 10x` gVCF coverage audit, categorized Clinical / CPIC Tier-1 PGx / Common Trait genotyping (phased- and unphased-safe *APOE*; sex-aware X-linked *G6PD* using `sex` from `config/samples.tsv`; *MT-RNR1* aminoglycoside variants with heteroplasmy fraction), and **CPIC Medication Card** (`pgx_medication_card.md`).
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

## 4. BAM / CRAM Alignment Workflow
Place alignment files (or symlinks to them) in `data/alignments/<sample>/` (gitignored). Steps 3 and 10 use the
first `*.bam`, else the first `*.cram`:
- `data/alignments/partner_a/<sample>.cram` + `.cram.crai` (or `.bam` + `.bam.bai`)
- `data/alignments/partner_b/<sample>.cram` + `.cram.crai`

Having BAM/CRAM alignments locally will unlock:
1. **Repeat Expansions (*FMR1* Fragile X, *HTT*, *C9orf72*, *DMPK*, *FXN*):** `ExpansionHunter`
2. **Spinal Muscular Atrophy (*SMN1*/*SMN2* exon 7 & 2+0 carriers):** `SMNCopyNumberCaller`
3. **Full *CYP2D6* Star-Allele & CNV/Hybrid Calling:** `Cyrius` or `Aldy`
4. **HLA Typing (*HLA-A*/*-B* CPIC hypersensitivity alleles `B*57:01`, `B*58:01`, `B*15:02`, `A*31:01`):** `T1K` (`--preset hla-wgs`), installed with its IPD-IMGT/HLA index by `bash scripts/setup_bam_tools.sh`. The same typing also gives celiac HLA-DQ2/DQ8 and HLA-B27 associations (Step 10 report, section 7).
5. ***UGT1A1\*28* Promoter TA Repeat (rs3064744; irinotecan, atazanavir):** counted in Step 10 (section 8) from read pairs that carry both flanking anchors (TA5 = \*36, TA6 = \*1, TA7 = \*28, TA8 = \*37). The diplotype is called only with >= 10 spanning read pairs and a log10 likelihood ratio >= 2 over the next-best pair (2% stutter model); the medication card combines it with the VCF's \*6, and Step 11 compares it with PharmCAT, whose UGT1A1 call otherwise rests on a gVCF reference block across the repeat.
6. **Visual Read-Level Verification & Cis/Trans Phasing:** `samtools tview` / IGV for any shared exonic SV calls (Section 2b of the couple report) or potential compound heterozygotes.

**CRAM input.** `scripts/bam_to_cram.sh <bam>` writes a no-reference CRAM 3.0 next to the BAM's real file and checks
it against the BAM (quickcheck, flagstat, idxstats, `samtools checksum -a`, region counts); it never changes the BAM.
Every Step 3 and Step 10 tool gives the same results on that CRAM as on the BAM. Two need workarounds
(`scripts/lib/alignment.py`):
- ExpansionHunter needs CRAM 3.0 rather than 3.1.
- T1K cannot read CRAM, so Step 10 passes it a temporary BAM of the HLA region, the decoy contigs and the unmapped
  reads.

Why the CRAM is not reference-based, and the backup checklist: [storage_backup.md](storage_backup.md).

## 5. PharmCAT Cross-Check of the Medication Card (Step 11)

One-time install (Java, bcftools 1.24 and PharmCAT 3.4.0 under `~/.local/opt/wgs_tools`):

```bash
bash scripts/setup_pharmcat.sh
```

Step 11 then runs at the end of `run_sample.sh` (or `.venv-wsl/bin/python scripts/11_pharmcat_crosscheck.py --sample <label>`). It:
- builds a PharmCAT input from the sample's variants plus explicit 0/0 calls where the gVCF has a PASS reference block of at least 10x (6x on male non-PAR chrX); every other absent position stays a no-call, never assumed reference;
- passes Step 10's Cyrius *CYP2D6* and untied T1K *HLA-A*/*HLA-B* calls as outside calls;
- compares PharmCAT gene by gene with `pgx_medication_card.tsv` and writes `results/reports/<label>/pharmcat_crosscheck.{md,tsv}`, with PharmCAT's own HTML/JSON report in `pharmcat/`.

PharmCAT also reports genes the card lacks (for example *NAT2*, *IFNL3*, *RYR1*/*CACNA1S*, *CFTR*, *ABCG2*, *CYP3A4*).
