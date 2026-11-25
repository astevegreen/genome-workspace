# Follow-ups

Build items agreed on 25 Nov 2026, in order. Tick each off with the commit that ships it.

1. [x] **HLA disease associations from the T1K typing.** Celiac disease (HLA-DQ2.5 / DQ2.2 / DQ8, half-DQ2) and
   HLA-B27, reported as risk associations, not diagnoses. Step 10 report section 7, `hla_disease_associations.tsv`,
   and two rows in the couple plain-language summary.
2. [x] **PharmCAT cross-check of the medication card.** Run PharmCAT on each VCF (with Cyrius CYP2D6 and T1K HLA
   passed in as outside calls), compare its CPIC calls gene by gene with `pgx_medication_card.tsv`, and list
   disagreements. PharmCAT also covers genes the card lacks (NAT2, IFNL3, RYR1/CACNA1S, CFTR for ivacaftor).
   Step 11, `results/reports/<label>/pharmcat_crosscheck.md`. First run on 25 Nov (PharmCAT 3.4.0): no
   disagreements for either sample.
3. [x] **Sample identity QC.** Check that each BAM and its VCF come from the same person (`bcftools gtcheck` or
   `somalier`) and that genetic sex matches `config/samples.tsv`. A swapped symlink would otherwise put one person's
   BAM-based calls on the other's card.
   Step 3, `scripts/03_sample_identity_qc.py`, writes `sample_identity_qc.md`. It uses bcftools, which is already
   installed (somalier would need a downloaded binary and sites file): about 2,700 of the VCF's own SNVs are
   re-genotyped from the BAM with `mpileup | call`, and the agreement is counted directly, because `gtcheck` drops
   sites where the BAM call is 0/0, the very sites a swapped BAM fails on. Sex comes from chrX/chrY depth in the BAM
   index and from the VCF's chrX calls. A failure skips Step 10, keeps BAM results off the Step 7/8/11 reports, and
   makes `run_sample.sh` exit 1. First run on 25 Nov: both BAMs match their VCFs (under 0.5% discordant; a
   simulated swap was 57%) and both sexes match.
4. [x] **"What changed since last year" report.** Diff each annual run against the previous one: ClinVar
   reclassifications of carried variants, new or dropped findings, and medication-card changes.
   Step 13 (`scripts/13_annual_snapshot.py`, `scripts/lib/snapshot.py`), run by `run_annual_checkup.sh`, writes one
   snapshot per sample per run to `data/snapshots/<label>/<run_date>.json`. Each holds the carried P/LP, conflicting and
   VUS variants (looked up in the full ClinVar) with class, review status and stars, plus the LoF candidates,
   medication card, HLA associations and PharmCAT calls. It diffs each snapshot against the previous one into
   `changes_since_last_run.{md,tsv}`. ClinVar is still refreshed in place, since each snapshot keeps that year's
   classes. The first snapshots were written on 25 Nov (about 460 tracked variants per sample); the diff is tested on
   synthetic snapshots.
5. [x] **UGT1A1\*28 from the BAM.** The TA-repeat promoter allele is missing from the VCF-based card (irinotecan,
   atazanavir). PharmCAT may cover it if the VCF has it; otherwise count TA repeats from the reads.
   Since Step 11, PharmCAT calls UGT1A1 from a gVCF reference block (>= 10x) across the repeat. That is no
   variant call there, not a read-level repeat count, so the read count is still worth doing.
   Step 10 section 8 (`scripts/lib/ugt1a1.py`) counts the repeat in read pairs that carry both flanking anchors,
   calls the most likely diplotype under a 2% stutter model (at least 10 spanning pairs, log10 likelihood ratio
   >= 2), and the card combines it with the VCF's \*6. Step 11 notes what each side rests on and flags any
   disagreement. A 30x genome gives only 10-25 spanning pairs, so some samples stay unresolved and need a clinical
   UGT1A1 test.
6. [ ] **Storage and backup (owner).** Convert the BAMs to CRAM (about half the size) and keep a second copy off the
   external D: drive. `scripts/bam_to_cram.sh` now converts and verifies each BAM (no-reference CRAM 3.0, about
   two-thirds of the size), and Steps 3 and 10 read CRAMs. The run checklist and the second-copy steps are in
   [storage_backup.md](storage_backup.md). Tick this once both CRAMs are verified and a verified copy exists off D:.
7. [x] **Stop tracking `pharmcat.log`.** PharmCAT appends to it at the repo root on every Step 11 run (sample labels
   and local paths) and it is tracked in this public repo: `git rm --cached pharmcat.log`, add it to `.gitignore`,
   or have Step 11 run PharmCAT with its log inside `data/processed/<label>/pharmcat/`.
   Done all three ways. PharmCAT's bundled `logback.xml` writes `pharmcat.log` to a relative path (the working
   directory), so Step 11 now runs it with `cwd` set to `data/processed/<label>/pharmcat/`. The file is untracked
   and gitignored. It held only sample labels and local paths, so history was not rewritten.
8. [x] **Re-render the PDFs.** The `.pdf` copies in `results/reports/` date from 22 Nov and are made outside the
   pipeline, so they lack Steps 3 and 11 and the UGT1A1 row.
   The 22 Nov copies came from a one-off scratch script (a hand-written Markdown-to-HTML converter, then Windows
   headless Chrome `--print-to-pdf`) for six reports. Step 12 (`scripts/12_render_pdfs.py`) now renders every report
   `.md` the same way at the end of `run_sample.sh` and for the couple in `run_annual_checkup.sh`. Re-rendered on
   25 Nov for both samples and the couple (16 PDFs).
