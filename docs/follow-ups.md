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
4. [ ] **"What changed since last year" report.** Diff each annual run against the previous one: ClinVar
   reclassifications of carried variants, new or dropped findings, and medication-card changes.
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
   external D: drive.
