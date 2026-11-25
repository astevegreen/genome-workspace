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
3. [ ] **Sample identity QC.** Check that each BAM and its VCF come from the same person (`bcftools gtcheck` or
   `somalier`) and that genetic sex matches `config/samples.tsv`. A swapped symlink would otherwise put one person's
   BAM-based calls on the other's card.
4. [ ] **"What changed since last year" report.** Diff each annual run against the previous one: ClinVar
   reclassifications of carried variants, new or dropped findings, and medication-card changes.
5. [ ] **UGT1A1\*28 from the BAM.** The TA-repeat promoter allele is missing from the VCF-based card (irinotecan,
   atazanavir). PharmCAT may cover it if the VCF has it; otherwise count TA repeats from the reads.
   Since Step 11, PharmCAT calls UGT1A1 from a gVCF reference block (>= 10x) across the repeat. That is no
   variant call there, not a read-level repeat count, so the read count is still worth doing.
6. [ ] **Storage and backup (owner).** Convert the BAMs to CRAM (about half the size) and keep a second copy off the
   external D: drive.
