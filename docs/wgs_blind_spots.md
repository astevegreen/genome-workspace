# Short-Read 30x WGS Blind Spots & Technical Limitations

Standard short-read (150 bp paired-end) 30x Whole Genome Sequencing has high sensitivity for single-nucleotide variants (SNVs) and small indels in unique regions of the genome, but has **known blind spots** that must be noted alongside any personal carrier screening report. A negative WGS result in these loci does **not** rule out carrier status.

## High-Priority Carrier Screening Blind Spots

| Gene / Locus | Associated Condition | Why Short-Read 30x WGS Struggles | Recommended Clinical Confirmation |
|---|---|---|---|
| **SMN1 / SMN2** | Spinal Muscular Atrophy (SMA) | >99.9% sequence homology between *SMN1* and pseudogene *SMN2*; standard SNV/CNV callers miss exon 7 deletions and "2+0" silent carriers | Dedicated *SMN1* copy-number assay (MLPA / qPCR) or specialized caller (e.g., SMNCopyNumberCaller on BAM/CRAM) |
| **FMR1** | Fragile X Syndrome | CGG trinucleotide repeat expansion in 5' UTR exceeds 150 bp read length | PCR / Southern blot or ExpansionHunter on BAM/CRAM (still needs clinical confirmation) |
| **HBA1 / HBA2** | Alpha-Thalassemia | High homology across the alpha-globin gene cluster; common 3.7 kb / 4.2 kb deletions and triplications are frequently missed by standard VCFs | Gap-PCR / MLPA |
| **CYP21A2** | Congenital Adrenal Hyperplasia (21-hydroxylase deficiency) | Highly homologous pseudogene *CYP21A1P* causes gene conversions and chimeric deletions | Targeted clinical sequencing + MLPA |
| **GBA1 (GBA)** | Gaucher Disease / Parkinson's risk | Pseudogene *GBAP1* (96% homology) causes read multimapping and recombinant alleles | Long-range PCR / clinical panel |
| **PMS2** | Lynch Syndrome | Pseudogene *PMS2CL* overlaps exons 12–15 | Long-range PCR + MLPA |
| **PKD1** | Autosomal Dominant Polycystic Kidney Disease | Exons 1–33 duplicated across 6 pseudogenes | Long-range PCR / clinical panel |
| **CYP2D6** | Pharmacogenomics (drug metabolism) | Complex star-allele copy-number variations, hybrids with *CYP2D7* | Dedicated PGx caller (e.g., Cyrius/Aldy on BAM) or clinical PGx assay |
| **CFTR (Poly-T / TG tract & large del/dup)** | Cystic Fibrosis | Intron 9 poly-T/TG repeat and exon-level deletions can be miscalled without specialized genotyping | Clinical *CFTR* panel |
| **Repeat Expansions (HTT, DMPK, FXN, C9orf72, ATXN*)** | Huntington's, Myotonic Dystrophy, Friedreich Ataxia, ALS/FTD, SCAs | Expansions larger than read length are invisible in standard SNV/indel VCFs | Repeat-primed PCR or ExpansionHunter on BAM/CRAM |

## VCF-Only vs. BAM/CRAM Limitations
- **Phase (cis vs. trans):** Two heterozygous pathogenic variants in the same recessive gene in an unphased VCF could be on the same chromosome (*cis*, carrier only) or opposite chromosomes (*trans*, compound heterozygous). Read-backed phasing (requires BAM/CRAM) or parental testing is needed to resolve phase.
- **No-call / Low-coverage regions:** Absence of a variant in a variant-only VCF does not prove reference genotype (`0/0`); checking reference blocks in the `.genome.vcf` (gVCF) or BAM depth is required to confirm adequate coverage (`DP >= 15–20`).
- **Clinical Disclaimer:** Research-grade bioinformatics pipelines are for educational/informational exploration only. Any actionable or carrier finding should be confirmed via a CLIA-certified clinical laboratory and discussed with a genetic counselor or physician.
