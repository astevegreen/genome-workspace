"""Shared clinical gene panels, targeted markers, and QC/filtering constants."""

from __future__ import annotations

ACMG_SF_VERSION = "v3.3"

# ACMG Secondary Findings v3.3 (84 medically actionable genes; 2025 update)
# Adds ABCD1, CYP27A1, and PLN to ACMG SF v3.2 (81 genes).
ACMG_SF_GENES: frozenset[str] = frozenset({
    "ABCD1", "ACTA2", "ACTC1", "ACVRL1", "APC", "APOB", "ATP7B", "BAG3",
    "BMPR1A", "BRCA1", "BRCA2", "BTD", "CACNA1S", "CALM1", "CALM2", "CALM3",
    "CASQ2", "COL3A1", "CYP27A1", "DES", "DSC2", "DSG2", "DSP", "ENG",
    "FBN1", "FLNC", "GAA", "GLA", "HFE", "HNF1A", "KCNH2", "KCNQ1",
    "LDLR", "LMNA", "MAX", "MEN1", "MLH1", "MSH2", "MSH6", "MUTYH",
    "MYBPC3", "MYH11", "MYH7", "MYL2", "MYL3", "NF2", "OTC", "PALB2",
    "PCSK9", "PKP2", "PLN", "PMS2", "PRKAG2", "PTEN", "RB1", "RBM20",
    "RET", "RPE65", "RYR1", "RYR2", "SCN5A", "SDHAF2", "SDHB", "SDHC",
    "SDHD", "SMAD3", "SMAD4", "STK11", "TGFBR1", "TGFBR2", "TMEM127",
    "TMEM43", "TNNC1", "TNNI3", "TNNT2", "TP53", "TPM1", "TRDN", "TSC1",
    "TSC2", "TTN", "TTR", "VHL", "WT1",
})

PRIMARY_CHROMS: tuple[str, ...] = tuple(
    [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
)
PRIMARY_CHROMS_SET: frozenset[str] = frozenset(PRIMARY_CHROMS)

CLINVAR_STARS: dict[str, int] = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "no_assertion_criteria_provided": 0,
    "no_classification_provided": 0,
    "no_classification_for_the_individual_variant": 0,
    "no_assertion_provided": 0,
}

# Core reproductive carrier & short-read WGS blind-spot genes for coverage auditing.
# Sample-specific carrier hits are added dynamically at runtime.
CARRIER_AND_BLINDSPOT_GENES: dict[str, str] = {
    "CFTR": "Carrier_Core",
    "GJB2": "Carrier_Core",
    "PAH": "Carrier_Core",
    "HEXA": "Carrier_Core",
    "HBB": "Carrier_Core",
    "SERPINA1": "Carrier_Core",
    "ASPA": "Carrier_Core",
    "BLM": "Carrier_Core",
    "DHCR7": "Carrier_Core",
    "FANCC": "Carrier_Core",
    "ELP1": "Carrier_Core",
    "MCOLN1": "Carrier_Core",
    "MEFV": "Carrier_Core",
    "ABCC8": "Carrier_Core",
    "ACADM": "Carrier_Core",
    "SMN1": "BlindSpot_Homology",
    "HBA1": "BlindSpot_Homology",
    "HBA2": "BlindSpot_Homology",
    "FMR1": "BlindSpot_Repeat",
    "CYP21A2": "BlindSpot_Pseudogene",
    "GBA1": "BlindSpot_Pseudogene",
    "PKD1": "BlindSpot_Pseudogene",
    "CYP2D6": "BlindSpot_PGx",
    "PRSS1": "BlindSpot_Pseudogene",
}

# Automated BA1 / BS1 population frequency thresholds for Mendelian P/LP flagging
BA1_COMMON_AF = 0.05       # >= 5% allele frequency: stand-alone benign / artifact for Mendelian P/LP
AD_BS1_AF = 0.01           # >= 1% allele frequency for Autosomal Dominant P/LP: likely benign/artifact/low-penetrance
NOVEL_LOF_MAX_AF = 0.001   # < 0.1% gnomAD AF for novel unreported ACMG LoF rescue

# Known benign/polymorphic SV/CNV loci frequently called in short-read WGS
KNOWN_POLYMORPHIC_SV_GENES: dict[str, str] = {
    "KANSL1": "Common 17q21.31 H2 inversion/duplication polymorphism (~20% European AF)",
    "GSTT2B": "Common polymorphic glutathione S-transferase deletion",
    "HLA-A": "Hypervariable MHC/HLA region (short-read SV calls typically reflect normal HLA haplotype diversity)",
    "HLA-DRB1": "Hypervariable MHC/HLA region (normal HLA-DRB haplotype copy-number polymorphism)",
    "HLA-DQB1": "Hypervariable MHC/HLA region (normal HLA haplotype diversity)",
    "CLCNKB": "Adjacent highly homologous CLCNKA/CLCNKB segmental duplication region",
    "PRSS1": "Highly homologous PRSS1/PRSS2/PRSS3P2 trypsinogen locus (short-read mapping artifact)",
}

# Canonical clinical, pharmacogenomic (CPIC Tier-1), and common-trait markers.
# Tuple format: (category, gene, rsid, chrom, pos_1based, ref, alt, description)
# All GRCh38 coordinates and REF/ALT alleles are forward-strand primary assembly.
TARGETED_MARKERS: list[tuple[str, str, str, str, int, str, str, str]] = [
    # --- Clinical Risk / Carrier Markers ---
    ("clinical", "APOE", "rs429358", "19", 44908684, "T", "C", "APOE ε4 defining (p.Cys130Arg)"),
    ("clinical", "APOE", "rs7412", "19", 44908822, "C", "T", "APOE ε2 defining (p.Arg176Cys)"),
    ("clinical", "HFE", "rs1800562", "6", 26092913, "G", "A", "Hereditary hemochromatosis C282Y"),
    ("clinical", "HFE", "rs1799945", "6", 26090951, "C", "G", "Hereditary hemochromatosis H63D (mild modifier)"),
    ("clinical", "F5", "rs6025", "1", 169549811, "C", "T", "Factor V Leiden thrombophilia"),
    ("clinical", "F2", "rs1799963", "11", 46739505, "G", "A", "Prothrombin G20210A thrombophilia"),
    ("clinical", "SERPINA1", "rs28929474", "14", 94378610, "C", "T", "Alpha-1 antitrypsin deficiency PiZ"),
    ("clinical", "SERPINA1", "rs17580", "14", 94380925, "T", "A", "Alpha-1 antitrypsin deficiency PiS"),
    ("clinical", "CFTR", "rs113993960", "7", 117559590, "ATCT", "A", "Cystic fibrosis ΔF508"),
    ("clinical", "GJB2", "rs80338939", "13", 20189546, "AC", "A", "Nonsyndromic hearing loss 35delG"),
    # --- Pharmacogenomics (CPIC Tier-1 & Key Drug Safety) ---
    ("pgx", "SLCO1B1", "rs4149056", "12", 21178615, "T", "C", "Statin-induced myopathy (*5/*15)"),
    ("pgx", "VKORC1", "rs9923231", "16", 31096368, "C", "T", "Warfarin sensitivity (-1639G>A)"),
    ("pgx", "CYP2C19", "rs4244285", "10", 94781859, "G", "A", "CYP2C19*2 (clopidogrel/PPI no function)"),
    ("pgx", "CYP2C19", "rs4986893", "10", 94780653, "G", "A", "CYP2C19*3 (no function)"),
    ("pgx", "CYP2C19", "rs28399504", "10", 94762706, "A", "G", "CYP2C19*4 (no function)"),
    ("pgx", "CYP2C19", "rs12248560", "10", 94761900, "C", "T", "CYP2C19*17 (increased function)"),
    ("pgx", "CYP2C9", "rs1799853", "10", 94942290, "C", "T", "CYP2C9*2 (warfarin/NSAID reduced function)"),
    ("pgx", "CYP2C9", "rs1057910", "10", 94981296, "A", "C", "CYP2C9*3 (warfarin/NSAID reduced function)"),
    ("pgx", "CYP4F2", "rs2108622", "19", 15879621, "C", "T", "CYP4F2*3 (warfarin higher-dose modifier)"),
    ("pgx", "CYP3A5", "rs776746", "7", 99672916, "T", "C", "CYP3A5*3 (tacrolimus non-expressor / no function)"),
    ("pgx", "CYP2B6", "rs3745274", "19", 41006936, "G", "T", "CYP2B6*9/*6 (efavirenz/methadone decreased function)"),
    ("pgx", "UGT1A1", "rs4148323", "2", 233760498, "G", "A", "UGT1A1*6 (irinotecan/atazanavir decreased function; *28 TA-repeat needs BAM/clinical PGx)"),
    ("pgx", "CYP2D6", "rs3892097", "22", 42128945, "C", "T", "CYP2D6*4 (no function; verify CNVs/hybrids with BAM/clinical PGx)"),
    ("pgx", "DPYD", "rs3918290", "1", 97450058, "C", "T", "DPYD*2A (fluoropyrimidine severe toxicity)"),
    ("pgx", "DPYD", "rs55886062", "1", 97515787, "A", "C", "DPYD*13 (fluoropyrimidine severe toxicity)"),
    ("pgx", "DPYD", "rs67376798", "1", 97082391, "T", "A", "DPYD c.2846A>T (fluoropyrimidine toxicity)"),
    ("pgx", "DPYD", "rs75017182", "1", 97579893, "G", "C", "DPYD HapB3 causal splice variant (fluoropyrimidine toxicity)"),
    ("pgx", "TPMT", "rs1800462", "6", 18143724, "C", "G", "TPMT*2 (thiopurine toxicity)"),
    ("pgx", "TPMT", "rs1800460", "6", 18138997, "C", "T", "TPMT*3B (thiopurine toxicity)"),
    ("pgx", "TPMT", "rs1142345", "6", 18130687, "T", "C", "TPMT*3C (thiopurine toxicity)"),
    ("pgx", "NUDT15", "rs116855232", "13", 48045719, "C", "T", "NUDT15*3 (thiopurine toxicity)"),
    ("pgx", "BCHE", "rs1803274", "3", 165773492, "C", "T", "BCHE K-variant (succinylcholine/mivacurium)"),
    ("pgx", "BCHE", "rs1799807", "3", 165830741, "T", "C", "BCHE A-variant (succinylcholine apnea risk)"),
    # --- Common / Non-Actionable Traits (Informational Only) ---
    ("common_trait", "MTHFR", "rs1801133", "1", 11796321, "G", "A", "MTHFR C677T (common folate variant; ACMG considers non-actionable alone)"),
    ("common_trait", "MTHFR", "rs1801131", "1", 11794419, "T", "G", "MTHFR A1298C (common folate variant; ACMG considers non-actionable alone)"),
    ("common_trait", "PNPLA3", "rs738409", "22", 43928847, "C", "G", "PNPLA3 I148M (common hepatic steatosis / fatty liver modifier)"),
]
