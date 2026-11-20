"""CPIC Tier-1 pharmacogenomics (PGx) star-allele interpretation & medication card renderer."""

from __future__ import annotations

from typing import Any

from .genomics import normalize_gt_unphased


def _alt_dosage(row: dict[str, Any] | None) -> int | None:
    """Return ALT allele count (0, 1, 2) for a targeted marker row, or None if not callable."""
    if not row:
        return None
    status = str(row.get("status") or "")
    if status in ("NOT_COVERED", "VARIANT_LOW_CONFIDENCE"):
        return None
    gt = normalize_gt_unphased(str(row.get("gt") or ""))
    if gt.startswith("0/0"):
        return 0
    if gt in ("0/1", "1/0"):
        return 1
    if gt == "1/1":
        return 2
    return None


def interpret_pgx_profile(marker_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """
    Map targeted PGx marker rows (from 07_deep_analysis.py) into structured CPIC
    diplotype, phenotype, actionability tier, affected drug classes, and clinical notes.
    """
    by_rs = {str(r["rsid"]): r for r in marker_rows}
    results: list[dict[str, Any]] = []

    def _d(rsid: str) -> int | None:
        return _alt_dosage(by_rs.get(rsid))

    # 1. CYP2C19 (*2, *3, *4 = No Function; *17 = Increased Function)
    d2, d3, d4, d17 = _d("rs4244285"), _d("rs4986893"), _d("rs28399504"), _d("rs12248560")
    if None in (d2, d3, d4, d17):
        results.append({
            "gene": "CYP2C19",
            "diplotype": "Indeterminate (low coverage / QC)",
            "phenotype": "Indeterminate",
            "tier": "INCOMPLETE",
            "drugs": "Clopidogrel, PPIs (omeprazole/pantoprazole), SSRIs (citalopram/escitalopram/sertraline), voriconazole",
            "guidance": "One or more defining CYP2C19 loci lacked PASS >=10x coverage; confirm clinically if prescribing.",
        })
    else:
        lof_alleles: list[str] = ["*2"] * d2 + ["*3"] * d3 + ["*4"] * d4
        inc_alleles: list[str] = ["*17"] * d17
        alleles = (lof_alleles + inc_alleles)[:2]
        while len(alleles) < 2:
            alleles.insert(0, "*1")
        dip = "/".join(sorted(alleles))
        n_lof = len(lof_alleles)
        n_inc = len(inc_alleles)
        if n_lof >= 2:
            pheno, tier = "Poor Metabolizer (PM)", "ACTIONABLE"
            note = "Avoid clopidogrel (use prasugrel/ticagrelor); consider 50% dose reduction for citalopram/escitalopram/sertraline; monitor voriconazole."
        elif n_lof == 1 and n_inc >= 1:
            pheno, tier = "Intermediate Metabolizer (IM) (*17 + LoF)", "ACTIONABLE"
            note = "Avoid clopidogrel for ACS/PCI (alternative antiplatelet recommended); standard/cautious dosing for PPIs & SSRIs."
        elif n_lof == 1:
            pheno, tier = "Intermediate Metabolizer (IM)", "ACTIONABLE"
            note = "Avoid clopidogrel for ACS/PCI (reduced activation; consider prasugrel/ticagrelor); lower PPI maintenance dose may suffice."
        elif n_inc >= 2:
            pheno, tier = "Ultrarapid Metabolizer (UM)", "ACTIONABLE"
            note = "Increased PPI clearance (consider higher initial dose for H. pylori/ulcer); consider alternative to citalopram/escitalopram or voriconazole."
        elif n_inc == 1:
            pheno, tier = "Rapid Metabolizer (RM)", "MODERATE"
            note = "Faster clearance of PPIs and voriconazole; standard clopidogrel activation; monitor response to escitalopram/citalopram."
        else:
            pheno, tier = "Normal Metabolizer (NM)", "NORMAL"
            note = "Standard CPIC dosing for clopidogrel, PPIs, SSRIs, and voriconazole."
        results.append({
            "gene": "CYP2C19",
            "diplotype": dip,
            "phenotype": pheno,
            "tier": tier,
            "drugs": "Clopidogrel, PPIs (omeprazole/pantoprazole), SSRIs (citalopram/escitalopram/sertraline), voriconazole",
            "guidance": note,
        })

    # 2. SLCO1B1 (c.521T>C, rs4149056 -> *5/*15 Decreased Function)
    d_slco = _d("rs4149056")
    if d_slco is None:
        results.append({
            "gene": "SLCO1B1",
            "diplotype": "Indeterminate",
            "phenotype": "Indeterminate",
            "tier": "INCOMPLETE",
            "drugs": "Statins (simvastatin, atorvastatin, rosuvastatin, pravastatin)",
            "guidance": "Locus lacked PASS >=10x coverage.",
        })
    elif d_slco == 2:
        results.append({
            "gene": "SLCO1B1",
            "diplotype": "*5/*5 (or *15/*15; c.521C/C)",
            "phenotype": "Poor Function",
            "tier": "ACTIONABLE",
            "drugs": "Statins (simvastatin, atorvastatin, lovastatin, rosuvastatin)",
            "guidance": "High myopathy risk with simvastatin (avoid; use rosuvastatin/pravastatin/fluvastatin at lower dose per CPIC 2022).",
        })
    elif d_slco == 1:
        results.append({
            "gene": "SLCO1B1",
            "diplotype": "*1/*5 (or *1/*15; c.521T/C)",
            "phenotype": "Decreased Function",
            "tier": "ACTIONABLE",
            "drugs": "Statins (simvastatin, atorvastatin, lovastatin)",
            "guidance": "Increased myopathy risk on simvastatin >20 mg (prefer lower dose or alternative statin such as rosuvastatin/pravastatin).",
        })
    else:
        results.append({
            "gene": "SLCO1B1",
            "diplotype": "*1/*1 (c.521T/T)",
            "phenotype": "Normal Function",
            "tier": "NORMAL",
            "drugs": "Statins (simvastatin, atorvastatin, rosuvastatin)",
            "guidance": "Standard statin starting doses per clinical lipid guidelines.",
        })

    # 3. Warfarin triad: CYP2C9 (*2, *3) + VKORC1 (-1639G>A) + CYP4F2 (*3)
    c2, c3 = _d("rs1799853"), _d("rs1057910")
    if None in (c2, c3):
        results.append({
            "gene": "CYP2C9",
            "diplotype": "Indeterminate",
            "phenotype": "Indeterminate",
            "tier": "INCOMPLETE",
            "drugs": "Warfarin, NSAIDs (celecoxib, ibuprofen, meloxicam, piroxicam), phenytoin",
            "guidance": "Locus lacked PASS >=10x coverage.",
        })
    else:
        c_alleles = (["*2"] * c2 + ["*3"] * c3)[:2]
        while len(c_alleles) < 2:
            c_alleles.insert(0, "*1")
        dip = "/".join(sorted(c_alleles))
        score = sum(1.0 if a == "*1" else (0.5 if a == "*2" else 0.0) for a in c_alleles)
        if score <= 0.5:
            pheno, tier = f"Poor Metabolizer (AS={score:.1f})", "ACTIONABLE"
            note = "Markedly reduced clearance: major warfarin dose reduction; avoid piroxicam/meloxicam/tenoxicam and use lowest effective NSAID dose; 50% phenytoin maintenance reduction."
        elif score <= 1.5:
            pheno, tier = f"Intermediate Metabolizer (AS={score:.1f})", "MODERATE" if score == 1.5 else "ACTIONABLE"
            note = "Reduced warfarin clearance (lower starting dose); for AS=1.0 consider 25-50% lower NSAID/phenytoin starting dose."
        else:
            pheno, tier = "Normal Metabolizer (AS=2.0)", "NORMAL"
            note = "Standard CPIC dosing for warfarin (combine with VKORC1), NSAIDs, and phenytoin."
        results.append({
            "gene": "CYP2C9",
            "diplotype": dip,
            "phenotype": pheno,
            "tier": tier,
            "drugs": "Warfarin, NSAIDs (celecoxib, ibuprofen, meloxicam, piroxicam), phenytoin",
            "guidance": note,
        })

    d_vk = _d("rs9923231")
    if d_vk is not None:
        if d_vk == 2:
            vk_dip, vk_ph, vk_t, vk_n = "A/A (-1639 T/T fwd)", "High Warfarin Sensitivity", "ACTIONABLE", "Requires substantially lower warfarin dose (~40-50% reduction; use CPIC/IWPC warfarin dosing algorithm with CYP2C9)."
        elif d_vk == 1:
            vk_dip, vk_ph, vk_t, vk_n = "G/A (-1639 C/T fwd)", "Intermediate Warfarin Sensitivity", "MODERATE", "Moderately lower warfarin dose requirement (~20-25% reduction; combine with CYP2C9 and clinical factors)."
        else:
            vk_dip, vk_ph, vk_t, vk_n = "G/G (-1639 C/C fwd)", "Normal Warfarin Sensitivity", "NORMAL", "Standard warfarin sensitivity (combine with CYP2C9)."
        results.append({
            "gene": "VKORC1",
            "diplotype": vk_dip,
            "phenotype": vk_ph,
            "tier": vk_t,
            "drugs": "Warfarin",
            "guidance": vk_n,
        })

    d_4f2 = _d("rs2108622")
    if d_4f2 is not None:
        if d_4f2 >= 1:
            f_dip = "*3/*3" if d_4f2 == 2 else "*1/*3"
            results.append({
                "gene": "CYP4F2",
                "diplotype": f_dip,
                "phenotype": "Reduced Vitamin K Metabolism",
                "tier": "MODERATE",
                "drugs": "Warfarin",
                "guidance": "Associated with ~5-10% higher warfarin dose requirement (CPIC optional modifier alongside CYP2C9/VKORC1).",
            })
        else:
            results.append({
                "gene": "CYP4F2",
                "diplotype": "*1/*1",
                "phenotype": "Normal Function",
                "tier": "NORMAL",
                "drugs": "Warfarin",
                "guidance": "No CYP4F2 warfarin dose adjustment indicated.",
            })

    # 4. DPYD (Fluoropyrimidines: 5-FU, capecitabine — CPIC Activity Score)
    d_2a, d_13, d_2846, d_hapb3 = _d("rs3918290"), _d("rs55886062"), _d("rs67376798"), _d("rs75017182")
    if None in (d_2a, d_13, d_2846, d_hapb3):
        results.append({
            "gene": "DPYD",
            "diplotype": "Indeterminate",
            "phenotype": "Indeterminate",
            "tier": "INCOMPLETE",
            "drugs": "Fluoropyrimidines (5-fluorouracil / 5-FU, capecitabine)",
            "guidance": "One or more DPYD loci lacked PASS >=10x coverage; clinical DPYD testing required before 5-FU/capecitabine.",
        })
    else:
        # Each no-function allele (*2A, *13) subtracts 1.0; each decreased-function (c.2846A>T, HapB3) subtracts 0.5
        reduction = min(2.0, 1.0 * (d_2a + d_13) + 0.5 * (d_2846 + d_hapb3))
        as_score = max(0.0, 2.0 - reduction)
        found_vars = (
            (["*2A"] * d_2a) + (["*13"] * d_13) + (["c.2846A>T"] * d_2846) + (["HapB3"] * d_hapb3)
        )
        dip = "/".join(found_vars[:2]) if len(found_vars) >= 2 else (f"Normal/{found_vars[0]}" if found_vars else "Normal/Normal")
        if as_score <= 0.5:
            pheno, tier = f"Poor Metabolizer (AS={as_score:.1f})", "ACTIONABLE"
            note = "CRITICAL SAFETY ALERT: Avoid 5-FU and capecitabine (life-threatening toxicity risk)."
        elif as_score < 2.0:
            pheno, tier = f"Intermediate Metabolizer (AS={as_score:.1f})", "ACTIONABLE"
            note = "Reduce 5-FU / capecitabine starting dose by 50% and titrate based on toxicity/TDM per CPIC."
        else:
            pheno, tier = "Normal Metabolizer (AS=2.0)", "NORMAL"
            note = "Standard fluoropyrimidine dosing (note: does not exclude ultra-rare non-panel DPYD variants)."
        results.append({
            "gene": "DPYD",
            "diplotype": dip,
            "phenotype": pheno,
            "tier": tier,
            "drugs": "Fluoropyrimidines (5-FU, capecitabine)",
            "guidance": note,
        })

    # 5. TPMT (*2, *3B, *3C -> *3A when *3B + *3C in cis)
    t2, t3b, t3c = _d("rs1800462"), _d("rs1800460"), _d("rs1142345")
    if None in (t2, t3b, t3c):
        results.append({
            "gene": "TPMT",
            "diplotype": "Indeterminate",
            "phenotype": "Indeterminate",
            "tier": "INCOMPLETE",
            "drugs": "Thiopurines (azathioprine, 6-mercaptopurine, thioguanine)",
            "guidance": "One or more TPMT loci lacked PASS >=10x coverage.",
        })
    else:
        # Standard CPIC phasing convention: *3B + *3C het is overwhelmingly *1/*3A (cis)
        if t3b == 1 and t3c == 1 and t2 == 0:
            dip, n_def = "*1/*3A (assumed cis; phase unverified)", 1
        elif t3b == 2 and t3c == 2:
            dip, n_def = "*3A/*3A", 2
        else:
            t_alleles = (["*2"] * t2 + ["*3B"] * t3b + ["*3C"] * t3c)[:2]
            n_def = len(t_alleles)
            while len(t_alleles) < 2:
                t_alleles.insert(0, "*1")
            dip = "/".join(sorted(t_alleles))
        if n_def >= 2:
            pheno, tier, note = "Poor Metabolizer (PM)", "ACTIONABLE", "CRITICAL SAFETY ALERT: Avoid thiopurines for non-malignant conditions or reduce dose 10-fold (3x/week) for malignancy."
        elif n_def == 1:
            pheno, tier, note = "Intermediate Metabolizer (IM)", "ACTIONABLE", "Start thiopurines at 30-80% of standard dose per CPIC to prevent myelosuppression."
        else:
            pheno, tier, note = "Normal Metabolizer (NM)", "NORMAL", "Standard thiopurine dosing (also check NUDT15)."
        results.append({
            "gene": "TPMT",
            "diplotype": dip,
            "phenotype": pheno,
            "tier": tier,
            "drugs": "Thiopurines (azathioprine, 6-MP, thioguanine)",
            "guidance": note,
        })

    # 6. NUDT15 (*3)
    d_nudt = _d("rs116855232")
    if d_nudt is not None:
        if d_nudt == 2:
            results.append({
                "gene": "NUDT15", "diplotype": "*3/*3", "phenotype": "Poor Metabolizer (PM)",
                "tier": "ACTIONABLE", "drugs": "Thiopurines (azathioprine, 6-MP, thioguanine)",
                "guidance": "CRITICAL SAFETY ALERT: Avoid thiopurines or drastically reduce dose (severe leukopenia risk).",
            })
        elif d_nudt == 1:
            results.append({
                "gene": "NUDT15", "diplotype": "*1/*3", "phenotype": "Intermediate Metabolizer (IM)",
                "tier": "ACTIONABLE", "drugs": "Thiopurines (azathioprine, 6-MP, thioguanine)",
                "guidance": "Start at 30-80% of standard thiopurine dose per CPIC.",
            })
        else:
            results.append({
                "gene": "NUDT15", "diplotype": "*1/*1", "phenotype": "Normal Metabolizer (NM)",
                "tier": "NORMAL", "drugs": "Thiopurines (azathioprine, 6-MP, thioguanine)",
                "guidance": "Standard thiopurine dosing (combine with TPMT).",
            })

    # 7. CYP3A5 (*3)
    d_3a5 = _d("rs776746")
    if d_3a5 is not None:
        if d_3a5 == 2:
            results.append({
                "gene": "CYP3A5", "diplotype": "*3/*3", "phenotype": "Poor Metabolizer (Non-Expressor)",
                "tier": "NORMAL", "drugs": "Tacrolimus",
                "guidance": "Standard tacrolimus starting dose (most common genotype in European ancestry; use therapeutic drug monitoring).",
            })
        elif d_3a5 == 1:
            results.append({
                "gene": "CYP3A5", "diplotype": "*1/*3", "phenotype": "Intermediate Metabolizer (Expressor)",
                "tier": "ACTIONABLE", "drugs": "Tacrolimus",
                "guidance": "Faster tacrolimus clearance: CPIC recommends 1.5-2x higher starting dose (not exceeding 0.3 mg/kg/day) with TDM.",
            })
        else:
            results.append({
                "gene": "CYP3A5", "diplotype": "*1/*1", "phenotype": "Normal Metabolizer (Expressor)",
                "tier": "ACTIONABLE", "drugs": "Tacrolimus",
                "guidance": "Faster tacrolimus clearance: CPIC recommends 1.5-2x higher starting dose (not exceeding 0.3 mg/kg/day) with TDM.",
            })

    # 8. CYP2B6 (*9/*6 c.516G>T)
    d_2b6 = _d("rs3745274")
    if d_2b6 is not None:
        if d_2b6 == 2:
            results.append({
                "gene": "CYP2B6", "diplotype": "*6/*6 (or *9/*9; c.516T/T)", "phenotype": "Poor Metabolizer (PM)",
                "tier": "ACTIONABLE", "drugs": "Efavirenz, methadone",
                "guidance": "Increased CNS toxicity risk on standard efavirenz 600 mg (CPIC: consider 400 mg or alternative); monitor methadone QTc/sedation.",
            })
        elif d_2b6 == 1:
            results.append({
                "gene": "CYP2B6", "diplotype": "*1/*6 (or *1/*9; c.516G/T)", "phenotype": "Intermediate Metabolizer (IM)",
                "tier": "MODERATE", "drugs": "Efavirenz, methadone",
                "guidance": "Moderately reduced clearance; standard efavirenz or 400 mg if CNS side effects occur.",
            })
        else:
            results.append({
                "gene": "CYP2B6", "diplotype": "*1/*1 (c.516G/G)", "phenotype": "Normal Metabolizer (NM)",
                "tier": "NORMAL", "drugs": "Efavirenz, methadone",
                "guidance": "Standard dosing.",
            })

    # 9. UGT1A1 (*6 only; *28 TA-repeat requires BAM/clinical assay)
    d_ugt = _d("rs4148323")
    if d_ugt is not None:
        if d_ugt >= 1:
            u_dip = "*6/*6" if d_ugt == 2 else "*1/*6"
            results.append({
                "gene": "UGT1A1", "diplotype": f"{u_dip} (*28 unassessed)",
                "phenotype": "Poor Metabolizer" if d_ugt == 2 else "Intermediate Metabolizer",
                "tier": "ACTIONABLE" if d_ugt == 2 else "MODERATE",
                "drugs": "Irinotecan, atazanavir",
                "guidance": "Increased neutropenia/hyperbilirubinemia risk (note: *28 promoter TA-repeat is not resolved by standard VCFs).",
            })
        else:
            results.append({
                "gene": "UGT1A1", "diplotype": "No *6 (*28 unassessed)",
                "phenotype": "Partial — *6 Negative",
                "tier": "PARTIAL",
                "drugs": "Irinotecan, atazanavir",
                "guidance": "Negative for *6; clinical test needed to rule out *28 TA-repeat before irinotecan/atazanavir.",
            })

    # 10. BCHE (Succinylcholine / mivacurium)
    bk, ba = _d("rs1803274"), _d("rs1799807")
    if None not in (bk, ba):
        if ba >= 2 or (ba >= 1 and bk >= 1):
            results.append({
                "gene": "BCHE", "diplotype": f"A-var({ba}) / K-var({bk})",
                "phenotype": "Deficient / Prolonged Apnea Risk", "tier": "ACTIONABLE",
                "drugs": "Succinylcholine, mivacurium (anesthesia neuromuscular blockers)",
                "guidance": "Inform anesthesiologist prior to surgery: increased risk of prolonged neuromuscular blockade/apnea.",
            })
        elif ba == 1 or bk == 2:
            results.append({
                "gene": "BCHE", "diplotype": f"A-var({ba}) / K-var({bk})",
                "phenotype": "Mildly Reduced Activity", "tier": "MODERATE",
                "drugs": "Succinylcholine, mivacurium",
                "guidance": "Mildly prolonged recovery possible with succinylcholine/mivacurium; mention to anesthesia if relevant.",
            })
        else:
            results.append({
                "gene": "BCHE", "diplotype": f"Normal (K-var={bk}, A-var=0)",
                "phenotype": "Normal Activity", "tier": "NORMAL",
                "drugs": "Succinylcholine, mivacurium",
                "guidance": "Standard anesthesia risk for tested variants.",
            })

    # 11. CYP2D6 (*4 only — explicitly marked PARTIAL until BAM Cyrius/Aldy)
    d_2d6 = _d("rs3892097")
    if d_2d6 is not None:
        if d_2d6 >= 1:
            dip_2d6 = "*4/*4 (partial)" if d_2d6 == 2 else "*1/*4 (partial)"
            pheno_2d6 = "Poor Metabolizer (at least)" if d_2d6 == 2 else "At least Intermediate Metabolizer"
            results.append({
                "gene": "CYP2D6", "diplotype": dip_2d6, "phenotype": pheno_2d6,
                "tier": "ACTIONABLE",
                "drugs": "Codeine, tramadol, tamoxifen, TCAs (amitriptyline/nortriptyline), SSRIs (paroxetine/fluoxetine), ondansetron",
                "guidance": "Carries non-functional *4 allele (avoid codeine/tramadol if PM; confirm full star-allele/CNV/hybrid status via BAM Cyrius/Aldy or clinical PGx).",
            })
        else:
            results.append({
                "gene": "CYP2D6", "diplotype": "*4 Negative (CNVs/hybrids unassessed)",
                "phenotype": "Partial — Needs BAM/Clinical PGx", "tier": "PARTIAL",
                "drugs": "Codeine, tramadol, tamoxifen, TCAs, SSRIs, ondansetron",
                "guidance": "Negative for *4, but standard VCFs cannot detect *5 whole-gene deletion, duplications (UM), or CYP2D7 hybrids. Will resolve once BAM/CRAM is added.",
            })

    tier_order = {"ACTIONABLE": 0, "MODERATE": 1, "PARTIAL": 2, "INCOMPLETE": 3, "NORMAL": 4}
    results.sort(key=lambda r: (tier_order.get(r["tier"], 9), r["gene"]))
    return results


def render_pgx_medication_card(
    sample_label: str,
    ref_meta: dict[str, str],
    pgx_entries: list[dict[str, Any]],
) -> str:
    """Render a concise, clinician-friendly CPIC Tier-1 Pharmacogenomics Medication Card."""
    actionable = [e for e in pgx_entries if e["tier"] == "ACTIONABLE"]
    moderate = [e for e in pgx_entries if e["tier"] == "MODERATE"]

    md = [
        f"# Personal Pharmacogenomics (PGx) Medication Card — {sample_label.title()}",
        "",
        f"**Generated:** {ref_meta.get('run_date', '')} | **Reference Basis:** CPIC Tier-1 Guidelines (GRCh38 30x WGS)",
        "",
        "> **Clinical Notice:** Derived from research-grade 30x WGS VCF/gVCF data. Before changing any prescription or dose, "
        "confirm actionable findings via a CLIA-certified pharmacogenomics test and discuss with your prescribing physician or pharmacist.",
        "",
        "## 1. Key Actionable & Moderate Medication Alerts",
        "",
    ]

    if actionable or moderate:
        for e in actionable + moderate:
            badge = "ACTIONABLE" if e["tier"] == "ACTIONABLE" else "MODERATE"
            md.append(
                f"- **[{badge}] `{e['gene']}` — {e['phenotype']} (`{e['diplotype']}`):** "
                f"*{e['drugs']}* — {e['guidance']}"
            )
    else:
        md.append("_No actionable or moderate CPIC Tier-1 variants identified across the tested genes._")

    md.extend([
        "",
        "## 2. Complete CPIC Tier-1 Gene Summary",
        "",
        "| Status | Gene | Diplotype / Genotype | Predicted Phenotype | Key Affected Medications | CPIC Guidance Summary |",
        "|---|---|---|---|---|---|",
    ])
    for e in pgx_entries:
        md.append(
            f"| **{e['tier']}** | **{e['gene']}** | `{e['diplotype']}` | {e['phenotype']} | "
            f"{e['drugs']} | {e['guidance']} |"
        )

    md.extend([
        "",
        "## 3. Known Short-Read VCF Pharmacogenomics Limitations",
        "- ***CYP2D6* structural alleles:** Standard VCFs detect `*4` (`rs3892097`) but cannot reliably call whole-gene deletions (`*5`), gene duplications (ultrarapid metabolizers), or *CYP2D7* hybrids. Uploading BAM/CRAM files unlocks `Cyrius`/`Aldy`.",
        "- ***UGT1A1\\*28* TA-repeat:** Promoter `(TA)7` repeat expansion is not reliably called in standard short-read SNV VCFs.",
        "- ***HLA-A* / *HLA-B* drug hypersensitivity:** Severe hypersensitivity alleles (`HLA-B*57:01` abacavir, `HLA-B*15:02` carbamazepine, `HLA-B*58:01` allopurinol) require dedicated HLA typing on BAM/CRAM reads or clinical HLA testing.",
        "- **Rare non-panel alleles:** A `Normal` call confirms reference (`0/0`) at the major CPIC defining loci, not absence of ultra-rare private variants.",
    ])
    return "\n".join(md) + "\n"
