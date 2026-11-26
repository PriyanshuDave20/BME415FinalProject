# Mutation Impact Summary: R168C vs WT Across All RNA Targets

## Overview

This report summarizes protein-RNA contact analysis for KHDRBS2 WT and R168C mutant
across four RNA targets: NRXN1, NRXN2, NRXN3, and AR-V7.

- **Docking method**: Rigid-body docking (HDOCKlite-compatible format)
- **Models per condition**: 10
- **Cutoff distances analyzed**: 4.0 A and 5.0 A
- **Clash filter**: contacts with minimum distance < 2.0 A excluded
- **Primary comparison cutoff**: 5.0 A

---

## Per-Target Contact Summary (5.0 A cutoff, median across 10 models)

| Target | Condition | Median Viable Residues | Median Unique RNA | Median Viable Pairs | Total Unique Residues |
|--------|-----------|----------------------|-------------------|--------------------|----------------------|
| NRXN1 | WT | 6.5 | 8.5 | 13.0 | 86 |
| NRXN1 | R168C | 27.0 | 26.5 | 52.5 | 163 |
| NRXN2 | WT | 4.5 | 8.0 | 12.5 | 53 |
| NRXN2 | R168C | 6.0 | 8.5 | 13.0 | 55 |
| NRXN3 | WT | 4.5 | 8.0 | 10.5 | 37 |
| NRXN3 | R168C | 8.0 | 8.0 | 18.0 | 58 |
| AR-V7 | WT | 4.5 | 6.0 | 10.0 | 43 |
| AR-V7 | R168C | 6.0 | 7.0 | 11.5 | 53 |

---

## R168C / WT Fold Change (5.0 A cutoff)

| Target | Residues Fold Change | RNA Fold Change | Pairs Fold Change | Total Residues Fold Change |
|--------|---------------------|-----------------|-------------------|---------------------------|
| NRXN1 | 4.15x | 3.12x | 4.04x | 1.90x |
| NRXN2 | 1.33x | 1.06x | 1.04x | 1.04x |
| NRXN3 | 1.78x | 1.00x | 1.71x | 1.57x |
| AR-V7 | 1.33x | 1.17x | 1.15x | 1.23x |

---

## Residue 168 (Mutation Site) RNA Contact Analysis

| Target | Condition | Residue 168 Identity | Contacts RNA? | RNA Nucleotides Contacted |
|--------|-----------|---------------------|---------------|--------------------------|
| NRXN1 | WT | ARG | Yes | A18;U16;U17;U7;U8 |
| NRXN1 | R168C | CYS | Yes | A25;A77;G24 |
| NRXN2 | WT | N/A | No | — |
| NRXN2 | R168C | N/A | No | — |
| NRXN3 | WT | N/A | No | — |
| NRXN3 | R168C | N/A | No | — |
| AR-V7 | WT | N/A | No | — |
| AR-V7 | R168C | N/A | No | — |

---

## Overall Interpretation

Across the 4 RNA targets analyzed:

- **2 target(s)** show **increased** viable contacts in R168C (>1.2x WT)
- **0 target(s)** show **decreased** viable contacts in R168C (<0.8x WT)
- **2 target(s)** show **similar** viable contacts in R168C (0.8–1.2x WT)

- **Largest increase**: NRXN1 (4.04x)
- **Largest decrease**: NRXN2 (1.04x)

### Key observations:

- **NRXN1**: R168C shows increased contacts (4.04x WT for pairs)
- **NRXN2**: R168C shows similar contacts (1.04x WT for pairs)
- **NRXN3**: R168C shows increased contacts (1.71x WT for pairs)
- **AR-V7**: R168C shows similar contacts (1.15x WT for pairs)

### Residue 168 (mutation site):

- Residue 168 contacts RNA in **1/4** R168C models
- Residue 168 contacts RNA in **1/4** WT models

The R168C mutation does not substantially change the frequency of direct
RNA contact at residue 168; any functional impact likely arises from broader
structural or electrostatic changes rather than this single residue contact.

---

*Report generated from docking_contact_results/ across 4 targets, 10 models each.*