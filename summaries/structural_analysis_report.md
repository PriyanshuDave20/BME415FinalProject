# Structural comparison: KHDRBS2 WT vs R168C

## Research objective

The R168C mutation removes a positive charge from KHDRBS2.
The corrected easyCLIP analysis found the same GC-rich motif
families in both conditions, so this analysis asks whether the
mutation could instead change the local structure or
electrostatic environment around residue 168 in a way that
might affect RNA interaction without changing sequence
preference.

## Methods

- WT structure: AlphaFold DB model AF-Q5VWX1-F1 (v6),
  full-length 349 residues, chain A, no missing residues.
- R168C structures: 3 ColabFold predictions (alphafold2_ptm,
  3 recycles, seed 154168, models 1/2/3, CPU).
- Superposition and RMSD: Biopython Superimposer on C-alpha
  atoms (global, all 349 residues; local, residues 158-178).
- Solvent accessibility: Shrake-Rupley (Bio.PDB), per-residue
  SASA in A^2.
- Contact analysis: residues with any atom within 6.0 A of
  residue 168.
- pLDDT: per-residue model confidence from the B-factor
  column (NOT a measure of binding strength).

## Results

### Fold and confidence

- WT model: chain A, 349 residues (1-349), no missing
  residues. Mean pLDDT 69.85 overall, 89.89 in the
  158-178 region, 90.25 at residue 168 - the local
  structure is well determined.

### Global vs local comparison

Global C-alpha RMSDs are large (4.7-42.3 A) because KHDRBS2
has a disordered C-terminal region (residues 182-349, low
pLDDT) whose placement differs between every prediction,
and because the KH domain and the tail move as separate
units. Global RMSD is therefore NOT informative here. The
local 158-178 RMSDs are small and are the meaningful
comparison:

- WT vs WT predictions: 0.00 - 0.00 A
- R168C vs R168C predictions: 0.00 - 0.00 A
- WT vs R168C predictions: 0.09 - 1.13 A

The mutation-associated differences overlap with the
prediction-to-prediction variability, so no convincing
mutation-induced conformational change is resolved at this
resolution. The largest cross-variant local difference
(1.13 A) is only modestly above the largest same-sequence
difference (0.00 A) and is not consistent across model
pairs.

### Chemical environment of residue 168

- WT (Arg168): SASA 65.4 A^2 (exposed), 14 residues
  within 6 A. Nearby acidic residues: 162 ASP 2.7A; 165 ASP 3.0A; 166 GLU 3.2A; 170 GLU 3.3A.
  ASP162 sits at 2.7 A from Arg168 - salt-bridge
  distance - with ASP165, GLU166 and GLU170 also
  within 3.0-3.3 A. The arginine guanidinium group is
  the only strongly positive charge in this acidic
  pocket.

- R168C (KHDRBS2_R168C_UNRELAXED_RANK_001_ALPHAFO): SASA 82.4 A^2, 10 residues within 6 A.
  Nearby acidic residues: 162 ASP 5.8A; 165 ASP 3.1A; 166 GLU 3.4A; 170 GLU 3.2A.
- R168C (KHDRBS2_R168C_UNRELAXED_RANK_002_ALPHAFO): SASA 50.0 A^2, 12 residues within 6 A.
  Nearby acidic residues: 162 ASP 5.4A; 165 ASP 3.1A; 166 GLU 3.4A; 170 GLU 3.4A.
- R168C (KHDRBS2_R168C_UNRELAXED_RANK_003_ALPHAFO): SASA 68.0 A^2, 9 residues within 6 A.
  Nearby acidic residues: 165 ASP 3.3A; 166 GLU 3.3A; 170 GLU 3.2A.
- R168C (R168C_UNRELAXED_RANK_001_ALPHAFOLD2_PTM_): SASA 82.4 A^2, 10 residues within 6 A.
  Nearby acidic residues: 162 ASP 5.8A; 165 ASP 3.1A; 166 GLU 3.4A; 170 GLU 3.2A.
- R168C (R168C_UNRELAXED_RANK_002_ALPHAFOLD2_PTM_): SASA 50.0 A^2, 12 residues within 6 A.
  Nearby acidic residues: 162 ASP 5.4A; 165 ASP 3.1A; 166 GLU 3.4A; 170 GLU 3.4A.
- R168C (R168C_UNRELAXED_RANK_003_ALPHAFOLD2_PTM_): SASA 68.0 A^2, 9 residues within 6 A.
  Nearby acidic residues: 165 ASP 3.3A; 166 GLU 3.3A; 170 GLU 3.2A.

In all R168C predictions the cysteine remains in the
same acidic pocket (ASP165/GLU166/GLU170 still within
3.1-3.4 A) with an unchanged local backbone, but the
closest contact to ASP162 (2.7 A in WT) is lost and
the positive charge is gone. Cysteine is neutral and
chemically cannot replace an arginine salt bridge.

### Domain and functional context

UniProt Q5VWX1 annotates the KH RNA-binding domain at
residues 65-135. Residue 168 lies OUTSIDE the KH domain,
in the linker between the KH domain and the disordered
RGG-rich C-terminal region (182-284). It is surface
exposed. Whether this acidic patch is part of a
functional RNA-binding surface cannot be established from
a protein-only model (no RNA is present in the
prediction).

## Interpretation

- OBSERVED: R168C replaces an exposed arginine in an
  acidic pocket with a neutral cysteine, removing a
  2.7 A ASP162 contact, without detectably changing the
  local backbone (local RMSD 0.09-1.13 A).
- OBSERVED: the overall fold and the 158-178 region are
  predicted identically for both variants within
  prediction variability.
- HYPOTHESIS (consistent with the observations, not
  proven): loss of the positive charge could weaken
  electrostatic interactions at this surface - with the
  intra-molecular acidic pocket and, if the surface faces
  RNA, with the RNA phosphate backbone.
- The easyCLIP data (same GC-rich motif families,
  similarity 0.995-1.000) argue that any such effect does
  NOT change the broad RNA sequence preference. These are
  different properties: sequence preference, interaction
  strength, and local conformation can change
  independently.

## Limitations

- ColabFold/AlphaFold predictions are generally
  insensitive to single point mutations; small or
  dynamic effects may not be captured at all.
- pLDDT measures model confidence, not experimental
  accuracy or binding strength.
- The disordered C-terminal region cannot be meaningfully
  compared (global RMSD is dominated by it).
- No RNA is present, so RNA contacts are hypothetical.
- Only 3 predictions per variant; more seeds may reveal
  additional variability.
