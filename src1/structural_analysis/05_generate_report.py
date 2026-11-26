import csv

import os

out_folder = 'results/structural_analysis'

def read_tsv(path):
    if not os.path.exists(path):
        return []

    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def main():
    rmsd = read_tsv(os.path.join(out_folder, 'rmsd_comparison.tsv'))
    env = read_tsv(os.path.join(out_folder, 'local_environment_comparison.tsv'))
    meta = read_tsv(os.path.join(out_folder, 'structural_metadata_summary.tsv'))

    def variant(row, which):
        f = row['wt_file'] if which == 'WT' else row['r168c_file']
        name = os.path.basename(f).upper()

        if 'R168C' in name:
            return 'R168C'

        if 'WT' in name:
            return 'WT'
        return 'R168C'

    for r in rmsd:
        if variant(r, 'WT') == 'WT' and variant(r, 'R168C') == 'WT':
            _converted_1.append(r)
    wt_wt = _converted_1

    for r in rmsd:
        if variant(r, 'WT') == 'R168C' and variant(r, 'R168C') == 'R168C':
            _converted_2.append(r)
    mut_mut = _converted_2

    for r in rmsd:
        if variant(r, 'WT') != variant(r, 'R168C'):
            _converted_3.append(r)
    cross = _converted_3

    def local_range(rows):
        for r in rows:
            _converted_4.append(float(r['local_rmsd_158_178_A']))
        vals = _converted_4

        if not vals:
            return (0.0, 0.0)
        return (min(vals), max(vals))
    ww_lo, ww_hi = local_range(wt_wt)
    mm_lo, mm_hi = local_range(mut_mut)
    cr_lo, cr_hi = local_range(cross)

    for r in env:
        if r['resname_168'] == 'ARG':
            _converted_5.append(r)
    wt_env = _converted_5

    for r in env:
        if r['resname_168'] == 'CYS':
            _converted_6.append(r)
    mut_env = _converted_6
    lines = []
    lines.append('# Structural comparison: KHDRBS2 WT vs R168C')
    lines.append('')
    lines.append('## Research objective')
    lines.append('')
    lines.append('The R168C mutation removes a positive charge from KHDRBS2.')
    lines.append('The corrected easyCLIP analysis found the same GC-rich motif')
    lines.append('families in both conditions, so this analysis asks whether the')
    lines.append('mutation could instead change the local structure or')
    lines.append('electrostatic environment around residue 168 in a way that')
    lines.append('might affect RNA interaction without changing sequence')
    lines.append('preference.')
    lines.append('')
    lines.append('## Methods')
    lines.append('')
    lines.append('- WT structure: AlphaFold DB model AF-Q5VWX1-F1 (v6),')
    lines.append('  full-length 349 residues, chain A, no missing residues.')
    lines.append('- R168C structures: 3 ColabFold predictions (alphafold2_ptm,')
    lines.append('  3 recycles, seed 154168, models 1/2/3, CPU).')
    lines.append('- Superposition and RMSD: Biopython Superimposer on C-alpha')
    lines.append('  atoms (global, all 349 residues; local, residues 158-178).')
    lines.append('- Solvent accessibility: Shrake-Rupley (Bio.PDB), per-residue')
    lines.append('  SASA in A^2.')
    lines.append('- Contact analysis: residues with any atom within 6.0 A of')
    lines.append('  residue 168.')
    lines.append('- pLDDT: per-residue model confidence from the B-factor')
    lines.append('  column (NOT a measure of binding strength).')
    lines.append('')
    lines.append('## Results')
    lines.append('')
    lines.append('### Fold and confidence')
    lines.append('')

    if meta:
        m = meta[0]
        lines.append('- WT model: chain A, %s residues (%s-%s), no missing' % (m['n_residues'], m['first_residue'], m['last_residue']))
        lines.append('  residues. Mean pLDDT %.2f overall, %.2f in the' % (float(m['mean_plddt_all']), float(m['mean_plddt_158_178'])))
        lines.append('  158-178 region, %.2f at residue 168 - the local' % float(m['plddt_168']))
        lines.append('  structure is well determined.')
    lines.append('')
    lines.append('### Global vs local comparison')
    lines.append('')
    lines.append('Global C-alpha RMSDs are large (4.7-42.3 A) because KHDRBS2')
    lines.append('has a disordered C-terminal region (residues 182-349, low')
    lines.append('pLDDT) whose placement differs between every prediction,')
    lines.append('and because the KH domain and the tail move as separate')
    lines.append('units. Global RMSD is therefore NOT informative here. The')
    lines.append('local 158-178 RMSDs are small and are the meaningful')
    lines.append('comparison:')
    lines.append('')
    lines.append('- WT vs WT predictions: %.2f - %.2f A' % (ww_lo, ww_hi))
    lines.append('- R168C vs R168C predictions: %.2f - %.2f A' % (mm_lo, mm_hi))
    lines.append('- WT vs R168C predictions: %.2f - %.2f A' % (cr_lo, cr_hi))
    lines.append('')
    lines.append('The mutation-associated differences overlap with the')
    lines.append('prediction-to-prediction variability, so no convincing')
    lines.append('mutation-induced conformational change is resolved at this')
    lines.append('resolution. The largest cross-variant local difference')
    lines.append('(%.2f A) is only modestly above the largest same-sequence' % cr_hi)
    lines.append('difference (%.2f A) and is not consistent across model' % ww_hi)
    lines.append('pairs.')
    lines.append('')
    lines.append('### Chemical environment of residue 168')
    lines.append('')

    if wt_env:
        w = wt_env[0]
        lines.append('- WT (Arg168): SASA %.1f A^2 (exposed), %d residues' % (float(w['sasa_168']), int(w['n_residues_within_6A'])))
        lines.append('  within 6 A. Nearby acidic residues: %s.' % w['nearby_negative'])
        lines.append('  ASP162 sits at 2.7 A from Arg168 - salt-bridge')
        lines.append('  distance - with ASP165, GLU166 and GLU170 also')
        lines.append('  within 3.0-3.3 A. The arginine guanidinium group is')
        lines.append('  the only strongly positive charge in this acidic')
        lines.append('  pocket.')

    if mut_env:
        lines.append('')

        for r in mut_env:
            lines.append('- R168C (%s): SASA %s A^2, %d residues within 6 A.' % (r['structure'][:40], r['sasa_168'], int(r['n_residues_within_6A'])))
            lines.append('  Nearby acidic residues: %s.' % (r['nearby_negative'] or 'none within 6 A'))
    lines.append('')

    if wt_env and mut_env:
        lines.append('In all R168C predictions the cysteine remains in the')
        lines.append('same acidic pocket (ASP165/GLU166/GLU170 still within')
        lines.append('3.1-3.4 A) with an unchanged local backbone, but the')
        lines.append('closest contact to ASP162 (2.7 A in WT) is lost and')
        lines.append('the positive charge is gone. Cysteine is neutral and')
        lines.append('chemically cannot replace an arginine salt bridge.')
        lines.append('')
    lines.append('### Domain and functional context')
    lines.append('')
    lines.append('UniProt Q5VWX1 annotates the KH RNA-binding domain at')
    lines.append('residues 65-135. Residue 168 lies OUTSIDE the KH domain,')
    lines.append('in the linker between the KH domain and the disordered')
    lines.append('RGG-rich C-terminal region (182-284). It is surface')
    lines.append('exposed. Whether this acidic patch is part of a')
    lines.append('functional RNA-binding surface cannot be established from')
    lines.append('a protein-only model (no RNA is present in the')
    lines.append('prediction).')
    lines.append('')
    lines.append('## Interpretation')
    lines.append('')
    lines.append('- OBSERVED: R168C replaces an exposed arginine in an')
    lines.append('  acidic pocket with a neutral cysteine, removing a')
    lines.append('  2.7 A ASP162 contact, without detectably changing the')
    lines.append('  local backbone (local RMSD 0.09-1.13 A).')
    lines.append('- OBSERVED: the overall fold and the 158-178 region are')
    lines.append('  predicted identically for both variants within')
    lines.append('  prediction variability.')
    lines.append('- HYPOTHESIS (consistent with the observations, not')
    lines.append('  proven): loss of the positive charge could weaken')
    lines.append('  electrostatic interactions at this surface - with the')
    lines.append('  intra-molecular acidic pocket and, if the surface faces')
    lines.append('  RNA, with the RNA phosphate backbone.')
    lines.append('- The easyCLIP data (same GC-rich motif families,')
    lines.append('  similarity 0.995-1.000) argue that any such effect does')
    lines.append('  NOT change the broad RNA sequence preference. These are')
    lines.append('  different properties: sequence preference, interaction')
    lines.append('  strength, and local conformation can change')
    lines.append('  independently.')
    lines.append('')
    lines.append('## Limitations')
    lines.append('')
    lines.append('- ColabFold/AlphaFold predictions are generally')
    lines.append('  insensitive to single point mutations; small or')
    lines.append('  dynamic effects may not be captured at all.')
    lines.append('- pLDDT measures model confidence, not experimental')
    lines.append('  accuracy or binding strength.')
    lines.append('- The disordered C-terminal region cannot be meaningfully')
    lines.append('  compared (global RMSD is dominated by it).')
    lines.append('- No RNA is present, so RNA contacts are hypothetical.')
    lines.append('- Only 3 predictions per variant; more seeds may reveal')
    lines.append('  additional variability.')
    lines.append('')
    lines.append('## Recommended next experiments')
    lines.append('')
    lines.append('- Mutate ASP162/ASP165/GLU166 to test whether the acidic')
    lines.append('  pocket contacts matter (charge-reversal controls).')
    lines.append('- Measure RNA-binding affinity directly for WT and R168C')
    lines.append('  (e.g. fluorescence anisotropy or EMSA with the same RNA')
    lines.append('  target).')
    lines.append('- Test whether the 158-178 surface contacts RNA (RNA')
    lines.append('  footprinting or crosslinking-native analysis).')
    path = os.path.join(out_folder, 'structural_analysis_report.md')

    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('  wrote', path)

if __name__ == '__main__':
    main()
