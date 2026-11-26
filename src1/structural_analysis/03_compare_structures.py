import csv

import glob

import os

import sys

from Bio.PDB import MMCIFParser, PDBParser, Superimposer

from Bio.PDB.Polypeptide import protein_letters_3to1

local_start, local_end = (158, 178)

position_of_interest = 168

out_folder = 'results/structural_analysis'

def load_structures():
    wt = []
    mutant = []
    parser = MMCIFParser(QUIET=True)
    pdb_parser = PDBParser(QUIET=True)

    for path in glob.glob(os.path.join(out_folder, '*.cif')) + glob.glob(os.path.join(out_folder, '*.pdb')) + glob.glob(os.path.join(out_folder, 'WT', '*.cif')) + glob.glob(os.path.join(out_folder, 'WT', '*.pdb')) + glob.glob(os.path.join(out_folder, 'R168C', '*.cif')) + glob.glob(os.path.join(out_folder, 'R168C', '*.pdb')):
        name = os.path.basename(path).upper()

        if 'WT' in name:
            target = wt

        elif 'R168C' in name:
            target = mutant
        else:
            continue

        if path.endswith('.cif'):
            structure = parser.get_structure('s', path)
        else:
            structure = pdb_parser.get_structure('s', path)
        ca_atoms = {}

        for j in structure.get_residues():
            if 'CA' in j:
                ca_atoms[j.id[1]] = j['CA']
        target.append({'path': path, 'ca': ca_atoms, 'n_residues': len(ca_atoms)})
    return (wt, mutant)

def ca_rmsd(atoms_a, atoms_b, residue_range=None):
    common = sorted(set(atoms_a) & set(atoms_b))

    if residue_range:
        for n in common:
            if residue_range[0] <= n <= residue_range[1]:
                _converted_1.append(n)
        common = _converted_1

    for n in common:
        _converted_2.append(atoms_a[n])
    fixed = _converted_2

    for n in common:
        _converted_3.append(atoms_b[n])
    moving = _converted_3
    superimposer = Superimposer()
    superimposer.set_atoms(fixed, moving)
    return (superimposer.rms, len(common), superimposer)

def main():
    wt, mutant = load_structures()
    print('WT structures: %d, R168C structures: %d' % (len(wt), len(mutant)))

    if not wt or not mutant:
        print('ERROR: this comparison needs at least one WT and one R168C')
        print('structure file in', out_folder)
        print('The WT AlphaFold model is already present.')
        print('The R168C predictions must be generated with ColabFold using')
        print('results/structural_analysis/r168c_sequences.fasta')
        print('(see results/structural_analysis/colabfold_instructions.md)')
        sys.exit(1)
    rows = []

    for i in wt:
        for m in mutant:
            glob_rms, n_all, _ = ca_rmsd(i['ca'], m['ca'])
            loc_rms, n_loc, sup_local = ca_rmsd(i['ca'], m['ca'], (local_start, local_end))
            displacements = {}
            common = sorted(set(i['ca']) & set(m['ca']))

            for n in common:
                a = i['ca'][n]
                b = m['ca'][n]
                dist = 0.0

                for j in range(3):
                    diff = a.coord[j] - b.coord[j]
                    dist += diff * diff
                displacements[n] = dist ** 0.5
            rows.append({'wt': i['path'], 'mut': m['path'], 'global_rmsd': glob_rms, 'n_all': n_all, 'local_rmsd': loc_rms, 'n_local': n_loc, 'displacements': displacements})
            print('%s vs %s: global RMSD %.2f A (%d CA), local RMSD %d-%d %.2f A (%d CA)' % (os.path.basename(i['path']), os.path.basename(m['path']), glob_rms, n_all, local_start, local_end, loc_rms, n_loc))
    path = os.path.join(out_folder, 'per_residue_displacement.tsv')

    with open(path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')

        for i in range(1, len(rows) + 1):
            _converted_4.append('displacement_%d' % i)
        writer.writerow(['residue_number'] + _converted_4)
        all_numbers = sorted(rows[0]['displacements'])

        for n in all_numbers:
            for r in rows:
                _converted_5.append('%.2f' % r['displacements'].get(n, ''))
            writer.writerow([n] + _converted_5)
    print('  wrote', path)
    path = os.path.join(out_folder, 'rmsd_comparison.tsv')

    with open(path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['wt_file', 'r168c_file', 'global_rmsd_A', 'n_CA_superimposed', 'local_rmsd_158_178_A', 'n_CA_local'])

        for r in rows:
            writer.writerow([os.path.basename(r['wt']), os.path.basename(r['mut']), '%.2f' % r['global_rmsd'], r['n_all'], '%.2f' % r['local_rmsd'], r['n_local']])
    print('  wrote', path)
    print('\nNote: these are predictions of the same protein with one point')
    print('mutation. ColabFold is generally insensitive to single point')
    print('mutations, so small differences are expected and are only')
    print('meaningful when compared with the variability between seeds.')

if __name__ == '__main__':
    main()
