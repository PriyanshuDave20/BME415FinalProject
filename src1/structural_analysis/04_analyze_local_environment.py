import csv

import glob

import os

import statistics

import sys

from Bio.PDB import MMCIFParser, PDBParser, ShrakeRupley

from Bio.PDB.Polypeptide import protein_letters_3to1

out_folder = 'results/structural_analysis'

position_of_interest = 168

local_start, local_end = (158, 178)

contact_distance = 6.0

hbond_distance = 3.5

charged_positive = ['ARG', 'LYS', 'HIS']

charged_negative = ['ASP', 'GLU']

polar = ['ASN', 'GLN', 'SER', 'THR', 'TYR', 'CYS', 'TRP', 'HIS']

hydrophobic = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'PRO', 'GLY']

def load_structures():
    structures = []
    parser = MMCIFParser(QUIET=True)
    pdb_parser = PDBParser(QUIET=True)

    for path in glob.glob(os.path.join(out_folder, '*.cif')) + glob.glob(os.path.join(out_folder, '*.pdb')) + glob.glob(os.path.join(out_folder, 'WT', '*.cif')) + glob.glob(os.path.join(out_folder, 'WT', '*.pdb')) + glob.glob(os.path.join(out_folder, 'R168C', '*.cif')) + glob.glob(os.path.join(out_folder, 'R168C', '*.pdb')):
        name = os.path.basename(path).upper()

        if 'WT' not in name and 'R168C' not in name:
            continue

        if path.endswith('.cif'):
            structure = parser.get_structure('s', path)
        else:
            structure = pdb_parser.get_structure('s', path)
        structures.append((name, structure))
    return structures

def analyze_environment(structure, name):
    chain = None
    target = None

    for i in structure:
        for c in i:
            for r in c:
                if r.id[0] == ' ':
                    _converted_1.append(r)
            residues = _converted_1

            for r in residues:
                _converted_2.append(r.id[1])
            numbers = _converted_2

            if position_of_interest in numbers:
                chain = c

                for r in residues:
                    if r.id[1] == position_of_interest:
                        _converted_3.append(r)
                target = _converted_3[0]
                break

        if target is not None:
            break

    if target is None:
        return None
    sasa = ShrakeRupley()
    sasa.compute(chain, level='R')

    for r in chain:
        if r.id[0] == ' ':
            _converted_4.append(r)
    residues = _converted_4

    for r in residues:
        _converted_5[r.id[1]] = r
    by_number = _converted_5
    nearby = []

    for r in residues:
        if r.id[1] == position_of_interest:
            continue
        min_dist = None

        for atom_a in target:
            for j in r:
                dist = atom_a - j

                if min_dist is None or dist < min_dist:
                    min_dist = dist

        if min_dist is not None and min_dist <= contact_distance:
            nearby.append((r.id[1], r.resname, min_dist))

    for n, rn, d in nearby:
        if rn in charged_positive:
            _converted_6.append((n, rn, round(d, 1)))

    for n, rn, d in nearby:
        if rn in charged_negative:
            _converted_7.append((n, rn, round(d, 1)))

    for n, rn, d in nearby:
        if rn in polar:
            _converted_8.append((n, rn, round(d, 1)))

    for n, rn, d in nearby:
        if rn in hydrophobic:
            _converted_9.append((n, rn, round(d, 1)))
    result = {'name': name, 'resname_168': target.resname, 'sasa_168': round(target.sasa, 1), 'sasa_local': {}, 'nearby_positive': _converted_6, 'nearby_negative': _converted_7, 'nearby_polar': _converted_8, 'nearby_hydrophobic': _converted_9, 'n_contacts': len(nearby)}

    for n in range(local_start, local_end + 1):
        if n in by_number and 'sasa' in by_number[n].xtra:
            result['sasa_local'][n] = round(by_number[n].sasa, 1)
    return result

def main():
    structures = load_structures()

    for n, s in structures:
        if 'WT' in n.upper():
            _converted_10.append((n, s))
    wt = _converted_10

    for n, s in structures:
        if 'R168C' in n.upper():
            _converted_11.append((n, s))
    mut = _converted_11

    if not wt:
        print('ERROR: no WT structure found in', out_folder)
        sys.exit(1)
    results = []

    for name, structure in wt + mut:
        print('Analyzing', name, '...')
        result = analyze_environment(structure, name)

        if result is None:
            print('  ERROR: residue %d not found in this structure' % position_of_interest)
            sys.exit(1)
        results.append(result)
        print('  residue %d is %s, SASA %.1f A^2, %d residues within %.1f A' % (position_of_interest, result['resname_168'], result['sasa_168'], result['n_contacts'], contact_distance))

        for j in ('nearby_positive', 'nearby_negative', 'nearby_polar', 'nearby_hydrophobic'):
            print('    %s: %s' % (j, result[j]))
    comparison_pending = len(mut) == 0

    if comparison_pending:
        print('\nWARNING: no R168C structure available yet - the WT-vs-mutant')
        print('environment comparison will be done after ColabFold predictions')
        print('are placed in', out_folder)
        print('Submit results/structural_analysis/r168c_sequences.fasta to')
        print('ColabFold (see colabfold_instructions.md in the same folder)')
    path = os.path.join(out_folder, 'local_environment_comparison.tsv')

    with open(path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['structure', 'resname_168', 'sasa_168', 'n_residues_within_6A', 'nearby_positive', 'nearby_negative', 'nearby_polar', 'nearby_hydrophobic'])

        for r in results:
            for c in r['nearby_positive']:
                _converted_12.append('%d %s %.1fA' % c)

            for c in r['nearby_negative']:
                _converted_13.append('%d %s %.1fA' % c)

            for c in r['nearby_polar']:
                _converted_14.append('%d %s %.1fA' % c)

            for c in r['nearby_hydrophobic']:
                _converted_15.append('%d %s %.1fA' % c)
            writer.writerow([r['name'], r['resname_168'], r['sasa_168'], r['n_contacts'], '; '.join(_converted_12), '; '.join(_converted_13), '; '.join(_converted_14), '; '.join(_converted_15)])
    print('\n  wrote', path)

    if comparison_pending:
        print('\nSTATUS: WT analysis complete, R168C comparison pending')
        print('ColabFold predictions')

if __name__ == '__main__':
    main()
