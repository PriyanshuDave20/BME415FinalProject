import sys

from pathlib import Path

from collections import defaultdict

import csv

import re

import numpy as np

from scipy.spatial import cKDTree

CONDITIONS = {'WT': Path('wt/NRXN1'), 'R168C': Path('R168C/NRXN1')}

OUTPUT_DIR = Path('docking_contact_results')

CUTOFFS = [4.0, 5.0]

def is_hydrogen(atom_name, element):
    if element:
        return element.upper() in {'H', 'D'}
    cleaned = atom_name.strip().lstrip('0123456789')
    return cleaned.upper().startswith(('H', 'D'))

def read_complex(pdb_path):
    sections = defaultdict(list)
    section_number = 1

    with open(pdb_path, 'r') as file:
        for i in file:
            if i.startswith('TER'):
                section_number += 1
                continue

            if not i.startswith('ATOM'):
                continue

            if section_number not in (1, 2):
                continue
            atom_name = i[12:16].strip()
            residue_name = i[17:20].strip()
            residue_number = i[22:26].strip()
            insertion_code = i[26:27].strip()
            element = i[76:78].strip() if len(i) >= 78 else ''

            if is_hydrogen(atom_name, element):
                continue

            try:
                x = float(i[30:38])
                y = float(i[38:46])
                z = float(i[46:54])
            except ValueError:
                continue
            residue_key = (residue_name, residue_number, insertion_code)
            sections[section_number].append({'atom': atom_name, 'resname': residue_name, 'resnum': residue_number, 'icode': insertion_code, 'residue_key': residue_key, 'xyz': (x, y, z)})
    protein = sections[1]
    rna = sections[2]

    if not protein or not rna:
        print("Error: ValueError(f'Could not identify both protein and RNA in {pdb_path}. Expected protein before the first TER and RNA after it.')")
        sys.exit(1)
    return (protein, rna)

def analyze_contacts(protein, rna, cutoff):
    for atom in protein:
        _converted_1.append(atom['xyz'])
    protein_coords = np.array(_converted_1, dtype=float)

    for atom in rna:
        _converted_2.append(atom['xyz'])
    rna_coords = np.array(_converted_2, dtype=float)
    protein_tree = cKDTree(protein_coords)
    rna_tree = cKDTree(rna_coords)
    neighbors = protein_tree.query_ball_tree(rna_tree, cutoff)
    atom_pairs = []

    for protein_index, rna_indices in enumerate(neighbors):
        for j in rna_indices:
            atom_pairs.append((protein_index, j))

    if atom_pairs:
        for pair in atom_pairs:
            _converted_3.append(pair[0])
        protein_indices = np.array(_converted_3, dtype=int)

        for pair in atom_pairs:
            _converted_4.append(pair[1])
        rna_indices = np.array(_converted_4, dtype=int)
        distances = np.linalg.norm(protein_coords[protein_indices] - rna_coords[rna_indices], axis=1)
        valid = distances <= cutoff
        protein_indices = protein_indices[valid]
        rna_indices = rna_indices[valid]
        distances = distances[valid]
    else:
        protein_indices = np.array([], dtype=int)
        rna_indices = np.array([], dtype=int)
        distances = np.array([], dtype=float)
    residue_pairs = set()
    protein_residues = set()
    rna_residues = set()
    atom_pair_counts = defaultdict(int)
    protein_to_rna = defaultdict(set)
    protein_min_distance = {}

    for p_index, r_index, distance in zip(protein_indices, rna_indices, distances):
        p_atom = protein[p_index]
        r_atom = rna[r_index]
        p_key = p_atom['residue_key']
        r_key = r_atom['residue_key']
        residue_pairs.add((p_key, r_key))
        protein_residues.add(p_key)
        rna_residues.add(r_key)
        atom_pair_counts[p_key, r_key] += 1
        protein_to_rna[p_key].add(r_key)

        if p_key not in protein_min_distance or distance < protein_min_distance[p_key]:
            protein_min_distance[p_key] = float(distance)

    for atom in protein:
        if atom['resnum'] == '168':
            _converted_5.append(atom)
    residue168 = _converted_5

    if residue168:
        residue168_name = residue168[0]['resname']
        residue168_key = residue168[0]['residue_key']
        residue168_rna = protein_to_rna.get(residue168_key, set())
        residue168_contact = bool(residue168_rna)
    else:
        residue168_name = 'NOT_FOUND'
        residue168_key = None
        residue168_rna = set()
        residue168_contact = False
    return {'atom_pairs': len(protein_indices), 'protein_residues': protein_residues, 'rna_residues': rna_residues, 'residue_pairs': residue_pairs, 'protein_to_rna': protein_to_rna, 'atom_pair_counts': atom_pair_counts, 'protein_min_distance': protein_min_distance, 'residue168_name': residue168_name, 'residue168_contact': residue168_contact, 'residue168_rna': residue168_rna}

def format_residue(residue_key):
    name, number, insertion = residue_key
    return f'{name}{number}{insertion}'

def write_csv(filename, rows, fieldnames):
    with open(filename, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    summary_rows = []
    detail_rows = []

    for condition, folder in CONDITIONS.items():
        pdb_files = sorted(folder.glob('model_*.pdb'), key=lambda path: int(re.search('model_(\\d+)', path.stem).group(1)))

        if not pdb_files:
            print(f'WARNING: No PDB files found in {folder}')
            continue
        print(f'\nAnalyzing {condition}: {len(pdb_files)} structures')

        for pdb_path in pdb_files:
            model_number = int(re.search('model_(\\d+)', pdb_path.stem).group(1))
            protein, rna = read_complex(pdb_path)

            for atom in protein:
                _converted_6.append(atom['residue_key'])
            protein_residue_count = len(_converted_6)

            for atom in rna:
                _converted_7.append(atom['residue_key'])
            rna_residue_count = len(_converted_7)
            print(f'  {pdb_path.name}: {len(protein)} protein heavy atoms, {len(rna)} RNA heavy atoms')
            results = {}

            for cutoff in CUTOFFS:
                results[cutoff] = analyze_contacts(protein, rna, cutoff)

            for cutoff, result in results.items():
                for key in result['residue168_rna']:
                    _converted_8.append(format_residue(key))
                residue168_rna = sorted(_converted_8)
                summary_rows.append({'condition': condition, 'model': model_number, 'cutoff_A': cutoff, 'protein_residues_in_structure': protein_residue_count, 'rna_nucleotides_in_structure': rna_residue_count, 'atom_contact_pairs': result['atom_pairs'], 'unique_protein_residues_contacting_RNA': len(result['protein_residues']), 'unique_RNA_nucleotides_contacting_protein': len(result['rna_residues']), 'unique_residue_residue_pairs': len(result['residue_pairs']), 'residue168_identity': result['residue168_name'], 'residue168_contacts_RNA': result['residue168_contact'], 'residue168_RNA_nucleotides': ';'.join(residue168_rna)})

                for p_key in sorted(result['protein_to_rna'], key=lambda key: int(key[1])):
                    contacted_rna = sorted(result['protein_to_rna'][p_key], key=lambda key: int(key[1]))

                    for r_key in contacted_rna:
                        _converted_9.append(result['atom_pair_counts'][p_key, r_key])
                    atom_pair_count = sum(_converted_9)

                    for key in contacted_rna:
                        _converted_10.append(format_residue(key))
                    detail_rows.append({'condition': condition, 'model': model_number, 'cutoff_A': cutoff, 'protein_residue': format_residue(p_key), 'protein_residue_number': p_key[1], 'protein_residue_name': p_key[0], 'number_of_RNA_nucleotides_contacted': len(contacted_rna), 'RNA_nucleotides_contacted': ';'.join(_converted_10), 'atom_contact_pairs': atom_pair_count, 'minimum_distance_A': round(result['protein_min_distance'][p_key], 3)})
    write_csv(OUTPUT_DIR / 'contact_summary.csv', summary_rows, ['condition', 'model', 'cutoff_A', 'protein_residues_in_structure', 'rna_nucleotides_in_structure', 'atom_contact_pairs', 'unique_protein_residues_contacting_RNA', 'unique_RNA_nucleotides_contacting_protein', 'unique_residue_residue_pairs', 'residue168_identity', 'residue168_contacts_RNA', 'residue168_RNA_nucleotides'])
    write_csv(OUTPUT_DIR / 'protein_residue_contacts.csv', detail_rows, ['condition', 'model', 'cutoff_A', 'protein_residue', 'protein_residue_number', 'protein_residue_name', 'number_of_RNA_nucleotides_contacted', 'RNA_nucleotides_contacted', 'atom_contact_pairs', 'minimum_distance_A'])
    print('\nAnalysis complete.')
    print(f'Results saved in: {OUTPUT_DIR.resolve()}')
    print('  contact_summary.csv')
    print('  protein_residue_contacts.csv')

if __name__ == '__main__':
    main()
