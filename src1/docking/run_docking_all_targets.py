import os

import sys

import math

import random

import numpy as np

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

PROTEIN_DIRS = {'WT': PROJECT_ROOT / 'results' / 'structural_analysis' / 'WT', 'R168C': PROJECT_ROOT / 'results' / 'structural_analysis' / 'R168C'}

RNA_DIR = PROJECT_ROOT / 'results' / 'target_isoform_docking' / 'target_3d_structs'

RNA_FILES = {'NRXN1': RNA_DIR / 'NRXN1.pdb', 'NRXN2': RNA_DIR / 'NRXN2.pdb', 'NRXN3': RNA_DIR / 'NRXN3.pdb', 'AR_V7': RNA_DIR / 'AR-V7.pdb'}

PROTEIN_FILE = 'KHDRBS2_{}_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_154168.pdb'

OUTPUT_BASE = PROJECT_ROOT

N_ORIENTATIONS = 200

TOP_N = 10

CONTACT_CUTOFF = 8.0

SEED = 42

def read_pdb_atoms(pdb_path):
    atoms = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            residue_number = line[22:26].strip()
            insertion_code = line[26:27].strip()
            chain = line[21:22].strip()
            element = line[76:78].strip() if len(line) >= 78 else ''

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            atoms.append({'name': atom_name, 'resname': residue_name, 'resnum': residue_number, 'icode': insertion_code, 'chain': chain, 'element': element, 'xyz': (x, y, z), 'line': line.rstrip('\n')})
    return atoms

def center_atoms(atoms):
    for a in atoms:
        _converted_1.append(a['xyz'])
    coords = np.array(_converted_1)
    center = coords.mean(axis=0)

    for a in atoms:
        a['xyz'] = tuple(a['xyz'] - center)
    return (atoms, center)

def random_rotation_matrix(rng):
    A = rng.normal(size=(3, 3))
    Q, R = np.linalg.qr(A)

    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q

def rotate_atoms(atoms, R):
    for a in atoms:
        xyz = np.array(a['xyz'])
        a['xyz'] = tuple(R @ xyz)

def translate_atoms(atoms, t):
    for a in atoms:
        a['xyz'] = tuple(np.array(a['xyz']) + t)

def compute_contact_score(protein_coords, rna_coords, cutoff=8.0):
    cutoff_sq = cutoff ** 2
    count = 0
    grid = {}
    cell_size = cutoff

    for i, atom in enumerate(rna_coords):
        key = (int(atom[0] // cell_size), int(atom[1] // cell_size), int(atom[2] // cell_size))

        if key not in grid:
            grid[key] = []
        grid[key].append(i)

    for i in protein_coords:
        p_key = (int(i[0] // cell_size), int(i[1] // cell_size), int(i[2] // cell_size))

        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    key = (p_key[0] + dx, p_key[1] + dy, p_key[2] + dz)

                    if key in grid:
                        for j in grid[key]:
                            diff = i - rna_coords[j]

                            if (diff ** 2).sum() <= cutoff_sq:
                                count += 1
    return count

def dock_protein_rna(protein_atoms, rna_atoms, rng, n_orientations=500):
    protein_atoms, p_center = center_atoms(protein_atoms)
    rna_atoms, r_center = center_atoms(rna_atoms)

    for a in protein_atoms:
        _converted_2.append(a['xyz'])
    p_coords = np.array(_converted_2)

    for a in rna_atoms:
        _converted_3.append(a['xyz'])
    r_coords = np.array(_converted_3)
    p_radius = np.sqrt((p_coords ** 2).sum(axis=1)).max()
    r_radius = np.sqrt((r_coords ** 2).sum(axis=1)).max()
    surface_distance = p_radius * 0.85 + r_radius * 0.5
    candidates = []

    for trial in range(n_orientations):
        R = random_rotation_matrix(rng)
        rot_rna = r_coords @ R.T
        direction = rng.normal(size=3)
        direction /= np.linalg.norm(direction)
        translation = direction * surface_distance
        placed_rna = rot_rna + translation
        score = compute_contact_score(p_coords, placed_rna, CONTACT_CUTOFF)
        candidates.append({'score': score, 'R': R, 't': translation})
    candidates.sort(key=lambda c: c['score'], reverse=True)
    return candidates[:TOP_N]

def write_model_pdb(protein_atoms, rna_atoms, R, t, output_path, score, model_number):
    with open(output_path, 'w') as f:
        f.write(f'REMARK Number:     {model_number}\n')
        f.write(f'REMARK Ligand:     rna.pdb\n')
        f.write(f'REMARK Contact:     0     0     0\n')
        f.write(f'REMARK Score:  {score:.2f}\n')
        f.write(f'REMARK RMSD:     0.00\n')
        f.write('MODEL        1\n')
        atom_index = 1

        for a in protein_atoms:
            x, y, z = a['xyz']
            line = f"ATOM  {atom_index:5d} {a['name']:<4s}{a['resname']:>3s} {a['chain']}{a['resnum']:>4s}{a['icode']:1s}   {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {a['element']:>2s}\n"
            f.write(line)
            atom_index += 1
        f.write('TER\n')

        for a in rna_atoms:
            xyz = np.array(a['xyz'])
            transformed = R @ xyz + t
            x, y, z = transformed
            line = f"ATOM  {atom_index:5d} {a['name']:<4s}{a['resname']:>3s} {a['chain']}{a['resnum']:>4s}{a['icode']:1s}   {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {a['element']:>2s}\n"
            f.write(line)
            atom_index += 1
        f.write('TER\n')
        f.write('ENDMDL\n')
        f.write('END\n')

def main():
    rng = np.random.RandomState(SEED)

    for condition, protein_dir in PROTEIN_DIRS.items():
        protein_path = protein_dir / PROTEIN_FILE.format(condition)

        if not protein_path.exists():
            print(f'WARNING: Protein PDB not found: {protein_path}')
            continue
        protein_atoms = read_pdb_atoms(protein_path)
        print(f"\n{'=' * 60}")
        print(f'Condition: {condition}')
        print(f'  Protein: {protein_path.name} ({len(protein_atoms)} atoms)')

        for rna_name, rna_path in RNA_FILES.items():
            if not rna_path.exists():
                print(f'  WARNING: RNA PDB not found: {rna_path}')
                continue
            rna_atoms = read_pdb_atoms(rna_path)
            print(f'  RNA: {rna_name} ({len(rna_atoms)} atoms)')
            output_dir = OUTPUT_BASE / condition / rna_name
            output_dir.mkdir(parents=True, exist_ok=True)
            existing = list(output_dir.glob('model_*.pdb'))

            if len(existing) >= TOP_N:
                print(f'    Already docked ({len(existing)} models), skipping.')
                continue
            p_all = read_pdb_atoms(protein_path)
            r_all = read_pdb_atoms(rna_path)

            for a in p_all:
                if a['name'] == 'CA':
                    _converted_4.append(a)
            p_score = _converted_4

            for a in r_all:
                if a['name'] in ('P', "C4'", "C1'"):
                    _converted_5.append(a)
            r_score = _converted_5

            for a in p_score:
                _converted_6.append(a['xyz'])
            p_score_coords = np.array(_converted_6)

            for a in r_score:
                _converted_7.append(a['xyz'])
            r_score_coords = np.array(_converted_7)
            p_center = p_score_coords.mean(axis=0)
            r_center = r_score_coords.mean(axis=0)
            p_centered = p_score_coords - p_center
            r_centered = r_score_coords - r_center
            p_radius = np.sqrt((p_centered ** 2).sum(axis=1)).max()
            r_radius = np.sqrt((r_centered ** 2).sum(axis=1)).max()
            surface_distance = p_radius * 0.85 + r_radius * 0.5
            candidates = []

            for trial in range(N_ORIENTATIONS):
                R = random_rotation_matrix(rng)
                rot_rna = r_centered @ R.T
                direction = rng.normal(size=3)
                direction /= np.linalg.norm(direction)
                translation = direction * surface_distance
                placed_rna = rot_rna + translation
                score = compute_contact_score(p_centered, placed_rna, CONTACT_CUTOFF)
                candidates.append({'score': score, 'R': R, 't': translation})
            candidates.sort(key=lambda c: c['score'], reverse=True)
            results = candidates[:TOP_N]
            p_full = read_pdb_atoms(protein_path)
            r_full = read_pdb_atoms(rna_path)

            for a in p_full:
                _converted_8.append(a['xyz'])
            p_full_coords = np.array(_converted_8)

            for a in r_full:
                _converted_9.append(a['xyz'])
            r_full_coords = np.array(_converted_9)
            p_full_center = p_full_coords.mean(axis=0)
            r_full_center = r_full_coords.mean(axis=0)

            for a in p_full:
                a['xyz'] = tuple(np.array(a['xyz']) - p_full_center)

            for a in r_full:
                a['xyz'] = tuple(np.array(a['xyz']) - r_full_center)

            for i, result in enumerate(results, start=1):
                out_path = output_dir / f'model_{i}.pdb'
                write_model_pdb(p_full, r_full, result['R'], result['t'], out_path, result['score'], i)
            print(f'    Docked: top {len(results)} models saved in {output_dir}')

if __name__ == '__main__':
    main()
