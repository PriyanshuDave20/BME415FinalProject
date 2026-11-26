import csv

import glob

import os

import statistics

import sys

from Bio.PDB import MMCIFParser, PDBParser, ShrakeRupley

out_folder = 'results/structural_analysis'

local_start = 158

local_end = 178

position_of_interest = 168

reported_plddt = {156: 95.12, 157: 94.19, 158: 94.81, 159: 94.69, 160: 92.44, 161: 89.5, 162: 86.5, 163: 76.12, 164: 83.69, 165: 88.25, 166: 91.44, 167: 89.81, 168: 90.25, 169: 92.12, 170: 92.94}

def find_structure_files():
    found = []

    for i in ('*.cif', '*.pdb'):
        for path in glob.glob(i) + glob.glob(os.path.join('results', 'structural_analysis', i)):
            found.append(path)
    return found

def inspect_structure(path):
    if path.endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure('model', path)
    model = structure[0]
    chains = {}

    for chain in model:
        for r in chain:
            if r.id[0] == ' ':
                _converted_1.append(r)
        residues = _converted_1

        for r in residues:
            _converted_2.append(r.id[1])
        numbers = sorted(_converted_2)
        plddt = {}

        for r in residues:
            if 'CA' in r:
                plddt[r.id[1]] = r['CA'].bfactor
        missing = []

        for a, b in zip(numbers, numbers[1:]):
            if b - a > 1:
                missing.append('%d-%d' % (a, b - 1))
        chains[chain.id] = {'n_residues': len(residues), 'first': numbers[0], 'last': numbers[-1], 'missing': missing, 'plddt': plddt, 'residues': residues}
    return {'path': path, 'chains': chains}

UNIPROT_FEATURES = [('KH domain', 65, 135), ('Disordered region', 182, 284), ('Disordered region', 319, 349)]

FOCUS_REGIONS = [('KH domain (65-135)', 65, 135), ('KH C-terminal part (120-145)', 120, 145), ('Linker before disordered region (136-181)', 136, 181), ('Linker around R168 (150-181)', 150, 181), ('Disordered region start (182-200)', 182, 200), ('Annotated disordered region (182-284)', 182, 284), ('Annotated disordered region (319-349)', 319, 349)]

def region_context(chain):
    for r in chain:
        if r.id[0] == ' ':
            _converted_3.append(r)
    residues = _converted_3

    for r in residues:
        _converted_4[r.id[1]] = r
    by_number = _converted_4
    ShrakeRupley().compute(chain, level='R')
    region_rows = []

    for name, start, end in FOCUS_REGIONS + UNIPROT_FEATURES:
        for n in range(start, end + 1):
            if n in by_number:
                _converted_5.append(by_number[n]['CA'].bfactor)
        vals = _converted_5

        for n in range(start, end + 1):
            if n in by_number:
                _converted_6.append(by_number[n].sasa)
        sasa_vals = _converted_6

        if not vals:
            continue

        for s in sasa_vals:
            if s > 30:
                _converted_7.append(1)
        exposed = sum(_converted_7)
        region_rows.append([name, start, end, len(vals), round(statistics.mean(vals), 2), round(min(vals), 2), round(statistics.mean(sasa_vals), 1), round(100.0 * exposed / len(sasa_vals), 1)])
    return (region_rows, by_number, residues)

def linker_kh_contacts(chain):
    for r in chain:
        if r.id[0] == ' ':
            _converted_8.append(r)
    residues = _converted_8

    for r in residues:
        _converted_9[r.id[1]] = r
    by_number = _converted_9

    for r in residues:
        if 136 <= r.id[1] <= 181:
            _converted_10.append(r)
    linker = _converted_10

    for r in residues:
        if 65 <= r.id[1] <= 135:
            _converted_11.append(r)
    kh = _converted_11
    contacting = []

    for r in linker:
        touching = set()

        for atom_a in r:
            for k in kh:
                for atom_b in k:
                    if atom_a - atom_b <= 6.0:
                        touching.add(k.id[1])

        if touching:
            contacting.append((r.id[1], r.resname, min(touching), max(touching), len(touching)))
    r168 = by_number[position_of_interest]
    closest = None

    for k in kh:
        for atom_a in r168:
            for atom_b in k:
                _converted_12.append(atom_a - atom_b)
        dist = min(_converted_12)

        if closest is None or dist < closest[2]:
            closest = (k.id[1], k.resname, dist)
    return {'linker_n': len(linker), 'linker_contacting_kh': contacting, 'closest_kh_to_168': closest, 'plddt_168': by_number[position_of_interest]['CA'].bfactor}

def main():
    os.makedirs(out_folder, exist_ok=True)
    files = find_structure_files()

    if not files:
        print('ERROR: no structure files (.cif or .pdb) found in the project')
        sys.exit(1)
    print('Found %d structure file(s)' % len(files))
    summary_rows = []
    context_rows = []

    for path in files:
        variant = 'R168C' if 'R168C' in os.path.basename(path).upper() else 'WT'
        print('\nStructure:', path, '(treated as %s)' % variant)
        meta = inspect_structure(path)

        for chain_id, c in meta['chains'].items():
            print('  chain %s: %d residues, numbering %d-%d' % (chain_id, c['n_residues'], c['first'], c['last']))

            if c['missing']:
                print('  missing residues: %s' % ', '.join(c['missing']))
            else:
                print('  missing residues: none')
            plddt = c['plddt']

            for n in range(local_start, local_end + 1):
                if n in plddt:
                    _converted_13.append(plddt[n])
            local_vals = _converted_13
            print('  mean pLDDT (all residues): %.2f' % statistics.mean(plddt.values()))
            print('  mean pLDDT (residues %d-%d): %.2f' % (local_start, local_end, statistics.mean(local_vals)))
            print('  pLDDT at residue %d: %.2f' % (position_of_interest, plddt.get(position_of_interest, float('nan'))))
            mismatches = []

            for n in sorted(reported_plddt):
                model_value = round(plddt.get(n, -1), 2)

                if abs(model_value - reported_plddt[n]) > 0.01:
                    mismatches.append('%d: model %.2f vs reported %.2f' % (n, model_value, reported_plddt[n]))

            if mismatches:
                print('  WARNING: reported pLDDT values do not match: %s' % '; '.join(mismatches))
            else:
                print('  reported pLDDT values match the model exactly')
            summary_rows.append([variant, os.path.basename(path), chain_id, c['n_residues'], c['first'], c['last'], '%.2f' % statistics.mean(plddt.values()), '%.2f' % statistics.mean(local_vals), '%.2f' % plddt.get(position_of_interest, float('nan')), '; '.join(c['missing']) or 'none'])

            if path.endswith('.cif'):
                bio_parser = MMCIFParser(QUIET=True)
            else:
                bio_parser = PDBParser(QUIET=True)
            bio_chain = bio_parser.get_structure('s', path)[0][chain_id]
            region_rows, by_number, residues = region_context(bio_chain)

            for row in region_rows:
                context_rows.append([variant, os.path.basename(path)] + row)
            r168 = by_number[position_of_interest]

            for row in region_rows:
                if row[0].startswith('Linker before'):
                    _converted_14.append(row)
            linker_rows = _converted_14

            if linker_rows:
                linker_mean_plddt = linker_rows[0][4]
                print('  linker 136-181: mean pLDDT %.2f (R168 pLDDT %.2f)' % (linker_mean_plddt, plddt[position_of_interest]))
            contacts = linker_kh_contacts(bio_chain)
            print('  R168 is %s at position %d, %d residues before the annotated disordered region (182)' % (r168.resname, position_of_interest, 182 - position_of_interest))
            print('  linker residues in contact with the KH domain: %d of %d' % (len(contacts['linker_contacting_kh']), contacts['linker_n']))

            if contacts['linker_contacting_kh']:
                for n, rn, kmin, kmax, _c in contacts['linker_contacting_kh']:
                    _converted_15.append('%d (%s, KH %d-%d)' % (n, rn, kmin, kmax))
                spans = _converted_15
                print('    contacting linker residues: %s' % ', '.join(spans))
            ck = contacts['closest_kh_to_168']
            print('  closest KH-domain residue to R168: %d (%s) at %.1f A' % (ck[0], ck[1], ck[2]))
    path = os.path.join(out_folder, 'structural_metadata_summary.tsv')

    with open(path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['variant', 'file', 'chain', 'n_residues', 'first_residue', 'last_residue', 'mean_plddt_all', 'mean_plddt_158_178', 'plddt_168', 'missing_residues'])

        for row in summary_rows:
            writer.writerow(row)
    print('\n  wrote', path)
    path = os.path.join(out_folder, 'region_context.tsv')

    with open(path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['variant', 'file', 'region', 'start', 'end', 'n_residues', 'mean_plddt', 'min_plddt', 'mean_sasa_A2', 'pct_exposed_gt_30A2'])

        for row in context_rows:
            writer.writerow(row)
    print('  wrote', path)

if __name__ == '__main__':
    main()
