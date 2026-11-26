import csv

from pathlib import Path

from collections import defaultdict

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

DOCKING_RESULTS = PROJECT_ROOT / 'docking_contact_results'

TARGETS = ['NRXN1', 'NRXN2', 'NRXN3', 'AR_V7']

MIN_DISTANCE = 2.0

def parse_float(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return None

def count_contacts(input_file, output_file):
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        print(f'  WARNING: No data in {input_file}')
        return []
    groups = defaultdict(list)

    for row in rows:
        key = (row.get('target', ''), row['condition'], row['model'], row['cutoff_A'])
        groups[key].append(row)
    results = []

    for (target, condition, model, cutoff), group in sorted(groups.items(), key=lambda x: (x[0][0], x[0][1], int(x[0][2]), float(x[0][3]))):
        original_contacting_residues = len(group)
        valid = []

        for row in group:
            dist = parse_float(row.get('minimum_distance_A', ''))

            if dist is not None and dist >= MIN_DISTANCE:
                valid.append(row)
        excluded_clashing = original_contacting_residues - len(valid)

        for row in valid:
            _converted_1.append(row['protein_residue_number'])
        unique_protein_residues = len(set(_converted_1))
        rna_nucleotides = set()

        for row in valid:
            nucleotides = row.get('RNA_nucleotides_contacted', '')

            for n in nucleotides.split(';'):
                n = n.strip()

                if n:
                    rna_nucleotides.add(n)
        residue_rna_pairs = 0

        for row in valid:
            nucleotides = row.get('RNA_nucleotides_contacted', '')

            for n in nucleotides.split(';'):
                if n.strip():
                    _converted_2.append(n)
            pairs = _converted_2
            residue_rna_pairs += len(pairs)
        results.append({'target': target, 'condition': condition, 'model': model, 'cutoff_A': cutoff, 'original_contacting_residues': original_contacting_residues, 'excluded_clashing_residues': excluded_clashing, 'viable_protein_residues': unique_protein_residues, 'unique_RNA_nucleotides': len(rna_nucleotides), 'viable_residue_RNA_pairs': residue_rna_pairs})

    if results:
        fieldnames = list(results[0].keys())

        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        print(f'  Saved: {output_file} ({len(results)} rows)')
    return results

def main():
    all_results = []

    for target in TARGETS:
        target_dir = DOCKING_RESULTS / target
        input_file = target_dir / 'protein_residue_contacts.csv'
        output_file = target_dir / 'viable_contact_summary.csv'

        if not input_file.exists():
            print(f'  WARNING: {input_file} not found, skipping {target}')
            continue
        print(f'\nTarget: {target}')
        results = count_contacts(input_file, output_file)
        all_results.extend(results)
    combined_file = DOCKING_RESULTS / 'viable_contact_summary_all_targets.csv'

    if all_results:
        fieldnames = list(all_results[0].keys())

        with open(combined_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_results)
        print(f'\nCombined summary saved: {combined_file} ({len(all_results)} rows)')
    print('\n' + '=' * 80)
    print('MEDIAN VIABLE CONTACTS PER TARGET/CONDITION/CUTOFF')
    print('=' * 80)
    groups = defaultdict(list)

    for r in all_results:
        key = (r['target'], r['condition'], r['cutoff_A'])
        groups[key].append(r)
    print(f"{'Target':<10} {'Condition':<10} {'Cutoff':>6} {'Median Residues':>16} {'Median RNA':>12} {'Median Pairs':>14}")
    print('-' * 80)

    for (target, condition, cutoff), rows in sorted(groups.items()):
        for r in rows:
            _converted_3.append(r['viable_protein_residues'])
        residues = sorted(_converted_3)

        for r in rows:
            _converted_4.append(r['unique_RNA_nucleotides'])
        rnas = sorted(_converted_4)

        for r in rows:
            _converted_5.append(r['viable_residue_RNA_pairs'])
        pairs = sorted(_converted_5)

        def median(lst):
            n = len(lst)

            if n == 0:
                return 0

            if n % 2 == 1:
                return lst[n // 2]
            return (lst[n // 2 - 1] + lst[n // 2]) / 2
        print(f'{target:<10} {condition:<10} {cutoff:>6} {median(residues):>16} {median(rnas):>12} {median(pairs):>14}')

if __name__ == '__main__':
    main()
