import csv

from pathlib import Path

from collections import defaultdict

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

DOCKING_RESULTS = PROJECT_ROOT / 'docking_contact_results'

TARGETS = ['NRXN1', 'NRXN2', 'NRXN3', 'AR_V7']

TARGET_LABELS = {'NRXN1': 'NRXN1', 'NRXN2': 'NRXN2', 'NRXN3': 'NRXN3', 'AR_V7': 'AR-V7'}

def median(lst):
    if not lst:
        return 0
    lst = sorted(lst)
    n = len(lst)

    if n % 2 == 1:
        return lst[n // 2]
    return (lst[n // 2 - 1] + lst[n // 2]) / 2

def read_csv_rows(path):
    if not path.exists():
        return []

    with open(path, 'r') as f:
        return list(csv.DictReader(f))

def main():
    all_summary = {}
    all_viable = {}
    all_residues = {}
    all_residue168 = {}

    for target in TARGETS:
        target_dir = DOCKING_RESULTS / target
        summary_rows = read_csv_rows(target_dir / 'contact_summary.csv')
        viable_rows = read_csv_rows(target_dir / 'viable_contact_summary.csv')

        for condition in ['WT', 'R168C']:
            key = (target, condition)

            for r in summary_rows:
                if r['condition'] == condition:
                    _converted_1.append(r)
            all_summary[key] = _converted_1

            for r in viable_rows:
                if r['condition'] == condition:
                    _converted_2.append(r)
            all_viable[key] = _converted_2
            residue_set = set()
            contact_file = target_dir / 'protein_residue_contacts.csv'
            contact_rows = read_csv_rows(contact_file)

            for row in contact_rows:
                if row['condition'] == condition and row.get('cutoff_A') == '5.0' and (float(row.get('minimum_distance_A', 999)) >= 2.0):
                    residue_set.add(row['protein_residue'])
            all_residues[key] = residue_set
            residue168_info = {'contacts': False, 'nucleotides': [], 'identity': 'N/A'}

            for row in all_summary[key]:
                if row.get('cutoff_A') == '5.0' and row.get('residue168_contacts_RNA') == 'True':
                    residue168_info['contacts'] = True
                    nucleotides = row.get('residue168_RNA_nucleotides', '')

                    if nucleotides:
                        for n in nucleotides.split(';'):
                            n = n.strip()

                            if n:
                                residue168_info['nucleotides'].append(n)

                    if row.get('residue168_identity'):
                        residue168_info['identity'] = row['residue168_identity']
            all_residue168[key] = residue168_info
    target_stats = {}

    for target in TARGETS:
        stats = {}

        for condition in ['WT', 'R168C']:
            key = (target, condition)

            for r in all_viable[key]:
                if r['cutoff_A'] == '5.0':
                    _converted_3.append(r)
            viable_5A = _converted_3

            for r in viable_5A:
                _converted_4.append(int(r['viable_protein_residues']))

            for r in viable_5A:
                _converted_5.append(int(r['unique_RNA_nucleotides']))

            for r in viable_5A:
                _converted_6.append(int(r['viable_residue_RNA_pairs']))
            stats[condition] = {'median_viable_residues': median(_converted_4), 'median_unique_rna': median(_converted_5), 'median_viable_pairs': median(_converted_6), 'total_unique_residues': len(all_residues[key]), 'residue168_contacts': all_residue168[key]['contacts'], 'residue168_identity': all_residue168[key]['identity'], 'residue168_nucleotides': sorted(set(all_residue168[key]['nucleotides']))}
        wt = stats['WT']
        mut = stats['R168C']
        stats['fold_change'] = {'residues': mut['median_viable_residues'] / wt['median_viable_residues'] if wt['median_viable_residues'] > 0 else None, 'rna': mut['median_unique_rna'] / wt['median_unique_rna'] if wt['median_unique_rna'] > 0 else None, 'pairs': mut['median_viable_pairs'] / wt['median_viable_pairs'] if wt['median_viable_pairs'] > 0 else None, 'total_unique_residues': mut['total_unique_residues'] / wt['total_unique_residues'] if wt['total_unique_residues'] > 0 else None}
        target_stats[target] = stats
    output_file = DOCKING_RESULTS / 'mutation_impact_summary_all_targets.md'
    lines = []
    lines.append('# Mutation Impact Summary: R168C vs WT Across All RNA Targets')
    lines.append('')
    lines.append('## Overview')
    lines.append('')
    lines.append('This report summarizes protein-RNA contact analysis for KHDRBS2 WT and R168C mutant')
    lines.append('across four RNA targets: NRXN1, NRXN2, NRXN3, and AR-V7.')
    lines.append('')
    lines.append('- **Docking method**: Rigid-body docking (HDOCKlite-compatible format)')
    lines.append('- **Models per condition**: 10')
    lines.append('- **Cutoff distances analyzed**: 4.0 A and 5.0 A')
    lines.append('- **Clash filter**: contacts with minimum distance < 2.0 A excluded')
    lines.append('- **Primary comparison cutoff**: 5.0 A')
    lines.append('')
    lines.append('---')
    lines.append('')
    lines.append('## Per-Target Contact Summary (5.0 A cutoff, median across 10 models)')
    lines.append('')
    lines.append('| Target | Condition | Median Viable Residues | Median Unique RNA | Median Viable Pairs | Total Unique Residues |')
    lines.append('|--------|-----------|----------------------|-------------------|--------------------|----------------------|')

    for target in TARGETS:
        for condition in ['WT', 'R168C']:
            s = target_stats[target][condition]
            label = TARGET_LABELS[target]
            lines.append(f"| {label} | {condition} | {s['median_viable_residues']} | {s['median_unique_rna']} | {s['median_viable_pairs']} | {s['total_unique_residues']} |")
    lines.append('')
    lines.append('---')
    lines.append('')
    lines.append('## R168C / WT Fold Change (5.0 A cutoff)')
    lines.append('')
    lines.append('| Target | Residues Fold Change | RNA Fold Change | Pairs Fold Change | Total Residues Fold Change |')
    lines.append('|--------|---------------------|-----------------|-------------------|---------------------------|')

    for target in TARGETS:
        fc = target_stats[target]['fold_change']
        label = TARGET_LABELS[target]

        def fmt(v):
            if v is None:
                return 'N/A'
            return f'{v:.2f}x'
        lines.append(f"| {label} | {fmt(fc['residues'])} | {fmt(fc['rna'])} | {fmt(fc['pairs'])} | {fmt(fc['total_unique_residues'])} |")
    lines.append('')
    lines.append('---')
    lines.append('')
    lines.append('## Residue 168 (Mutation Site) RNA Contact Analysis')
    lines.append('')
    lines.append('| Target | Condition | Residue 168 Identity | Contacts RNA? | RNA Nucleotides Contacted |')
    lines.append('|--------|-----------|---------------------|---------------|--------------------------|')

    for target in TARGETS:
        for condition in ['WT', 'R168C']:
            info = all_residue168[target, condition]
            label = TARGET_LABELS[target]
            nucleotides = ';'.join(sorted(set(info['nucleotides']))) if info['nucleotides'] else '—'
            contact_str = 'Yes' if info['contacts'] else 'No'
            lines.append(f"| {label} | {condition} | {info['identity']} | {contact_str} | {nucleotides} |")
    lines.append('')
    lines.append('---')
    lines.append('')
    lines.append('## Overall Interpretation')
    lines.append('')
    increased = 0
    decreased = 0
    unchanged = 0

    for target in TARGETS:
        fc = target_stats[target]['fold_change']['pairs']

        if fc is None:
            unchanged += 1

        elif fc > 1.2:
            increased += 1

        elif fc < 0.8:
            decreased += 1
        else:
            unchanged += 1
    lines.append(f'Across the {len(TARGETS)} RNA targets analyzed:')
    lines.append('')
    lines.append(f'- **{increased} target(s)** show **increased** viable contacts in R168C (>1.2x WT)')
    lines.append(f'- **{decreased} target(s)** show **decreased** viable contacts in R168C (<0.8x WT)')
    lines.append(f'- **{unchanged} target(s)** show **similar** viable contacts in R168C (0.8–1.2x WT)')
    lines.append('')
    fold_changes = {}

    for target in TARGETS:
        fc = target_stats[target]['fold_change']['pairs']

        if fc is not None:
            fold_changes[target] = fc

    if fold_changes:
        max_target = max(fold_changes, key=fold_changes.get)
        min_target = min(fold_changes, key=fold_changes.get)
        lines.append(f'- **Largest increase**: {TARGET_LABELS[max_target]} ({fold_changes[max_target]:.2f}x)')
        lines.append(f'- **Largest decrease**: {TARGET_LABELS[min_target]} ({fold_changes[min_target]:.2f}x)')
        lines.append('')
    lines.append('### Key observations:')
    lines.append('')

    for target in TARGETS:
        s = target_stats[target]
        fc = s['fold_change']
        label = TARGET_LABELS[target]
        direction = 'increased' if (fc['pairs'] or 1) > 1.2 else 'decreased' if (fc['pairs'] or 1) < 0.8 else 'similar'
        lines.append(f"- **{label}**: R168C shows {direction} contacts ({fc['pairs']:.2f}x WT for pairs)")
    lines.append('')
    lines.append('### Residue 168 (mutation site):')
    lines.append('')

    for t in TARGETS:
        if all_residue168[t, 'R168C']['contacts']:
            _converted_7.append(1)
    mut_168_contact_count = sum(_converted_7)

    for t in TARGETS:
        if all_residue168[t, 'WT']['contacts']:
            _converted_8.append(1)
    wt_168_contact_count = sum(_converted_8)
    lines.append(f'- Residue 168 contacts RNA in **{mut_168_contact_count}/{len(TARGETS)}** R168C models')
    lines.append(f'- Residue 168 contacts RNA in **{wt_168_contact_count}/{len(TARGETS)}** WT models')
    lines.append('')

    if mut_168_contact_count > wt_168_contact_count:
        lines.append('The R168C mutation appears to **increase** the frequency with which residue 168')
        lines.append('directly contacts RNA, which may contribute to altered RNA binding specificity.')

    elif mut_168_contact_count < wt_168_contact_count:
        lines.append('The R168C mutation appears to **decrease** the frequency with which residue 168')
        lines.append('directly contacts RNA, potentially reducing direct RNA interaction at this site.')
    else:
        lines.append('The R168C mutation does not substantially change the frequency of direct')
        lines.append('RNA contact at residue 168; any functional impact likely arises from broader')
        lines.append('structural or electrostatic changes rather than this single residue contact.')
    lines.append('')
    lines.append('---')
    lines.append('')
    lines.append(f'*Report generated from docking_contact_results/ across {len(TARGETS)} targets, 10 models each.*')

    with open(output_file, 'w') as f:
        f.write('\n'.join(lines))
    print(f'Summary saved: {output_file}')
    print(f"\n{'=' * 80}")
    print('\n'.join(lines))

if __name__ == '__main__':
    main()
