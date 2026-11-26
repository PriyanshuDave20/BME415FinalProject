import pandas as pd

import sys

input_file = sys.argv[1]

MIN_DISTANCE = 2.0

df = pd.read_csv(input_file)

df['minimum_distance_A'] = pd.to_numeric(df['minimum_distance_A'], errors='coerce')

viable = df[df['minimum_distance_A'] >= MIN_DISTANCE].copy()

results = []

for (condition, model, cutoff), group in df.groupby(['condition', 'model', 'cutoff_A']):
    valid = group[group['minimum_distance_A'] >= MIN_DISTANCE]
    unique_protein_residues = valid['protein_residue_number'].nunique()
    rna_nucleotides = set()

    for nucleotides in valid['RNA_nucleotides_contacted'].dropna():
        for j in str(nucleotides).split(';'):
            if j.strip():
                rna_nucleotides.add(j.strip())
    residue_rna_pairs = 0

    for nucleotides in valid['RNA_nucleotides_contacted'].dropna():
        for n in str(nucleotides).split(';'):
            if n.strip():
                _converted_1.append(n)
        residue_rna_pairs += len(_converted_1)
    results.append({'condition': condition, 'model': model, 'cutoff_A': cutoff, 'original_contacting_residues': len(group), 'excluded_clashing_residues': len(group) - len(valid), 'viable_protein_residues': unique_protein_residues, 'unique_RNA_nucleotides': len(rna_nucleotides), 'viable_residue_RNA_pairs': residue_rna_pairs})

output = pd.DataFrame(results)

output.to_csv('viable_contact_summary.csv', index=False)

print('\nViable contact summary')

print(output.to_string(index=False))

print('\nSaved: viable_contact_summary.csv')

print('\nMedian viable contacts across docking models:')

summary = output.groupby(['condition', 'cutoff_A']).agg({'viable_protein_residues': 'median', 'unique_RNA_nucleotides': 'median', 'viable_residue_RNA_pairs': 'median'}).reset_index()

print(summary.to_string(index=False))
