import os

from collections import defaultdict

MATCHES_FILE = 'results/motif_target_analysis/02_motif_scanning/all_motif_matches.tsv'

ANNOTATED_FILE = 'results/motif_target_analysis/03_transcript_annotation/annotated_motif_sites.tsv'

ACCESS_FILE = 'results/motif_target_analysis/04_rna_accessibility/accessibility_results.tsv'

OUT_DIR = 'results/motif_target_analysis/05_functional_analysis'

def read_tsv(path):
    rows = []

    with open(path) as f:
        header = f.readline().strip().split('\t')

        for i in f:
            cols = i.rstrip('\n').split('\t')
            row = dict(zip(header, cols))
            rows.append(row)
    return rows

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.exists(ANNOTATED_FILE):
        print('ERROR: annotated file not found:', ANNOTATED_FILE)
        return
    annotated = read_tsv(ANNOTATED_FILE)
    acc_lookup = {}

    if os.path.exists(ACCESS_FILE):
        for row in read_tsv(ACCESS_FILE):
            key = (row['Transcript_ID'], row['Transcript_start'], row['Motif'])
            acc_lookup[key] = row.get('Accessibility_heuristic', '')
    isoform_counts = defaultdict(set)

    for row in annotated:
        key = (row['Gene'], row['Motif_name'], row['Matched_sequence'])
        isoform_counts[key].add(row['Transcript_ID'])
    candidates = []

    for row in annotated:
        gene = row['Gene']
        motif = row['Motif_name']
        tid = row['Transcript_ID']
        matched = row['Matched_sequence']
        feature = row.get('Feature', '')
        junction = row.get('Junction_proximity', 'none')
        score = 0
        reasons = []

        if feature in ('CDS', "3' UTR"):
            score += 1
            reasons.append('in_%s' % feature.replace(' ', '_').replace("'", ''))

        if junction != 'none':
            score += 1
            reasons.append('near_splice_junction')
        key = (gene, motif, matched)
        n_isoforms = len(isoform_counts.get(key, set()))

        if n_isoforms > 1:
            score += 1
            reasons.append('shared_across_%d_isoforms' % n_isoforms)
        acc_key = (tid, row['Transcript_start'], motif)
        acc_str = acc_lookup.get(acc_key, '')

        if acc_str:
            try:
                acc = float(acc_str)

                if acc >= 0.5:
                    score += 1
                    reasons.append('potentially_accessible')
            except ValueError:
                pass

        if score >= 3:
            confidence = 'high_computational'

        elif score >= 2:
            confidence = 'moderate_computational'
        else:
            confidence = 'exploratory'
        candidates.append({'Gene': gene, 'Motif': motif, 'Transcript_ID': tid, 'Transcript_start': row['Transcript_start'], 'Transcript_end': row['Transcript_end'], 'Matched_sequence': matched, 'Feature': feature, 'Junction_proximity': junction, 'N_isoforms': n_isoforms, 'Accessibility': acc_str, 'Priority_score': score, 'Max_score': 4, 'Confidence': confidence, 'Score_reasons': ';'.join(reasons)})
    candidates.sort(key=lambda r: (-int(r['Priority_score']), r['Gene'], r['Motif']))
    out_path = os.path.join(OUT_DIR, 'candidate_sites.tsv')
    columns = list(candidates[0].keys()) if candidates else ['no_data']

    with open(out_path, 'w') as f:
        f.write('\t'.join(columns) + '\n')

        for row in candidates:
            for c in columns:
                _converted_1.append(str(row[c]))
            f.write('\t'.join(_converted_1) + '\n')
    print('wrote', out_path, '(%d rows)' % len(candidates))
    gene_summary = defaultdict(lambda: {'WT5': 0, 'R5': 0, 'total': 0})

    for row in candidates:
        gene_summary[row['Gene']]['total'] += 1

        if row['Motif'] in ('WT5', 'R5'):
            gene_summary[row['Gene']][row['Motif']] += 1
    gene_path = os.path.join(OUT_DIR, 'target_gene_summary.tsv')

    with open(gene_path, 'w') as f:
        f.write('Gene\tWT5_matches\tR5_matches\tTotal\n')

        for gene in sorted(gene_summary):
            s = gene_summary[gene]
            f.write('%s\t%d\t%d\t%d\n' % (gene, s['WT5'], s['R5'], s['total']))
    print('wrote', gene_path)
    isoform_counts_out = defaultdict(int)

    for row in candidates:
        isoform_counts_out[row['Transcript_ID']] += 1
    isoform_path = os.path.join(OUT_DIR, 'isoform_summary.tsv')

    with open(isoform_path, 'w') as f:
        f.write('Transcript_ID\tN_candidate_sites\n')

        for tid in sorted(isoform_counts_out):
            f.write('%s\t%d\n' % (tid, isoform_counts_out[tid]))
    print('wrote', isoform_path)
    priority_counts = defaultdict(int)

    for row in candidates:
        priority_counts[row['Confidence']] += 1
    priority_path = os.path.join(OUT_DIR, 'target_prioritization.tsv')

    with open(priority_path, 'w') as f:
        f.write('Confidence\tN_sites\n')

        for level in ('high_computational', 'moderate_computational', 'exploratory'):
            f.write('%s\t%d\n' % (level, priority_counts.get(level, 0)))
    print('wrote', priority_path)
    interp_path = os.path.join(OUT_DIR, 'functional_interpretation.md')
    f = open(interp_path, 'w')
    f.write('# Functional interpretation of candidate sites\n\n')
    f.write('## Prioritization method\n\n')
    f.write('Each site was scored on a 0-4 scale:\n\n')
    f.write("1. +1 if the site is in a CDS or 3' UTR\n")
    f.write('2. +1 if the site is near a splice junction (within 50 bp)\n')
    f.write('3. +1 if the site is shared across multiple transcript isoforms\n')
    f.write('4. +1 if the heuristic accessibility score >= 0.5\n\n')
    f.write('A score of 3-4 is high computational confidence; 2 is moderate; 0-1 is exploratory.\n\n')
    f.write('## Results\n\n')
    f.write('| Confidence | N sites |\n')
    f.write('|------------|---------|\n')

    for level in ('high_computational', 'moderate_computational', 'exploratory'):
        f.write('| %s | %d |\n' % (level, priority_counts.get(level, 0)))
    f.write('\n## Top candidates\n\n')
    f.write('| Gene | Motif | Transcript | Score | Feature | Reasons |\n')
    f.write('|------|-------|------------|-------|---------|---------|\n')

    for row in candidates[:10]:
        f.write('| %s | %s | %s | %s | %s | %s |\n' % (row['Gene'], row['Motif'], row['Transcript_ID'], row['Priority_score'], row['Feature'], row['Score_reasons']))
    f.write('\n## Caveats\n\n')
    f.write('- The prioritization score is a simple heuristic and does not represent experimentally validated binding.\n')
    f.write('- The accessibility measure is a GC-content proxy, not a real RNA structure prediction.\n')
    f.write('- All results must be interpreted with extreme caution.\n')
    f.close()
    print('wrote', interp_path)

if __name__ == '__main__':
    main()
