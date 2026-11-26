import os

from collections import defaultdict

MATCHES_FILE = 'results/motif_target_analysis/02_motif_scanning/all_motif_matches.tsv'

OUT_DIR = 'results/motif_target_analysis/03_transcript_annotation'

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.exists(MATCHES_FILE):
        print('ERROR: matches file not found:', MATCHES_FILE)
        return
    matches = []

    with open(MATCHES_FILE) as f:
        header = f.readline().strip().split('\t')

        for i in f:
            cols = i.rstrip('\n').split('\t')
            row = dict(zip(header, cols))
            matches.append(row)
    sites_by_gene = defaultdict(lambda: defaultdict(list))

    for i in matches:
        gene = i['Gene']
        motif = i['Motif_name']
        matched_seq = i['Matched_sequence']
        key = (motif, matched_seq)
        sites_by_gene[gene][key].append(i['Transcript_ID'])
    rows = []

    for gene in sorted(sites_by_gene):
        for (motif, matched_seq), transcript_ids in sorted(sites_by_gene[gene].items()):
            unique_tids = sorted(set(transcript_ids))
            n_isoforms = len(unique_tids)

            if n_isoforms > 1:
                isoform_status = 'shared_across_isoforms'
            else:
                isoform_status = 'isoform_specific'
            rows.append({'Gene': gene, 'Motif': motif, 'Matched_sequence': matched_seq, 'N_isoforms_with_site': n_isoforms, 'Isoform_status': isoform_status, 'Transcript_IDs': ';'.join(unique_tids)})
    out_path = os.path.join(OUT_DIR, 'isoform_comparison.tsv')
    columns = list(rows[0].keys()) if rows else ['no_data']

    with open(out_path, 'w') as f:
        f.write('\t'.join(columns) + '\n')

        for row in rows:
            for c in columns:
                _converted_1.append(str(row[c]))
            f.write('\t'.join(_converted_1) + '\n')
    print('wrote', out_path, '(%d rows)' % len(rows))

    for r in rows:
        if r['Isoform_status'] == 'shared_across_isoforms':
            _converted_2.append(1)
    n_shared = sum(_converted_2)

    for r in rows:
        if r['Isoform_status'] == 'isoform_specific':
            _converted_3.append(1)
    n_specific = sum(_converted_3)
    print('  shared across isoforms: %d' % n_shared)
    print('  isoform-specific: %d' % n_specific)

if __name__ == '__main__':
    main()
