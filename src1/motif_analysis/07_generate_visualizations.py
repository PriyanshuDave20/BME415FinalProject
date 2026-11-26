import os

GENE_SUMMARY = 'results/motif_target_analysis/05_functional_analysis/target_gene_summary.tsv'

PRIORITIZATION = 'results/motif_target_analysis/05_functional_analysis/target_prioritization.tsv'

OUT_DIR = 'results/motif_target_analysis/06_visualizations'

def read_tsv(path):
    rows = []

    with open(path) as f:
        header = f.readline().strip().split('\t')

        for i in f:
            cols = i.rstrip('\n').split('\t')
            row = dict(zip(header, cols))
            rows.append(row)
    return rows

def text_bar(label, value, max_value, width=40):
    if max_value <= 0:
        filled = 0
    else:
        filled = int(width * value / max_value)
    bar = '#' * filled + '.' * (width - filled)
    return '%-10s|%s| %d' % (label, bar, value)

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    lines = []
    lines.append('=' * 60)
    lines.append('WT5 vs R5 motif analysis - text summary')
    lines.append('=' * 60)
    lines.append('')

    if os.path.exists(GENE_SUMMARY):
        gene_rows = read_tsv(GENE_SUMMARY)

        for r in gene_rows:
            _converted_1.append(max(int(r['WT5_matches']), int(r['R5_matches'])))
        max_val = max(_converted_1, default=1)
        lines.append('Matches by gene')
        lines.append('-' * 40)

        for row in gene_rows:
            lines.append(text_bar(row['Gene'] + ' WT5', int(row['WT5_matches']), max_val))
            lines.append(text_bar(row['Gene'] + ' R5 ', int(row['R5_matches']), max_val))
        lines.append('')

    if os.path.exists(PRIORITIZATION):
        pri_rows = read_tsv(PRIORITIZATION)

        for r in pri_rows:
            _converted_2.append(int(r['N_sites']))
        max_val = max(_converted_2, default=1)
        lines.append('Candidate sites by confidence')
        lines.append('-' * 40)

        for row in pri_rows:
            lines.append(text_bar(row['Confidence'][:15], int(row['N_sites']), max_val))
        lines.append('')
    lines.append('Note: matplotlib is not installed. No PNG plots were generated. Install matplotlib to generate proper figures.')
    out_path = os.path.join(OUT_DIR, 'text_summary.txt')

    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('wrote', out_path)

if __name__ == '__main__':
    main()
