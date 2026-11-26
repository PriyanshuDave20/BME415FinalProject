import os

from collections import defaultdict

SELECTED = 'results/target_isoform_docking/selected_isoforms/selected_isoforms.tsv'

MATCHES = 'results/target_isoform_docking/motif_scanning/all_motif_matches.tsv'

SUMMARY_BY_GENE = 'results/target_isoform_docking/motif_scanning/motif_summary_by_gene.tsv'

CANDIDATES = 'results/target_isoform_docking/docking_targets/candidate_regions.tsv'

OUT_DIR = 'results/target_isoform_docking/reports'

def read_tsv(path):
    rows = []

    if not os.path.exists(path):
        return rows

    with open(path) as f:
        header = f.readline().strip().split('\t')

        for i in f:
            cols = i.rstrip('\n').split('\t')
            row = dict(zip(header, cols))
            rows.append(row)
    return rows

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    selected = read_tsv(SELECTED)
    matches = read_tsv(MATCHES)
    gene_summary = read_tsv(SUMMARY_BY_GENE)
    candidates = read_tsv(CANDIDATES)
    gene_counts = defaultdict(int)

    for r in selected:
        gene_counts[r['Gene']] += 1

    for m in matches:
        if m.get('Motif_name') == 'WT5':
            _converted_1.append(m)
    wt5_matches = _converted_1

    for m in matches:
        if m.get('Motif_name') == 'R5':
            _converted_2.append(m)
    r5_matches = _converted_2
    p = os.path.join(OUT_DIR, 'target_isoform_selection_report.md')

    with open(p, 'w') as f:
        f.write('# Target isoform selection report\n\n')
        f.write('## Isoform counts\n\n')
        f.write('| Gene | N selected |\n|------|------------|\n')

        for i in ('NRXN1', 'NRXN2', 'NRXN3', 'AR_V7'):
            n = gene_counts.get(i, 0)
            min_req = 3 if i != 'AR_V7' else 1
            status = 'PASS' if i == 'AR_V7' and n >= 1 or (i != 'AR_V7' and n >= 3) else 'FAIL'
            f.write('| %s | %d (min %d required) %s |\n' % (i, n, min_req, status))
        f.write('\n## AR-V7\n\n')

        for r in selected:
            if r['Transcript_ID'] == 'ENST00000504326.5':
                _converted_3.append(r)
        arv7 = _converted_3
        f.write('- AR-V7 ENST00000504326.5 included: %s\n' % ('yes' if arv7 else 'no'))

        if arv7:
            f.write('- Length: %s nt\n' % arv7[0].get('Length', '?'))
        f.write('\n## Motif counts\n\n')
        f.write('| Gene | WT5 | R5 | Total |\n')
        f.write('|------|-----|-----|-------|\n')

        for row in gene_summary:
            f.write('| %s | %s | %s | %s |\n' % (row.get('Gene', ''), row.get('WT5_matches', ''), row.get('R5_matches', ''), row.get('Total', '')))
    print('wrote', p)
    p = os.path.join(OUT_DIR, 'motif_scanning_report.md')

    with open(p, 'w') as f:
        f.write('# Motif scanning report\n\n')
        f.write('## Method\n\n')
        f.write('IUPAC-aware consensus matching was performed on all accepted transcript sequences. Only the forward strand of each transcript was scanned.\n\n')
        f.write('- WT5: `CCGTCTCAA` (width 9, all bases defined)\n')
        f.write('- R5:  `GGTCTYGMWC` (width 10, Y=C/T, M=A/C, W=A/T)\n\n')
        f.write('## Results\n\n')
        f.write('Total WT5 matches: %d\n' % len(wt5_matches))
        f.write('Total R5 matches: %d\n\n' % len(r5_matches))
        f.write('### Matches by gene\n\n')
        f.write('| Gene | WT5 | R5 |\n|------|-----|-----|\n')

        for row in gene_summary:
            f.write('| %s | %s | %s |\n' % (row.get('Gene', ''), row.get('WT5_matches', ''), row.get('R5_matches', '')))

        if wt5_matches:
            f.write('\n### WT5 matched sequences\n\n')
            f.write('| Transcript | Start | End | Matched |\n')
            f.write('|-----------|-------|-----|---------|\n')

            for m in wt5_matches:
                f.write('| %s | %s | %s | `%s` |\n' % (m['Transcript_ID'], m['Transcript_start'], m['Transcript_end'], m['Matched_sequence']))

        if r5_matches:
            f.write('\n### R5 matched sequences\n\n')
            f.write('| Transcript | Start | End | Matched |\n')
            f.write('|-----------|-------|-----|---------|\n')

            for m in r5_matches:
                f.write('| %s | %s | %s | `%s` |\n' % (m['Transcript_ID'], m['Transcript_start'], m['Transcript_end'], m['Matched_sequence']))
        f.write('\n## Caveats\n\n')
        f.write('- Consensus matching is not PWM scoring.\n')
        f.write('- WT5 and R5 have different lengths and degeneracy; raw hit counts are not directly comparable.\n')
        f.write('- No background model was used.\n')
        f.write('- Only the forward strand was scanned.\n')
    print('wrote', p)
    p = os.path.join(OUT_DIR, 'docking_target_summary.md')

    with open(p, 'w') as f:
        f.write('# Docking target summary\n\n')
        f.write('## Candidate regions extracted\n\n')
        f.write('Total: %d\n\n' % len(candidates))

        if candidates:
            f.write('| Candidate_ID | Gene | Transcript | Motif | Matched | Length | Complete |\n')
            f.write('|-------------|------|-----------|-------|---------|--------|----------|\n')

            for i in candidates:
                f.write('| %s | %s | %s | %s | `%s` | %s | %s |\n' % (i.get('Candidate_ID', ''), i.get('Gene', ''), i.get('Transcript_ID', ''), i.get('Motif', ''), i.get('Matched_sequence', ''), i.get('Extracted_length', ''), i.get('Complete', '')))
        f.write('\n## Output files\n\n')
        f.write('- `docking_targets/candidate_regions.tsv`\n')
        f.write('- `docking_targets/NRXN1/*.fa`\n')
        f.write('- `docking_targets/NRXN2/*.fa`\n')
        f.write('- `docking_targets/NRXN3/*.fa`\n')
        f.write('- `docking_targets/AR_V7/*.fa`\n')
        f.write('- `rna_sequences/full_transcripts/all_full_transcripts_rna.fa`\n')
        f.write('- `rna_sequences/motif_positive_regions/motif_positive_regions_rna.fa`\n')
    print('wrote', p)

if __name__ == '__main__':
    main()
