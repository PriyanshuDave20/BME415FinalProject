import os

import re

TRANSCRIPT_DIR = 'reference/target_transcripts'

MOTIF_DIR = 'results/motif_target_analysis/01_motif_validation'

OUT_DIR = 'results/motif_target_analysis/02_motif_scanning'

IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}

def parse_fasta(path):
    records = []
    header = None
    parts = []

    for i in open(path):
        i = i.strip()

        if i.startswith('>'):
            if header is not None:
                records.append((header, ''.join(parts)))
            header = i[1:]
            parts = []

        elif i:
            parts.append(i)

    if header is not None:
        records.append((header, ''.join(parts)))
    return records

def parse_header(header):
    fields = {'transcript_id': '', 'gene_name': '', 'transcript_type': '', 'chrom': '', 'strand': '', 'length': ''}
    tokens = header.split()

    if tokens:
        fields['transcript_id'] = tokens[0]

    for i in tokens[1:]:
        if '=' in i:
            key, value = i.split('=', 1)

            if key == 'gene':
                fields['gene_name'] = value

            elif key == 'type':
                fields['transcript_type'] = value

            elif key == 'chrom':
                if ':' in value:
                    chrom, strand = value.rsplit(':', 1)
                    fields['chrom'] = chrom
                    fields['strand'] = strand
                else:
                    fields['chrom'] = value

            elif key == 'length':
                fields['length'] = value
    return fields

def build_iupac_regex(consensus):
    pattern = []

    for base in consensus.upper():
        if base in IUPAC:
            alternatives = IUPAC[base]

            if len(alternatives) == 1:
                pattern.append(alternatives)
            else:
                pattern.append('[' + alternatives + ']')
        else:
            pattern.append('.')
    return re.compile(''.join(pattern))

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'U': 'A', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}

    for base in reversed(seq):
        _converted_1.append(complement.get(base, 'N'))
    return ''.join(_converted_1)

def scan_transcript(sequence, regex, motif_name):
    matches = []

    for m in regex.finditer(sequence):
        start = m.start() + 1
        end = m.end()
        matches.append((start, end, m.group(), '+'))
    return matches

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    wt5_consensus = 'CCGTCTCAA'
    r5_consensus = 'GGTCTYGMWC'
    wt5_regex = build_iupac_regex(wt5_consensus)
    r5_regex = build_iupac_regex(r5_consensus)
    print('WT5 consensus regex:', wt5_regex.pattern)
    print('R5  consensus regex:', r5_regex.pattern)
    all_matches = []
    transcript_info = {}

    for i in sorted(os.listdir(TRANSCRIPT_DIR)):
        if not i.endswith('.fa'):
            continue
        path = os.path.join(TRANSCRIPT_DIR, i)

        for header, seq in parse_fasta(path):
            fields = parse_header(header)
            tid = fields['transcript_id']
            transcript_info[tid] = fields
            upper_seq = seq.upper().replace('U', 'T')

            for motif_name, consensus, regex in (('WT5', wt5_consensus, wt5_regex), ('R5', r5_consensus, r5_regex)):
                for start, end, matched, orientation in scan_transcript(upper_seq, regex, motif_name):
                    all_matches.append({'Gene': fields['gene_name'], 'Transcript_ID': tid, 'Transcript_version': tid, 'Motif_name': motif_name, 'Motif_type': 'consensus', 'Transcript_start': start, 'Transcript_end': end, 'Matched_sequence': matched, 'Orientation': orientation, 'Raw_score': 'n/a', 'Normalized_score': 'n/a', 'Threshold': 'n/a', 'Site_classification': 'consensus_match', 'Motif_consensus': consensus})
    print('Total matches found:', len(all_matches))

    if all_matches:
        columns = list(all_matches[0].keys())
        out_path = os.path.join(OUT_DIR, 'all_motif_matches.tsv')

        with open(out_path, 'w') as f:
            f.write('\t'.join(columns) + '\n')

            for row in all_matches:
                for c in columns:
                    _converted_2.append(str(row[c]))
                f.write('\t'.join(_converted_2) + '\n')
        print('wrote', out_path)
    else:
        print('No matches found')

    for motif_name in ('WT5', 'R5'):
        for m in all_matches:
            if m['Motif_name'] == motif_name:
                _converted_3.append(m)
        subset = _converted_3
        out_name = '%s_matches.tsv' % motif_name
        out_path = os.path.join(OUT_DIR, out_name)

        if subset:
            columns = list(subset[0].keys())

            with open(out_path, 'w') as f:
                f.write('\t'.join(columns) + '\n')

                for row in subset:
                    for c in columns:
                        _converted_4.append(str(row[c]))
                    f.write('\t'.join(_converted_4) + '\n')
        else:
            with open(out_path, 'w') as f:
                f.write('no_matches\n')
        print('wrote', out_path, '(%d rows)' % len(subset))
    by_position = {}

    for m in all_matches:
        key = (m['Transcript_ID'], m['Transcript_start'], m['Transcript_end'])
        by_position.setdefault(key, []).append(m['Motif_name'])
    diff_rows = []

    for key, motifs in sorted(by_position.items()):
        tid, start, end = key

        if 'WT5' in motifs and 'R5' in motifs:
            classification = 'shared'

        elif 'WT5' in motifs:
            classification = 'WT5_only'

        elif 'R5' in motifs:
            classification = 'R5_only'
        else:
            classification = 'uncertain'
        diff_rows.append({'Transcript_ID': tid, 'Transcript_start': start, 'Transcript_end': end, 'Motifs_found': ','.join(sorted(set(motifs))), 'Classification': classification, 'Gene': transcript_info[tid]['gene_name']})

    if diff_rows:
        out_path = os.path.join(OUT_DIR, 'differential_sites.tsv')
        columns = list(diff_rows[0].keys())

        with open(out_path, 'w') as f:
            f.write('\t'.join(columns) + '\n')

            for row in diff_rows:
                for c in columns:
                    _converted_5.append(str(row[c]))
                f.write('\t'.join(_converted_5) + '\n')
        print('wrote', out_path)
    summary_path = os.path.join(OUT_DIR, 'motif_scanning_summary.md')
    f = open(summary_path, 'w')
    f.write('# WT5 and R5 target-transcript scanning summary\n\n')
    f.write('## Method\n\n')
    f.write('IUPAC-aware consensus matching was used. Because full PWM-based scanning requires calibrated thresholds that are not available without FIMO, these results are **consensus matches** rather than statistically calibrated motif scores.\n\n')
    f.write('- WT5 consensus: `%s` (width 9)\n' % wt5_consensus)
    f.write('- R5 consensus:  `%s` (width 10)\n\n' % r5_consensus)

    for m in all_matches:
        if m['Motif_name'] == 'WT5':
            _converted_6.append(m)
    wt5_matches = _converted_6

    for m in all_matches:
        if m['Motif_name'] == 'R5':
            _converted_7.append(m)
    r5_matches = _converted_7
    f.write('## Results\n\n')
    f.write('| Gene | WT5 matches | R5 matches |\n')
    f.write('|------|------------|-----------|\n')

    for i in ('AR', 'NRXN1', 'NRXN2', 'NRXN3'):
        for m in wt5_matches:
            if m['Gene'] == i:
                _converted_8.append(1)
        nw = sum(_converted_8)

        for m in r5_matches:
            if m['Gene'] == i:
                _converted_9.append(1)
        nr = sum(_converted_9)
        f.write('| %s | %d | %d |\n' % (i, nw, nr))
    f.write('\nTotal WT5: %d, Total R5: %d\n\n' % (len(wt5_matches), len(r5_matches)))

    for r in diff_rows:
        if r['Classification'] == 'shared':
            _converted_10.append(1)
    n_shared = sum(_converted_10)

    for r in diff_rows:
        if r['Classification'] == 'WT5_only':
            _converted_11.append(1)
    n_wt_only = sum(_converted_11)

    for r in diff_rows:
        if r['Classification'] == 'R5_only':
            _converted_12.append(1)
    n_r_only = sum(_converted_12)
    f.write('## Differential classification\n\n')
    f.write('- Shared positions: %d\n' % n_shared)
    f.write('- WT5-only positions: %d\n' % n_wt_only)
    f.write('- R5-only positions: %d\n' % n_r_only)
    f.write('\n## Important caveat\n\n')
    f.write('These are consensus matches on transcript sequences. They do not represent experimentally validated binding sites and do not indicate binding affinity.\n')
    f.close()
    print('wrote', summary_path)

if __name__ == '__main__':
    main()
