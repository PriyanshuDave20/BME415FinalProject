import os

TRANSCRIPT_DIR = 'reference/target_transcripts'

MATCHES_FILE = 'results/motif_target_analysis/02_motif_scanning/all_motif_matches.tsv'

OUT_DIR = 'results/motif_target_analysis/04_rna_accessibility'

WINDOW = 50

def parse_fasta(path):
    records = []
    header = None
    parts = []

    for line in open(path):
        line = line.strip()

        if line.startswith('>'):
            if header is not None:
                records.append((header, ''.join(parts)))
            header = line[1:]
            parts = []

        elif line:
            parts.append(line)

    if header is not None:
        records.append((header, ''.join(parts)))
    return records

def parse_header_tid(header):
    return header.split()[0]

def gc_content(seq):
    if not seq:
        return 0.0

    for base in seq:
        if base in 'GC':
            _converted_1.append(1)
    gc = sum(_converted_1)
    return gc / len(seq)

def accessibility_heuristic(seq):
    gc = gc_content(seq)
    return 1.0 - gc

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.exists(MATCHES_FILE):
        print('ERROR: matches file not found:', MATCHES_FILE)
        return
    sequences = {}

    for i in sorted(os.listdir(TRANSCRIPT_DIR)):
        if not i.endswith('.fa'):
            continue
        path = os.path.join(TRANSCRIPT_DIR, i)

        for header, seq in parse_fasta(path):
            tid = parse_header_tid(header)
            sequences[tid] = seq.upper().replace('U', 'T')
    matches = []

    with open(MATCHES_FILE) as f:
        header = f.readline().strip().split('\t')

        for line in f:
            cols = line.rstrip('\n').split('\t')
            row = dict(zip(header, cols))
            matches.append(row)
    results = []

    for i in matches:
        tid = i['Transcript_ID']
        start = int(i['Transcript_start'])
        end = int(i['Transcript_end'])

        if tid not in sequences:
            window_seq = ''
            note = 'transcript_not_found'
        else:
            full_seq = sequences[tid]
            w_start = max(0, start - 1 - WINDOW)
            w_end = min(len(full_seq), end + WINDOW)
            window_seq = full_seq[w_start:w_end]
            note = ''

        if window_seq:
            acc_score = accessibility_heuristic(window_seq)
            gc = gc_content(window_seq)
        else:
            acc_score = ''
            gc = ''
        results.append({'Gene': i['Gene'], 'Transcript_ID': tid, 'Motif': i['Motif_name'], 'Transcript_start': start, 'Transcript_end': end, 'Window_start': w_start + 1 if window_seq else '', 'Window_end': w_end if window_seq else '', 'Window_GC_content': '%.3f' % gc if gc != '' else '', 'Accessibility_heuristic': '%.3f' % acc_score if acc_score != '' else '', 'Note': note})
    out_path = os.path.join(OUT_DIR, 'accessibility_results.tsv')
    columns = list(results[0].keys()) if results else ['no_data']

    with open(out_path, 'w') as f:
        f.write('\t'.join(columns) + '\n')

        for row in results:
            for c in columns:
                _converted_2.append(str(row[c]))
            f.write('\t'.join(_converted_2) + '\n')
    print('wrote', out_path, '(%d rows)' % len(results))
    comparison = {}

    for row in results:
        motif = row['Motif']

        if row['Accessibility_heuristic'] == '':
            continue
        val = float(row['Accessibility_heuristic'])

        if motif not in comparison:
            comparison[motif] = []
        comparison[motif].append(val)
    comp_path = os.path.join(OUT_DIR, 'motif_accessibility_comparison.tsv')

    with open(comp_path, 'w') as f:
        f.write('Motif\tN_sites\tMean_accessibility\tMin\tMax\n')

        for motif in sorted(comparison):
            vals = comparison[motif]

            if vals:
                mean_v = sum(vals) / len(vals)
                f.write('%s\t%d\t%.3f\t%.3f\t%.3f\n' % (motif, len(vals), mean_v, min(vals), max(vals)))
            else:
                f.write('%s\t0\t\t\t\n' % motif)
    print('wrote', comp_path)
    report_path = os.path.join(OUT_DIR, 'accessibility_report.md')
    f = open(report_path, 'w')
    f.write('# RNA accessibility report\n\n')
    f.write('## Important caveat\n\n')
    f.write('ViennaRNA is not installed on this machine. The values below are a simple GC-content-based heuristic, NOT a real RNA secondary structure prediction. They should be interpreted with extreme caution and must not be treated as binding accessibility estimates.\n\n')
    f.write('## Method\n\n')
    f.write('For each motif match, a %d nt window on each side was extracted and its GC content was computed. A heuristic accessibility score of `1 - GC_fraction` was used.\n\n' % WINDOW)
    f.write('## Results by motif\n\n')
    f.write('| Motif | N sites | Mean acc. | Min | Max |\n')
    f.write('|-------|---------|-----------|-----|-----|\n')

    for motif in sorted(comparison):
        vals = comparison[motif]

        if vals:
            mean_v = sum(vals) / len(vals)
            f.write('| %s | %d | %.3f | %.3f | %.3f |\n' % (motif, len(vals), mean_v, min(vals), max(vals)))
    f.write('\n## Limitations\n\n')
    f.write('- No real RNA secondary structure prediction was performed.\n')
    f.write('- GC content is a poor proxy for accessibility.\n')
    f.write('- Proper analysis requires ViennaRNA or a similar tool.\n')
    f.close()
    print('wrote', report_path)

if __name__ == '__main__':
    main()
