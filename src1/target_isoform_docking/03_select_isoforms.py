import sys

import os

TRANSCRIPT_DIR = 'results/target_isoform_docking/transcripts'

OUT_DIR = 'results/target_isoform_docking/selected_isoforms'

SCAN_PATH = 'results/target_isoform_docking/motif_scanning/all_motif_matches.tsv'

MIN_NRXN = 3

AR_V7_ID = 'ENST00000504326.5'

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

def parse_header(header):
    fields = {'transcript_id': '', 'gene': '', 'type': '', 'chrom': '', 'strand': '', 'length': 0, 'motif_status': 'motif_negative'}
    tokens = header.split()

    if tokens:
        fields['transcript_id'] = tokens[0]

    for i in tokens[1:]:
        if '=' in i:
            key, val = i.split('=', 1)

            if key == 'gene':
                fields['gene'] = val

            elif key == 'type':
                fields['type'] = val

            elif key == 'chrom':
                if ':' in val:
                    fields['chrom'], fields['strand'] = val.rsplit(':', 1)
                else:
                    fields['chrom'] = val

            elif key == 'length':
                fields['length'] = int(val)
    return fields

def validate_sequence(seq):
    valid = set('ACGTNacgtn')
    bad = set(seq.upper()) - valid
    return len(bad) == 0

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    problems = []
    all_isoforms = {}

    for i in sorted(os.listdir(TRANSCRIPT_DIR)):
        if not i.endswith('.fa'):
            continue
        gene_file = i.replace('_transcripts.fa', '')
        gene = gene_file

        if gene == 'AR':
            gene = 'AR_V7'
        path = os.path.join(TRANSCRIPT_DIR, i)
        records = parse_fasta(path)
        isoforms = []
        seen_seqs = {}

        for header, seq in records:
            info = parse_header(header)
            tid = info['transcript_id']
            upper = seq.upper()
            seq_ok = validate_sequence(upper)

            if not seq_ok:
                print('  WARNING: %s contains invalid characters' % tid)

            if upper in seen_seqs:
                info['duplicate_of'] = seen_seqs[upper]
            else:
                seen_seqs[upper] = tid
                info['duplicate_of'] = ''
            info['gene'] = gene
            info['seq'] = upper
            info['header'] = header
            info['motif_status'] = 'motif_negative'
            isoforms.append(info)
        all_isoforms[gene] = isoforms

    for gene in ('NRXN1', 'NRXN2', 'NRXN3'):
        n = len(all_isoforms.get(gene, []))

        if n < MIN_NRXN:
            problems.append('%s has only %d isoforms (need >= %d)' % (gene, n, MIN_NRXN))
    ar_list = all_isoforms.get('AR_V7', [])
    ar_v7_found = False

    for rec in ar_list:
        if rec['transcript_id'] == AR_V7_ID:
            ar_v7_found = True

    if not ar_v7_found:
        problems.append('AR-V7 %s NOT found in AR transcript set' % AR_V7_ID)

    if problems:
        print('VALIDATION PROBLEMS:')

        for i in problems:
            print('  -', i)
        print('Error: SystemExit(1)')
        sys.exit(1)
    match_scans = {}

    if os.path.exists(SCAN_PATH):
        with open(SCAN_PATH) as f:
            hdr = f.readline().strip().split('\t')

            for line in f:
                cols = line.rstrip('\n').split('\t')
                row = dict(zip(hdr, cols))
                tid = row.get('Transcript_ID', '')
                motif = row.get('Motif_name', '')

                if tid not in match_scans:
                    match_scans[tid] = set()
                match_scans[tid].add(motif)

    for gene in ('NRXN1', 'NRXN2', 'NRXN3', 'AR_V7'):
        for rec in all_isoforms.get(gene, []):
            tid = rec['transcript_id']
            motifs = match_scans.get(tid, set())

            if len(motifs) >= 2:
                rec['motif_status'] = 'motif_candidate'

            elif len(motifs) == 1:
                rec['motif_status'] = 'motif_positive'
            else:
                rec['motif_status'] = 'motif_negative'
    tsv_path = os.path.join(OUT_DIR, 'selected_isoforms.tsv')

    with open(tsv_path, 'w') as f:
        f.write('Gene\tTranscript_ID\tType\tChrom\tStrand\tLength\tDuplicate_of\tMotif_status\n')

        for gene in ('NRXN1', 'NRXN2', 'NRXN3', 'AR_V7'):
            for rec in all_isoforms.get(gene, []):
                f.write('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\n' % (rec['gene'], rec['transcript_id'], rec['type'], rec['chrom'], rec['strand'], len(rec['seq']), rec['duplicate_of'], rec['motif_status']))
    print('wrote', tsv_path)
    ranked_path = os.path.join(OUT_DIR, 'ranked_isoforms.tsv')

    with open(ranked_path, 'w') as f:
        f.write('Rank\tGene\tTranscript_ID\tMotif_Status\tType\tLength\n')

        for gene in ('NRXN1', 'NRXN2', 'NRXN3', 'AR_V7'):
            isoforms = all_isoforms.get(gene, [])
            rank_order = {'motif_candidate': 0, 'motif_positive': 1, 'motif_negative': 2}
            sorted_iso = sorted(isoforms, key=lambda r: rank_order.get(r['motif_status'], 99))

            for i, rec in enumerate(sorted_iso, 1):
                f.write('%d\t%s\t%s\t%s\t%s\t%d\n' % (i, gene, rec['transcript_id'], rec['motif_status'], rec['type'], len(rec['seq'])))
    print('wrote', ranked_path)
    fa_path = os.path.join(OUT_DIR, 'selected_isoforms.fa')

    with open(fa_path, 'w') as f:
        for gene in ('NRXN1', 'NRXN2', 'NRXN3', 'AR_V7'):
            for rec in all_isoforms.get(gene, []):
                f.write('>%s gene=%s type=%s chrom=%s:%s length=%d motif_status=%s\n' % (rec['transcript_id'], rec['gene'], rec['type'], rec['chrom'], rec['strand'], len(rec['seq']), rec['motif_status']))
                s = rec['seq']

                for i in range(0, len(s), 60):
                    f.write(s[i:i + 60] + '\n')
    print('wrote', fa_path)
    neg_path = os.path.join(OUT_DIR, 'motif_negative_isoforms.tsv')
    neg_count = 0

    with open(neg_path, 'w') as f:
        f.write('Gene\tTranscript_ID\tType\tLength\n')

        for gene in ('NRXN1', 'NRXN2', 'NRXN3', 'AR_V7'):
            for rec in all_isoforms.get(gene, []):
                if rec['motif_status'] == 'motif_negative':
                    f.write('%s\t%s\t%s\t%d\n' % (rec['gene'], rec['transcript_id'], rec['type'], len(rec['seq'])))
                    neg_count += 1
    print('wrote', neg_path, '(%d motif-negative isoforms)' % neg_count)
    val_path = os.path.join(OUT_DIR, 'isoform_validation.tsv')

    with open(val_path, 'w') as f:
        f.write('Gene\tTranscript_ID\tValid_sequence\tDuplicate\tIs_AR_V7\tCount_check\tMotif_status\n')

        for gene in ('NRXN1', 'NRXN2', 'NRXN3', 'AR_V7'):
            isoforms = all_isoforms.get(gene, [])

            for r in isoforms:
                if validate_sequence(r['seq']):
                    _converted_1.append(1)
            n_valid = sum(_converted_1)

            for rec in isoforms:
                valid = validate_sequence(rec['seq'])
                is_arv7 = 'yes' if rec['transcript_id'] == AR_V7_ID else 'no'
                count_ok = 'yes' if gene not in ('NRXN1', 'NRXN2', 'NRXN3') or len(isoforms) >= MIN_NRXN else 'no'
                f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (rec['gene'], rec['transcript_id'], 'yes' if valid else 'no', rec['duplicate_of'] or 'no', is_arv7, count_ok, rec['motif_status']))

            if gene in ('NRXN1', 'NRXN2', 'NRXN3'):
                print('  %s: %d isoforms (min %d required) -> %s' % (gene, len(isoforms), MIN_NRXN, 'PASS' if len(isoforms) >= MIN_NRXN else 'FAIL'))
        print('  AR_V7: ENST00000504326.5 found: %s' % ('yes' if ar_v7_found else 'no'))
    print('wrote', val_path)

    for gene in ('NRXN1', 'NRXN2', 'NRXN3', 'AR_V7'):
        isoforms = all_isoforms.get(gene, [])

        for r in isoforms:
            if r['motif_status'] == 'motif_positive':
                _converted_2.append(1)
        pos = sum(_converted_2)

        for r in isoforms:
            if r['motif_status'] == 'motif_negative':
                _converted_3.append(1)
        neg = sum(_converted_3)

        for r in isoforms:
            if r['motif_status'] == 'motif_candidate':
                _converted_4.append(1)
        cand = sum(_converted_4)
        print('  %s: %d positive, %d negative, %d candidate' % (gene, pos, neg, cand))
    summary = os.path.join(OUT_DIR, 'isoform_validation_summary.txt')

    with open(summary, 'w') as f:
        for gene in ('NRXN1', 'NRXN2', 'NRXN3', 'AR_V7'):
            n = len(all_isoforms.get(gene, []))
            f.write('%s: %d isoforms\n' % (gene, n))
        f.write('AR_V7 found: %s\n' % ar_v7_found)
        f.write('All minimum requirements: %s\n' % ('PASS' if not problems else 'FAIL'))
    print('wrote', summary)

if __name__ == '__main__':
    main()
