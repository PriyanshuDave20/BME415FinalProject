import sys

import os

import re

SELECTED_FA = 'results/target_isoform_docking/selected_isoforms/selected_isoforms.fa'

SELECTED_TSV = 'results/target_isoform_docking/selected_isoforms/selected_isoforms.tsv'

MATCHES = 'results/target_isoform_docking/motif_scanning/all_motif_matches.tsv'

CANDIDATES = 'results/target_isoform_docking/docking_targets/candidate_regions.tsv'

RNA_FA = 'results/target_isoform_docking/rna_sequences/motif_positive_regions/motif_positive_regions_rna.fa'

LOG_DIR = 'results/target_isoform_docking/logs'

WT5 = 'CCGTCTCAA'

R5 = 'GGTCTYGMWC'

IUPAC = {'Y': 'CT', 'M': 'AC', 'W': 'AT'}

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

def matches_consensus(sub, consensus):
    for a, b in zip(sub, consensus):
        if b in IUPAC:
            if a not in IUPAC[b]:
                return False

        elif a != b:
            return False
    return True

def main():
    os.makedirs(LOG_DIR, exist_ok=True)
    log_lines = []
    results = []
    all_ok = True

    def log(msg):
        log_lines.append(msg)
        print(msg)

    def check(name, ok, detail=''):
        nonlocal all_ok
        status = 'PASS' if ok else 'FAIL'
        results.append((name, status, detail))
        log('[%s] %s %s' % (status, name, detail))

        if not ok:
            all_ok = False
    isoforms = parse_fasta(SELECTED_FA)
    isoform_dict = {}

    for header, seq in isoforms:
        tid = header.split()[0]
        isoform_dict[tid] = (header, seq.upper())
    gene_counts = {}

    for header, seq in isoforms:
        gene = 'unknown'

        for j in header.split():
            if j.startswith('gene='):
                gene = j.split('=', 1)[1]
        gene_counts[gene] = gene_counts.get(gene, 0) + 1

    for gene in ('NRXN1', 'NRXN2', 'NRXN3'):
        n = gene_counts.get(gene, 0)
        check('%s isoform count >= 3' % gene, n >= 3, '(%d accepted, need >= 3)' % n)
    check('AR-V7 ENST00000504326.5 included', 'ENST00000504326.5' in isoform_dict)
    bad_seqs = []

    for tid, (header, seq) in isoform_dict.items():
        invalid = set(seq) - set('ACGTN')

        if invalid:
            bad_seqs.append((tid, invalid))
    check('All accepted transcripts have valid sequences', len(bad_seqs) == 0, '%d bad' % len(bad_seqs) if bad_seqs else 'all OK')
    matches = []

    with open(MATCHES) as f:
        header = f.readline().strip().split('\t')

        for line in f:
            cols = line.rstrip('\n').split('\t')
            row = dict(zip(header, cols))
            matches.append(row)
    bad_matches = 0

    for m in matches:
        tid = m['Transcript_ID']
        start = int(m['Transcript_start'])
        end = int(m['Transcript_end'])
        expected = m['Expected_consensus']
        matched = m['Matched_sequence']

        if tid not in isoform_dict:
            bad_matches += 1
            continue
        sub = isoform_dict[tid][1][start - 1:end]

        if not matches_consensus(sub, expected):
            bad_matches += 1
            continue

        if sub != matched:
            bad_matches += 1
    check('Every motif match is a genuine sequence match', bad_matches == 0, '(%d bad)' % bad_matches if bad_matches else 'all OK')
    candidates = []

    with open(CANDIDATES) as f:
        header = f.readline().strip().split('\t')

        for line in f:
            cols = line.rstrip('\n').split('\t')
            row = dict(zip(header, cols))
            candidates.append(row)
    fa_headers = set()
    fa_records = parse_fasta(RNA_FA)

    for header, seq in fa_records:
        fa_headers.add(header.split()[0])
    cand_ids_in_fa = set()

    for i in fa_headers:
        parts = i.split('|')

        if parts:
            cand_ids_in_fa.add(parts[-1])
    cand_ids_in_tsv = set()
    bad_candidates = 0

    for i in candidates:
        cid = i['Candidate_ID']
        cand_ids_in_tsv.add(cid)
        tid = i['Transcript_ID']

        if tid not in isoform_dict:
            bad_candidates += 1

        if cid not in cand_ids_in_fa:
            bad_candidates += 1
    check('Every docking region maps back to its original transcript', bad_candidates == 0, '(%d bad)' % bad_candidates if bad_candidates else 'all OK')
    check('Every candidate-region metadata record has a FASTA record', cand_ids_in_tsv == cand_ids_in_fa, 'TSV=%d FASTA=%d' % (len(cand_ids_in_tsv), len(cand_ids_in_fa)))
    bad_rna = 0

    for header, seq in fa_records:
        invalid = set(seq.upper()) - set('AUCGN')

        if invalid:
            bad_rna += 1
    check('Every RNA sequence contains only A/U/G/C/N', bad_rna == 0, '(%d bad)' % bad_rna if bad_rna else 'all OK')
    check('Every FASTA record has a unique candidate identifier', len(cand_ids_in_fa) == len(fa_records), '%d unique ids / %d records' % (len(cand_ids_in_fa), len(fa_records)))
    val_path = os.path.join(LOG_DIR, 'validation_results.tsv')

    with open(val_path, 'w') as f:
        f.write('Check\tStatus\tDetail\n')

        for name, status, detail in results:
            f.write('%s\t%s\t%s\n' % (name, status, detail))
    log('wrote ' + val_path)
    log_path = os.path.join(LOG_DIR, 'pipeline.log')

    with open(log_path, 'w') as f:
        f.write('Pipeline run summary\n')
        f.write('=' * 40 + '\n\n')

        for gene in ('NRXN1', 'NRXN2', 'NRXN3', 'AR_V7'):
            f.write('%s: %d isoforms accepted\n' % (gene, gene_counts.get(gene, 0)))

        for m in matches:
            if m.get('Motif_name') == 'WT5':
                _converted_1.append(m)
        f.write('WT5 matches: %d\n' % len(_converted_1))

        for m in matches:
            if m.get('Motif_name') == 'R5':
                _converted_2.append(m)
        f.write('R5 matches: %d\n' % len(_converted_2))
        f.write('Candidate docking regions: %d\n' % len(candidates))
        f.write('Overall: %s\n' % ('PASS' if all_ok else 'FAIL'))
        f.write('\nValidation details:\n')

        for name, status, detail in results:
            f.write('  [%s] %s %s\n' % (status, name, detail))
    log('wrote ' + log_path)

    if not all_ok:
        print('Error: SystemExit(1)')
        sys.exit(1)
    log('\nALL VALIDATION CHECKS PASSED')

if __name__ == '__main__':
    main()
