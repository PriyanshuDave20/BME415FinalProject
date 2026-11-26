import os

SELECTED_FA = 'results/target_isoform_docking/selected_isoforms/selected_isoforms.fa'

CANDIDATES = 'results/target_isoform_docking/docking_targets/candidate_regions.tsv'

OUT_FULL = 'results/target_isoform_docking/rna_sequences/full_transcripts'

OUT_POS = 'results/target_isoform_docking/rna_sequences/motif_positive_regions'

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

def dna_to_rna(seq):
    return seq.upper().replace('T', 'U')

def write_fasta(path, records):
    with open(path, 'w') as f:
        for header, seq in records:
            f.write('>%s\n' % header)

            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + '\n')

def main():
    os.makedirs(OUT_FULL, exist_ok=True)
    os.makedirs(OUT_POS, exist_ok=True)
    full_rna = []

    for header, seq in parse_fasta(SELECTED_FA):
        rna = dna_to_rna(seq)
        invalid = set(rna) - set('AUCGN')

        if invalid:
            print('  WARNING: invalid RNA characters %s in %s' % (invalid, header.split()[0]))
        tid = header.split()[0]
        gene = 'unknown'

        for j in header.split():
            if j.startswith('gene='):
                gene = j.split('=', 1)[1]
        full_rna.append(('%s RNA full_transcript gene=%s length=%d' % (tid, gene, len(rna)), rna))
    out_path = os.path.join(OUT_FULL, 'all_full_transcripts_rna.fa')
    write_fasta(out_path, full_rna)
    print('wrote', out_path, '(%d records)' % len(full_rna))
    pos_rna = []

    if os.path.exists(CANDIDATES):
        with open(CANDIDATES) as f:
            header_line = f.readline().strip().split('\t')

            for line in f:
                cols = line.rstrip('\n').split('\t')
                row = dict(zip(header_line, cols))
                dna_seq = row.get('Region_sequence', '')
                rna = dna_to_rna(dna_seq)
                invalid = set(rna) - set('AUCGN')

                if invalid:
                    print('  WARNING: invalid RNA chars %s in %s' % (invalid, row.get('Candidate_ID', '?')))
                gene_folder = row.get('Gene_folder', row.get('Gene', ''))
                tid = row.get('Transcript_ID', '')
                motif = row.get('Motif', '')
                matched = row.get('Matched_sequence', '')
                cid = row.get('Candidate_ID', '')
                pos_rna.append(('%s|%s|%s|%s|%s length=%d' % (gene_folder, tid, motif, matched, cid, len(rna)), rna))
    out2 = os.path.join(OUT_POS, 'motif_positive_regions_rna.fa')
    write_fasta(out2, pos_rna)
    print('wrote', out2, '(%d records)' % len(pos_rna))

if __name__ == '__main__':
    main()
