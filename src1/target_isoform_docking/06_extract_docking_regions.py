import os

SELECTED_FA = 'results/target_isoform_docking/selected_isoforms/selected_isoforms.fa'

MATCHES = 'results/target_isoform_docking/motif_scanning/all_motif_matches.tsv'

OUT_DIR = 'results/target_isoform_docking/docking_targets'

FLANK = 50

R5_CONSENSUS = 'GGTCTYGMWC'

WT5_CONSENSUS = 'CCGTCTCAA'

IUPAC = {'Y': 'CT', 'M': 'AC', 'W': 'AT'}

MAX_MISMATCH = 2

AR_V7_ID = 'ENST00000504326.5'

GENE_IDS = {'NRXN1': 'ENSG00000179915', 'NRXN2': 'ENSG00000110076', 'NRXN3': 'ENSG00000021645', 'AR': 'ENSG00000169083', 'AR_V7': 'ENSG00000169083'}

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

def count_mismatches(sub, consensus):
    mm = 0

    for a, b in zip(sub, consensus):
        if b in IUPAC:
            if a not in IUPAC[b]:
                mm += 1

        elif a != b:
            mm += 1

        if mm > MAX_MISMATCH:
            break
    return mm

def find_closest_match(seq, consensus, tid, exclude_exact=True):
    best = None

    for i in range(len(seq) - len(consensus) + 1):
        sub = seq[i:i + len(consensus)]

        if exclude_exact and count_mismatches(sub, consensus) == 0:
            continue
        mm = count_mismatches(sub, consensus)

        if mm <= MAX_MISMATCH and (best is None or mm < best[0]):
            best = (mm, sub, i + 1, i + len(consensus))
    return best

def make_candidate_row(candidate_id, gene, gene_folder, tid, motif, matched, start, end, full_seq, gene_id):
    seq_len = len(full_seq)
    up = start - 1 - FLANK
    down = end + FLANK
    upstream_available = FLANK
    downstream_available = FLANK

    if up < 0:
        upstream_available = start - 1
        up = 0

    if down > seq_len:
        downstream_available = seq_len - end
        down = seq_len
    region_seq = full_seq[up:down]
    return ({'Candidate_ID': candidate_id, 'Gene': gene, 'Gene_folder': gene_folder, 'Gene_ID': gene_id, 'Transcript_ID': tid, 'Motif': motif, 'Matched_sequence': matched, 'Transcript_start': start, 'Transcript_end': end, 'Region_start_in_transcript': up + 1, 'Region_end_in_transcript': down, 'Upstream_flank': upstream_available, 'Downstream_flank': downstream_available, 'Extracted_length': len(region_seq), 'Complete': 'yes' if upstream_available == FLANK and downstream_available == FLANK else 'truncated', 'Region_sequence': region_seq}, region_seq)

def write_candidate_fa(gene_dir, candidate_id, gene_folder, tid, motif, matched, region_seq):
    os.makedirs(gene_dir, exist_ok=True)
    fa_path = os.path.join(gene_dir, candidate_id + '.fa')

    with open(fa_path, 'w') as f:
        f.write('>%s|%s|%s|%s\n' % (gene_folder, tid, motif, matched))

        for i in range(0, len(region_seq), 60):
            f.write(region_seq[i:i + 60] + '\n')

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    seqs = {}

    for header, seq in parse_fasta(SELECTED_FA):
        tid = parse_header_tid(header)
        seqs[tid] = seq.upper()
    matches = []

    with open(MATCHES) as f:
        header = f.readline().strip().split('\t')

        for line in f:
            cols = line.rstrip('\n').split('\t')
            row = dict(zip(header, cols))
            matches.append(row)
    print('  %d matches loaded' % len(matches))
    candidate_id_counter = 0
    rows = []
    gene_folders = {'NRXN1': 'NRXN1', 'NRXN2': 'NRXN2', 'NRXN3': 'NRXN3', 'AR': 'AR_V7', 'AR_V7': 'AR_V7'}
    gene_id_map = {'NRXN1': 'ENSG00000179915', 'NRXN2': 'ENSG00000110076', 'NRXN3': 'ENSG00000021645', 'AR': 'ENSG00000169083', 'AR_V7': 'ENSG00000169083'}
    genes_with_exact = set()

    for i in matches:
        tid = i['Transcript_ID']
        start = int(i['Transcript_start'])
        end = int(i['Transcript_end'])
        gene = i['Gene']
        motif = i['Motif_name']
        matched = i['Matched_sequence']

        if tid not in seqs:
            print('  WARNING: transcript %s not in sequences' % tid)
            continue
        full_seq = seqs[tid]
        seq_len = len(full_seq)
        up = start - 1 - FLANK
        down = end + FLANK
        upstream_available = FLANK
        downstream_available = FLANK

        if up < 0:
            upstream_available = start - 1
            up = 0

        if down > seq_len:
            downstream_available = seq_len - end
            down = seq_len
        region_seq = full_seq[up:down]
        candidate_id_counter += 1
        candidate_id = 'candidate_%03d' % candidate_id_counter
        gene_folder = gene_folders.get(gene, gene)
        genes_with_exact.add(gene_folder)
        gene_dir = os.path.join(OUT_DIR, gene_folder)
        gid = gene_id_map.get(gene, '')
        write_candidate_fa(gene_dir, candidate_id, gene_folder, tid, motif, matched, region_seq)
        row, region_seq2 = make_candidate_row(candidate_id, gene, gene_folder, tid, motif, matched, start, end, full_seq, gid)
        rows.append(row)

    if 'NRXN1' not in genes_with_exact:
        print('  NRXN1: no exact motif match, finding closest approximate...')
        best = None
        best_tid = None

        for tid, seq in seqs.items():
            if not tid.startswith('ENST') or 'NRXN1' not in str(tid):
                continue
            pass

        for header, seq in parse_fasta(SELECTED_FA):
            tid = parse_header_tid(header)
            gene_tag = ''

            for t in header.split():
                if t.startswith('gene='):
                    gene_tag = t.split('=', 1)[1]

            if gene_tag != 'NRXN1':
                continue
            upper = seq.upper()
            result = find_closest_match(upper, R5_CONSENSUS, tid)

            if result and (best is None or result[0] < best[0]):
                best = result
                best_tid = tid
                best_full_seq = upper

        if best:
            mm, matched, start, end = best
            candidate_id_counter += 1
            cid = 'candidate_%03d' % candidate_id_counter
            gene_folder = 'NRXN1'
            gene_dir = os.path.join(OUT_DIR, gene_folder)
            row, region_seq = make_candidate_row(cid, 'NRXN1', gene_folder, best_tid, 'R5_approx', matched, start, end, best_full_seq, 'ENSG00000179915')
            row['Matched_sequence'] = matched
            rows.append(row)
            write_candidate_fa(gene_dir, cid, gene_folder, best_tid, 'R5_approx', matched, region_seq)
            print('  NRXN1 closest: %s (mm=%d) at %d-%d in %s' % (matched, mm, start, end, best_tid))
        else:
            print('  NRXN1: no approximate match within %d mismatches' % MAX_MISMATCH)

    if 'AR_V7' not in genes_with_exact:
        print('  AR_V7: no exact motif match, adding comparison region...')

        for header, seq in parse_fasta(SELECTED_FA):
            tid = parse_header_tid(header)
            gene_tag = ''

            for t in header.split():
                if t.startswith('gene='):
                    gene_tag = t.split('=', 1)[1]

            if gene_tag != 'AR_V7' or tid != AR_V7_ID:
                continue
            upper = seq.upper()
            mid = len(upper) // 2
            start = max(1, mid - 4)
            end = min(len(upper), start + 9)
            candidate_id_counter += 1
            cid = 'candidate_%03d' % candidate_id_counter
            gene_folder = 'AR_V7'
            gene_dir = os.path.join(OUT_DIR, gene_folder)
            row, region_seq = make_candidate_row(cid, 'AR_V7', gene_folder, tid, 'none', '', start, end, upper, 'ENSG00000169083')
            rows.append(row)
            write_candidate_fa(gene_dir, cid, gene_folder, tid, 'none', '', region_seq)
            print('  AR_V7 comparison: %s at %d-%d (no motif)' % (tid, start, end))
            break
    tsv_path = os.path.join(OUT_DIR, 'candidate_regions.tsv')

    if rows:
        cols = list(rows[0].keys())

        with open(tsv_path, 'w') as f:
            f.write('\t'.join(cols) + '\n')

            for row in rows:
                for c in cols:
                    _converted_1.append(str(row[c]))
                f.write('\t'.join(_converted_1) + '\n')
        print('wrote', tsv_path, '(%d candidate regions)' % len(rows))

if __name__ == '__main__':
    main()
