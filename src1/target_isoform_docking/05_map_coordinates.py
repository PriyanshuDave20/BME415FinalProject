import os

GTF = 'reference/gencode.v44.annotation.gtf'

MATCHES = 'results/target_isoform_docking/motif_scanning/all_motif_matches.tsv'

OUT = 'results/target_isoform_docking/motif_scanning/annotated_matches.tsv'

TARGET_GENES = {'ENSG00000179915': 'NRXN1', 'ENSG00000110076': 'NRXN2', 'ENSG00000021645': 'NRXN3', 'ENSG00000169083': 'AR'}

def parse_attributes(text):
    attrs = {}

    for i in text.strip().split(';'):
        i = i.strip()

        if not i:
            continue
        pieces = i.split(' ', 1)

        if len(pieces) == 2:
            attrs[pieces[0]] = pieces[1].strip('"')
    return attrs

def read_gtf():
    transcripts = {}

    for line in open(GTF):
        if line.startswith('#'):
            continue
        cols = line.rstrip('\n').split('\t')

        if len(cols) < 9:
            continue
        feature = cols[2]

        if feature not in ('exon', 'CDS', 'five_prime_utr', 'three_prime_utr'):
            continue
        attrs = parse_attributes(cols[8])
        gid = attrs.get('gene_id', '').split('.')[0]

        if gid not in TARGET_GENES:
            continue
        tid = attrs.get('transcript_id', '')

        if not tid:
            continue

        if tid not in transcripts:
            transcripts[tid] = {'gene_id': gid, 'strand': cols[6], 'exons': [], 'cds': [], 'utr5': [], 'utr3': []}
        s, e = (int(cols[3]), int(cols[4]))

        if feature == 'exon':
            transcripts[tid]['exons'].append((s, e))

        elif feature == 'CDS':
            transcripts[tid]['cds'].append((s, e))

        elif feature == 'five_prime_utr':
            transcripts[tid]['utr5'].append((s, e))

        elif feature == 'three_prime_utr':
            transcripts[tid]['utr3'].append((s, e))

    for tid in transcripts:
        for j in ('exons', 'cds', 'utr5', 'utr3'):
            transcripts[tid][j].sort()
    return transcripts

def transcript_to_genomic(tid, t_start, t_end, transcripts):
    if tid not in transcripts:
        return (None, None, None, False)
    t = transcripts[tid]
    strand = t['strand']
    exons = sorted(t['exons'])

    if strand == '+':
        cum = 0

        for exon_num, (ex_s, ex_e) in enumerate(exons, 1):
            ex_len = ex_e - ex_s + 1

            if t_start <= cum + ex_len:
                off_s = t_start - cum - 1
                off_e = t_end - cum - 1
                return (ex_s + off_s, ex_s + off_e, exon_num, True)
            cum += ex_len
    else:
        exons_rev = list(reversed(exons))
        cum = 0

        for exon_num, (ex_s, ex_e) in enumerate(exons_rev, 1):
            ex_len = ex_e - ex_s + 1

            if t_start <= cum + ex_len:
                off_s = t_start - cum - 1
                off_e = t_end - cum - 1
                g_s = ex_e - off_s
                g_e = ex_e - off_e

                if g_s > g_e:
                    g_s, g_e = (g_e, g_s)
                return (g_s, g_e, exon_num, True)
            cum += ex_len
    return (None, None, None, False)

def classify(tid, g_s, g_e, transcripts):
    if tid not in transcripts:
        return 'unknown'
    t = transcripts[tid]

    def ov(lst, s, e):
        for a, b in lst:
            if s <= b and e >= a:
                return True
        return False

    if ov(t['utr5'], g_s, g_e):
        return "5' UTR"

    if ov(t['utr3'], g_s, g_e):
        return "3' UTR"

    if ov(t['cds'], g_s, g_e):
        return 'CDS'

    if ov(t['exons'], g_s, g_e):
        return 'exon'
    return 'intronic_or_other'

def nearest_junction(tid, g_s, g_e, transcripts):
    if tid not in transcripts:
        return ('unknown', -1)
    t = transcripts[tid]
    exons = sorted(t['exons'])
    junctions = []

    for i in range(len(exons) - 1):
        junctions.append(exons[i][1])
        junctions.append(exons[i + 1][0])

    for j in junctions:
        _converted_1.append(abs(g_s - j))
    dist = min(_converted_1) if junctions else -1

    for j in junctions:
        _converted_2.append(abs(g_e - j))
    dist2 = min(_converted_2) if junctions else -1

    for x in (dist, dist2):
        if x >= 0:
            _converted_3.append(x)
    d = min(_converted_3) if dist >= 0 or dist2 >= 0 else -1
    near = 'yes' if 0 <= d <= 50 else 'no'
    return (near, d)

def main():
    print('Reading GENCODE GTF (this may take a while)...')
    transcripts = read_gtf()
    print('  %d transcripts loaded' % len(transcripts))
    matches = []

    with open(MATCHES) as f:
        header = f.readline().strip().split('\t')

        for line in f:
            cols = line.rstrip('\n').split('\t')
            row = dict(zip(header, cols))
            matches.append(row)
    print('  %d matches loaded' % len(matches))
    rows = []

    for i in matches:
        tid = i['Transcript_ID']
        t_start = int(i['Transcript_start'])
        t_end = int(i['Transcript_end'])
        g_s, g_e, exon_num, ok = transcript_to_genomic(tid, t_start, t_end, transcripts)

        if ok:
            feature = classify(tid, g_s, g_e, transcripts)
            near, dist = nearest_junction(tid, g_s, g_e, transcripts)
            gene_id = transcripts[tid]['gene_id']
            strand = transcripts[tid]['strand']
        else:
            feature = 'transcript_not_in_gtf'
            near, dist = ('unknown', -1)
            gene_id, strand, exon_num = ('', '', '')
        row = dict(i)
        row['Genomic_start'] = g_s if g_s else ''
        row['Genomic_end'] = g_e if g_e else ''
        row['Exon_number'] = exon_num
        row['Feature'] = feature
        row['Near_splice_junction'] = near
        row['Distance_to_junction'] = dist
        row['Gene_ID'] = gene_id
        row['Strand'] = strand
        rows.append(row)

    if rows:
        cols = list(rows[0].keys())

        with open(OUT, 'w') as f:
            f.write('\t'.join(cols) + '\n')

            for row in rows:
                for c in cols:
                    _converted_4.append(str(row[c]))
                f.write('\t'.join(_converted_4) + '\n')
        print('wrote', OUT, '(%d rows)' % len(rows))

if __name__ == '__main__':
    main()
