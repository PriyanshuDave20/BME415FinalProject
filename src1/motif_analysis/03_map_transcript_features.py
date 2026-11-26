import os

GTF_FILE = 'reference/gencode.v44.annotation.gtf'

MATCHES_FILE = 'results/motif_target_analysis/02_motif_scanning/all_motif_matches.tsv'

OUT_DIR = 'results/motif_target_analysis/03_transcript_annotation'

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

def read_gtf_features():
    transcripts = {}

    for line in open(GTF_FILE):
        if line.startswith('#'):
            continue
        cols = line.rstrip('\n').split('\t')

        if len(cols) < 9:
            continue
        feature = cols[2]

        if feature not in ('exon', 'CDS', 'five_prime_utr', 'three_prime_utr', 'transcript'):
            continue
        attrs = parse_attributes(cols[8])
        gene_id = attrs.get('gene_id', '').split('.')[0]

        if gene_id not in TARGET_GENES:
            continue
        tid = attrs.get('transcript_id', '')

        if not tid:
            continue

        if tid not in transcripts:
            transcripts[tid] = {'gene_id': gene_id, 'gene_name': attrs.get('gene_name', ''), 'transcript_type': attrs.get('transcript_type', ''), 'chrom': cols[0], 'strand': cols[6], 'exons': [], 'cds': [], 'utr5': [], 'utr3': []}
        start = int(cols[3])
        end = int(cols[4])

        if feature == 'exon':
            transcripts[tid]['exons'].append((start, end))

        elif feature == 'CDS':
            transcripts[tid]['cds'].append((start, end))

        elif feature == 'five_prime_utr':
            transcripts[tid]['utr5'].append((start, end))

        elif feature == 'three_prime_utr':
            transcripts[tid]['utr3'].append((start, end))

    for tid in transcripts:
        for j in ('exons', 'cds', 'utr5', 'utr3'):
            transcripts[tid][j].sort()
    return transcripts

def transcript_to_genomic(tid, transcript_start, transcript_end, transcripts):
    if tid not in transcripts:
        return (None, None, None, False)
    t = transcripts[tid]
    strand = t['strand']
    exons = sorted(t['exons'])

    if strand == '+':
        cum = 0

        for exon_num, (exon_start, exon_end) in enumerate(exons, 1):
            exon_len = exon_end - exon_start + 1

            if transcript_start <= cum + exon_len:
                offset_start = transcript_start - cum - 1
                genomic_start = exon_start + offset_start
                offset_end = transcript_end - cum - 1
                genomic_end = exon_start + offset_end
                return (genomic_start, genomic_end, exon_num, True)
            cum += exon_len
    else:
        exons_reversed = list(reversed(exons))
        cum = 0

        for exon_num, (exon_start, exon_end) in enumerate(exons_reversed, 1):
            exon_len = exon_end - exon_start + 1

            if transcript_start <= cum + exon_len:
                offset_start = transcript_start - cum - 1
                genomic_start = exon_end - offset_start
                offset_end = transcript_end - cum - 1
                genomic_end = exon_end - offset_end

                if genomic_start > genomic_end:
                    genomic_start, genomic_end = (genomic_end, genomic_start)
                return (genomic_start, genomic_end, exon_num, True)
            cum += exon_len
    return (None, None, None, False)

def classify_feature(tid, g_start, g_end, transcripts):
    if tid not in transcripts:
        return 'unknown'
    t = transcripts[tid]

    def overlaps(interval_list, s, e):
        for istart, iend in interval_list:
            if s <= iend and e >= istart:
                return True
        return False

    if overlaps(t['utr5'], g_start, g_end):
        return "5' UTR"

    if overlaps(t['utr3'], g_start, g_end):
        return "3' UTR"

    if overlaps(t['cds'], g_start, g_end):
        return 'CDS'

    if overlaps(t['exons'], g_start, g_end):
        return 'exon'
    return 'intronic_or_other'

def find_nearby_junction(tid, g_start, g_end, transcripts, max_distance=50):
    if tid not in transcripts:
        return 'none'
    t = transcripts[tid]
    exons = sorted(t['exons'])
    junctions = []

    for i in range(len(exons) - 1):
        donor = exons[i][1]
        acceptor = exons[i + 1][0]
        junctions.append((donor, acceptor))
    near = []

    for donor, acceptor in junctions:
        if abs(g_start - donor) <= max_distance or abs(g_start - acceptor) <= max_distance:
            near.append('near_junction_%d_%d' % (donor, acceptor))

        elif abs(g_end - donor) <= max_distance or abs(g_end - acceptor) <= max_distance:
            near.append('near_junction_%d_%d' % (donor, acceptor))

    if near:
        return ';'.join(near)
    return 'none'

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.exists(GTF_FILE):
        print('ERROR: GTF file not found:', GTF_FILE)
        return

    if not os.path.exists(MATCHES_FILE):
        print('ERROR: Matches file not found:', MATCHES_FILE)
        return
    print('Reading the GTF (this may take a while)...')
    transcripts = read_gtf_features()
    print('  %d transcripts for the target genes' % len(transcripts))
    matches = []

    with open(MATCHES_FILE) as f:
        header = f.readline().strip().split('\t')

        for line in f:
            cols = line.rstrip('\n').split('\t')
            row = dict(zip(header, cols))
            matches.append(row)
    print('  %d motif matches loaded' % len(matches))
    annotated = []
    exon_rows = []
    junction_rows = []

    for i in matches:
        tid = i['Transcript_ID']
        t_start = int(i['Transcript_start'])
        t_end = int(i['Transcript_end'])
        g_start, g_end, exon_num, ok = transcript_to_genomic(tid, t_start, t_end, transcripts)

        if not ok:
            feature = 'transcript_not_in_gtf'
            g_start, g_end, exon_num = ('', '', '')
        else:
            feature = classify_feature(tid, g_start, g_end, transcripts)
        junction_info = find_nearby_junction(tid, g_start, g_end, transcripts) if ok else 'unknown'
        row = dict(i)
        row['Genomic_start'] = g_start
        row['Genomic_end'] = g_end
        row['Exon_number'] = exon_num
        row['Feature'] = feature
        row['Junction_proximity'] = junction_info

        if ok and tid in transcripts:
            row['Gene_ID'] = transcripts[tid]['gene_id']
            row['Transcript_biotype'] = transcripts[tid]['transcript_type']
            row['Strand'] = transcripts[tid]['strand']
        else:
            row['Gene_ID'] = ''
            row['Transcript_biotype'] = ''
            row['Strand'] = ''
        annotated.append(row)

        if ok:
            exon_rows.append({'Transcript_ID': tid, 'Gene': i['Gene'], 'Transcript_start': t_start, 'Transcript_end': t_end, 'Genomic_start': g_start, 'Genomic_end': g_end, 'Exon_number': exon_num, 'Feature': feature})
            junction_rows.append({'Transcript_ID': tid, 'Gene': i['Gene'], 'Motif': i['Motif_name'], 'Junction_proximity': junction_info})

    if annotated:
        columns = list(annotated[0].keys())
        out_path = os.path.join(OUT_DIR, 'annotated_motif_sites.tsv')

        with open(out_path, 'w') as f:
            f.write('\t'.join(columns) + '\n')

            for row in annotated:
                for c in columns:
                    _converted_1.append(str(row[c]))
                f.write('\t'.join(_converted_1) + '\n')
        print('wrote', out_path, '(%d rows)' % len(annotated))

    if exon_rows:
        columns = list(exon_rows[0].keys())
        out_path = os.path.join(OUT_DIR, 'exon_mapping.tsv')

        with open(out_path, 'w') as f:
            f.write('\t'.join(columns) + '\n')

            for row in exon_rows:
                for c in columns:
                    _converted_2.append(str(row[c]))
                f.write('\t'.join(_converted_2) + '\n')
        print('wrote', out_path, '(%d rows)' % len(exon_rows))

    if junction_rows:
        columns = list(junction_rows[0].keys())
        out_path = os.path.join(OUT_DIR, 'splice_junction_mapping.tsv')

        with open(out_path, 'w') as f:
            f.write('\t'.join(columns) + '\n')

            for row in junction_rows:
                for c in columns:
                    _converted_3.append(str(row[c]))
                f.write('\t'.join(_converted_3) + '\n')
        print('wrote', out_path)
    report_path = os.path.join(OUT_DIR, 'transcript_annotation_report.md')
    f = open(report_path, 'w')
    f.write('# Transcript feature annotation report\n\n')
    f.write('## Method\n\n')
    f.write("Motif matches in transcript coordinates were mapped back to genomic coordinates using the exon structure from GENCODE v44. Each match was classified as 5' UTR, CDS, 3' UTR, exon, or intronic/other.\n\n")
    feature_counts = {}

    for row in annotated:
        feat = row.get('Feature', 'unknown')
        feature_counts[feat] = feature_counts.get(feat, 0) + 1
    f.write('## Feature distribution\n\n')
    f.write('| Feature | Count |\n|---------|-------|\n')

    for feat in sorted(feature_counts):
        f.write('| %s | %d |\n' % (feat, feature_counts[feat]))

    for r in junction_rows:
        if r['Junction_proximity'] != 'none':
            _converted_4.append(1)
    n_near = sum(_converted_4)
    f.write('\n## Splice junction proximity\n\n')
    f.write('- Matches near a splice junction (within 50 bp): %d\n' % n_near)
    f.write('- Matches not near a junction: %d\n' % (len(junction_rows) - n_near))
    f.write('\n## Caveats\n\n')
    f.write('- The GTF coordinates are 1-based and inclusive; the mapping accounts for strand.\n')
    f.write("- 'Exon' classification means the match falls in an exon that is not 5' UTR, CDS or 3' UTR in this transcript.\n")
    f.close()
    print('wrote', report_path)

if __name__ == '__main__':
    main()
