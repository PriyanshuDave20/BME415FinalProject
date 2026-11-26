import os

gtf_file = 'reference/gencode.v44.annotation.gtf'

fasta_file = 'reference/GRCh38.primary_assembly.genome.fa'

out_folder = 'reference/target_transcripts'

target_genes = {'ENSG00000179915': 'NRXN1', 'ENSG00000110076': 'NRXN2', 'ENSG00000021645': 'NRXN3', 'ENSG00000169083': 'AR'}

def parse_attributes(text):
    attrs = {}

    for i in text.strip().split(';'):
        i = i.strip()

        if i == '':
            continue
        pieces = i.split(' ', 1)

        if len(pieces) == 2:
            attrs[pieces[0]] = pieces[1].strip('"')
    return attrs

def read_gtf_exons():
    transcripts = {}
    gtf = open(gtf_file)

    for line in gtf:
        if line.startswith('#'):
            continue
        columns = line.rstrip('\n').split('\t')

        if len(columns) < 9:
            continue

        if columns[2] != 'exon':
            continue
        attrs = parse_attributes(columns[8])
        gene_id = attrs.get('gene_id', '').split('.')[0]

        if gene_id not in target_genes:
            continue
        transcript_id = attrs.get('transcript_id', '')

        if transcript_id == '':
            continue

        if transcript_id not in transcripts:
            transcripts[transcript_id] = {'gene_id': gene_id, 'gene_name': attrs.get('gene_name', ''), 'transcript_type': attrs.get('transcript_type', ''), 'chrom': columns[0], 'strand': columns[6], 'exons': []}
        start = int(columns[3])
        end = int(columns[4])
        transcripts[transcript_id]['exons'].append((start, end))
    gtf.close()
    return transcripts

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    out = []

    for i in reversed(seq):
        out.append(complement[i])
    return ''.join(out)

def load_chromosomes(chroms_needed):
    chromosomes = {}
    current_name = None
    current_parts = []
    f = open(fasta_file)

    for line in f:
        if line.startswith('>'):
            if current_name in chroms_needed:
                chromosomes[current_name] = ''.join(current_parts)
            current_name = line[1:].split()[0]
            current_parts = []

        elif current_name in chroms_needed:
            current_parts.append(line.strip())

    if current_name in chroms_needed:
        chromosomes[current_name] = ''.join(current_parts)
    f.close()
    return chromosomes

def build_mrna(transcript, chromosomes):
    chrom = transcript['chrom']

    if chrom not in chromosomes:
        return (None, 'chromosome ' + chrom + ' not found in the genome file')
    chrom_seq = chromosomes[chrom]
    exons = sorted(transcript['exons'])
    pieces = []

    for start, end in exons:
        if start < 1 or start > end or end > len(chrom_seq):
            return (None, 'invalid exon coordinates %d-%d' % (start, end))
        pieces.append(chrom_seq[start - 1:end])
    mrna = ''.join(pieces)

    if transcript['strand'] == '-':
        mrna = reverse_complement(mrna)
    return (mrna, None)

def check_ar_v7(ar_transcripts, chromosomes):
    candidates = []

    for tid in sorted(ar_transcripts):
        t = ar_transcripts[tid]

        if t['transcript_type'] != 'protein_coding':
            continue

        if t['strand'] != '+':
            continue
        exons = sorted(t['exons'])

        if len(exons) < 5:
            continue
        exon3_end = exons[2][1]
        exon4_start = exons[3][0]
        intron_start = exon3_end + 1
        intron_end = exon4_start - 1

        if intron_end - intron_start < 117:
            continue

        for other_tid in sorted(ar_transcripts):
            for start, end in ar_transcripts[other_tid]['exons']:
                if start >= intron_start and end <= intron_end:
                    length = end - start + 1
                    matches_ce3 = length == 117
                    candidates.append((tid, other_tid, intron_start, intron_end, start, end, length, matches_ce3))
    return candidates

def main():
    if not os.path.exists(gtf_file):
        print('ERROR: GTF file not found:', gtf_file)
        return

    if not os.path.exists(fasta_file):
        print('ERROR: genome file not found:', fasta_file)
        return
    os.makedirs(out_folder, exist_ok=True)
    print('Reading the GTF annotation ...')
    transcripts = read_gtf_exons()
    print('  found %d transcripts for the target genes' % len(transcripts))
    chroms_needed = set()

    for tid in transcripts:
        chroms_needed.add(transcripts[tid]['chrom'])
    print('  chromosomes needed:', sorted(chroms_needed))
    print('Reading the genome file (this takes a few minutes) ...')
    chromosomes = load_chromosomes(chroms_needed)
    print('Building transcript sequences ...')
    gene_counts = {}
    gene_fasta = {}
    problems = []

    for tid in sorted(transcripts):
        t = transcripts[tid]
        gene_name = t['gene_name']
        mrna, error = build_mrna(t, chromosomes)

        if error is not None:
            problems.append(tid + ': ' + error)
            print('  WARNING:', tid, error)
            continue

        if gene_name not in gene_fasta:
            gene_fasta[gene_name] = []
        gene_fasta[gene_name].append((tid, t, mrna))
        gene_counts[gene_name] = gene_counts.get(gene_name, 0) + 1

    for gene_name in sorted(gene_fasta):
        out_path = os.path.join(out_folder, gene_name + '_transcripts.fa')
        f = open(out_path, 'w')

        for tid, t, mrna in gene_fasta[gene_name]:
            header = '>%s gene=%s type=%s chrom=%s:%s length=%d' % (tid, gene_name, t['transcript_type'], t['chrom'], t['strand'], len(mrna))
            f.write(header + '\n')

            for i in range(0, len(mrna), 60):
                f.write(mrna[i:i + 60] + '\n')
        f.close()
        print('  wrote', out_path, '(%d transcripts)' % len(gene_fasta[gene_name]))
    print('\nChecking for an AR-V7 transcript in the annotation ...')
    ar_transcripts = {}

    for tid, t, mrna in gene_fasta.get('AR', []):
        ar_transcripts[tid] = {'transcript_type': t['transcript_type'], 'strand': t['strand'], 'exons': t['exons']}
    candidates = check_ar_v7(ar_transcripts, chromosomes)

    if len(candidates) == 0:
        print('  No AR-V7 transcript was found in the GENCODE v44 annotation.')
        print('  None of the annotated AR transcripts contains an exon inside')
        print('  intron 3 (where the AR-V7 cryptic exon is located in cancer')
        print('  cells). The cryptic exon is not part of the normal genome')
        print('  annotation, so AR_V7.fa was NOT written.')
    else:
        ce3_matches = []

        for tid, other_tid, istart, iend, start, end, length, matches_ce3 in candidates:
            print('    transcript %s has an exon at %d-%d (%d bp) inside' % (other_tid, start, end, length))
            print('    intron 3 (%d-%d) of %s' % (istart, iend, tid))

            if matches_ce3:
                ce3_matches.append((other_tid, start, end, length))

        if len(ce3_matches) == 0:
            print('  None of these exons has the published 117 bp length of')
            print('  the AR-V7 cryptic exon, so none can be confidently')
            print('  identified as AR-V7. AR_V7.fa was NOT written.')
        else:
            v7_path = os.path.join(out_folder, 'AR_V7.fa')
            f = open(v7_path, 'w')

            for other_tid, start, end, length in ce3_matches:
                f.write('>%s AR-V7_candidate length=%d\n' % (other_tid, length))
            f.close()
            print('  wrote', v7_path)
    print('\nSummary')
    print('  transcripts extracted:')

    for gene_name in sorted(gene_counts):
        print('    %s: %d transcripts' % (gene_name, gene_counts[gene_name]))

    if len(problems) > 0:
        print('  problems:')

        for i in problems:
            print('    ' + i)

if __name__ == '__main__':
    main()
