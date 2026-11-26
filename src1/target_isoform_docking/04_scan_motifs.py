import os

import re

from collections import defaultdict

SELECTED_FA = 'results/target_isoform_docking/selected_isoforms/selected_isoforms.fa'

OUT_DIR = 'results/target_isoform_docking/motif_scanning'

WT5_CONSENSUS = 'CCGTCTCAA'

R5_CONSENSUS = 'GGTCTYGMWC'

IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}

def build_regex(consensus):
    pattern = []

    for i in consensus.upper():
        if i in IUPAC:
            alts = IUPAC[i]

            if len(alts) == 1:
                pattern.append(alts)
            else:
                pattern.append('[' + alts + ']')
        else:
            pattern.append('.')
    return re.compile(''.join(pattern))

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
    fields = {'transcript_id': '', 'gene': '', 'type': '', 'chrom': '', 'strand': '', 'length': 0}
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

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    wt5_regex = build_regex(WT5_CONSENSUS)
    r5_regex = build_regex(R5_CONSENSUS)
    all_matches = []
    isoform_match_counts = defaultdict(lambda: {'WT5': 0, 'R5': 0})
    gene_match_counts = defaultdict(lambda: {'WT5': 0, 'R5': 0})
    isoform_seqs = {}

    for header, seq in parse_fasta(SELECTED_FA):
        info = parse_header(header)
        tid = info['transcript_id']
        gene = info['gene']
        upper = seq.upper()

        for motif_name, consensus, regex in (('WT5', WT5_CONSENSUS, wt5_regex), ('R5', R5_CONSENSUS, r5_regex)):
            for j in regex.finditer(upper):
                start = j.start() + 1
                end = j.end()
                matched = j.group()
                all_matches.append({'Gene': gene, 'Transcript_ID': tid, 'Transcript_version': tid, 'Motif_name': motif_name, 'Expected_consensus': consensus, 'Matched_sequence': matched, 'Transcript_start': start, 'Transcript_end': end, 'Transcript_length': len(upper), 'Orientation': '+'})
                isoform_match_counts[tid][motif_name] += 1
                gene_match_counts[gene][motif_name] += 1
                key = (gene, motif_name, matched)

                if key not in isoform_seqs:
                    isoform_seqs[key] = set()
                isoform_seqs[key].add(tid)

    for row in all_matches:
        key = (row['Gene'], row['Motif_name'], row['Matched_sequence'])
        row['N_isoforms_sharing'] = len(isoform_seqs.get(key, set()))
    out1 = os.path.join(OUT_DIR, 'all_motif_matches.tsv')
    cols = list(all_matches[0].keys()) if all_matches else []

    with open(out1, 'w') as f:
        if cols:
            f.write('\t'.join(cols) + '\n')

            for row in all_matches:
                for c in cols:
                    _converted_1.append(str(row[c]))
                f.write('\t'.join(_converted_1) + '\n')
        else:
            f.write('no_matches\n')
    print('wrote', out1, '(%d matches)' % len(all_matches))
    out2 = os.path.join(OUT_DIR, 'motif_summary_by_gene.tsv')

    with open(out2, 'w') as f:
        f.write('Gene\tWT5_matches\tR5_matches\tTotal\n')

        for gene in ('NRXN1', 'NRXN2', 'NRXN3', 'AR_V7'):
            w = gene_match_counts[gene]['WT5']
            r = gene_match_counts[gene]['R5']
            f.write('%s\t%d\t%d\t%d\n' % (gene, w, r, w + r))
    print('wrote', out2)
    out3 = os.path.join(OUT_DIR, 'motif_summary_by_isoform.tsv')

    with open(out3, 'w') as f:
        f.write('Transcript_ID\tGene\tWT5_matches\tR5_matches\tTotal\n')

        for tid, counts in sorted(isoform_match_counts.items()):
            w = counts['WT5']
            r = counts['R5']
            f.write('%s\t\t%d\t%d\t%d\n' % (tid, w, r, w + r))
    print('wrote', out3)
    out4 = os.path.join(OUT_DIR, 'motif_negative_isoforms.tsv')

    with open(out4, 'w') as f:
        f.write('Transcript_ID\tGene\n')

        for header, seq in parse_fasta(SELECTED_FA):
            info = parse_header(header)
            tid = info['transcript_id']
            counts = isoform_match_counts.get(tid, {'WT5': 0, 'R5': 0})

            if counts['WT5'] == 0 and counts['R5'] == 0:
                f.write('%s\t%s\n' % (tid, info['gene']))
    print('wrote', out4)

if __name__ == '__main__':
    main()
