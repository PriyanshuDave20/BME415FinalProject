import sys

import os

import shutil

TRANSCRIPT_DIR = 'reference/target_transcripts'

OUT_DIR = 'results/target_isoform_docking/transcripts'

REQUIRED_FILES = {'NRXN1_transcripts.fa': 'NRXN1', 'NRXN2_transcripts.fa': 'NRXN2', 'NRXN3_transcripts.fa': 'NRXN3', 'AR_transcripts.fa': 'AR_V7'}

AR_V7_SOURCE_FILENAME = 'AR_transcripts.fa'

AR_V7_OUTPUT_FILENAME = 'AR_V7_transcript.fa'

REQUIRED_SINGLE_FILES = {'AR_transcripts.fa': 'AR_V7'}

def parse_header_gene(header):
    for i in header.split():
        if i.startswith('gene='):
            return i.split('=', 1)[1]
    return ''

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    problems = []

    for filename, gene in REQUIRED_FILES.items():
        src = os.path.join(TRANSCRIPT_DIR, filename)
        dst = os.path.join(OUT_DIR, filename)

        if not os.path.exists(src):
            print('  MISSING:', src)
            problems.append(filename + ' not found')
            continue

        if not os.path.exists(dst):
            shutil.copy(src, dst)
            print('  copied', src, '->', dst)
        else:
            print('  already exists:', dst)
    ar_v7_src = os.path.join(TRANSCRIPT_DIR, AR_V7_SOURCE_FILENAME)
    ar_v7_dst = os.path.join(OUT_DIR, AR_V7_OUTPUT_FILENAME)

    if os.path.exists(ar_v7_src):
        arv7_written = False

        if not os.path.exists(ar_v7_dst):
            with open(ar_v7_src) as inf, open(ar_v7_dst, 'w') as outf:
                header = None
                parts = []

                for line in inf:
                    line = line.strip()

                    if line.startswith('>'):
                        if header is not None:
                            gene_tag = parse_header_gene(header)
                            tid = header.split()[0]

                            if gene_tag == 'AR_V7' or tid == 'ENST00000504326.5':
                                outf.write('>' + header + '\n')

                                for i in range(0, len(''.join(parts)), 60):
                                    outf.write(''.join(parts)[i:i + 60] + '\n')
                                arv7_written = True
                        header = line[1:]
                        parts = []

                    elif line:
                        parts.append(line)

                if header is not None:
                    gene_tag = parse_header_gene(header)
                    tid = header.split()[0]

                    if gene_tag == 'AR_V7' or tid == 'ENST00000504326.5':
                        outf.write('>' + header + '\n')

                        for i in range(0, len(''.join(parts)), 60):
                            outf.write(''.join(parts)[i:i + 60] + '\n')
                        arv7_written = True
            print('  wrote AR-V7 transcript:', ar_v7_dst)
        else:
            print('  already exists:', ar_v7_dst)
    else:
        print('  MISSING:', ar_v7_src)
        problems.append(AR_V7_SOURCE_FILENAME + ' not found')

    if problems:
        print('\nPROBLEMS:')

        for p in problems:
            print('  -', p)
        print('Error: SystemExit(1)')
        sys.exit(1)
    print('\nAll transcript files copied/validated.')

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

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    problems = []

    for filename, gene in REQUIRED_FILES.items():
        src = os.path.join(TRANSCRIPT_DIR, filename)
        dst = os.path.join(OUT_DIR, filename)

        if not os.path.exists(src):
            print('  MISSING:', src)
            problems.append(filename + ' not found')
            continue

        if not os.path.exists(dst):
            shutil.copy(src, dst)
            print('  copied', src, '->', dst)
        else:
            print('  already exists:', dst)

    if problems:
        print('\nPROBLEMS:')

        for p in problems:
            print('  -', p)
        print('Error: SystemExit(1)')
        sys.exit(1)
    print('\nAll transcript files copied/validated.')

if __name__ == '__main__':
    main()
