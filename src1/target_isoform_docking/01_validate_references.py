import sys

import os

GTF = 'reference/gencode.v44.annotation.gtf'

FASTA = 'reference/GRCh38.primary_assembly.genome.fa'

TRANSCRIPT_DIR = 'reference/target_transcripts'

OUT_DIR = 'results/target_isoform_docking/reports'

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    problems = []
    lines = []
    lines.append('# Reference validation report\n')

    for i in (GTF, FASTA):
        label = 'GTF' if i == GTF else 'FASTA'

        if os.path.exists(i):
            size = os.path.getsize(i)
            lines.append('- `%s`: **OK** (%d bytes)\n' % (i, size))
        else:
            lines.append('- `%s`: **MISSING**\n' % i)
            problems.append(label + ' file missing: ' + i)

    if os.path.exists(TRANSCRIPT_DIR):
        for f in sorted(os.listdir(TRANSCRIPT_DIR)):
            if f.endswith('.fa'):
                _converted_1.append(f)
        fa_files = _converted_1
        lines.append('- Transcript directory: **OK** (%d FASTA files)\n' % len(fa_files))

        for f in fa_files:
            lines.append('  - `%s`\n' % f)
    else:
        lines.append('- Transcript directory: **MISSING**\n')
        problems.append('Transcript directory missing: ' + TRANSCRIPT_DIR)
    lines.append('\n## Summary\n\n')

    if problems:
        for p in problems:
            lines.append('- PROBLEM: %s\n' % p)
    else:
        lines.append('All references validated OK.\n')
    out = os.path.join(OUT_DIR, 'reference_validation.md')

    with open(out, 'w') as f:
        f.write(''.join(lines))
    print('wrote', out)

    for p in problems:
        print('  PROBLEM:', p)

    if problems:
        print('Error: SystemExit(1)')
        sys.exit(1)

if __name__ == '__main__':
    main()
