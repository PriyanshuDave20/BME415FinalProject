import csv

import os

import subprocess

import sys

sample_file = 'results/sample_summary.tsv'

read_folder = 'results/preprocessed'

output_folder = 'results/alignments'

summary_file = 'results/alignment_summary.tsv'

reference = 'reference/GRCh38.primary_assembly.genome.fa.gz'

threads = 8

minimap2 = '/root/miniconda3/envs/bme415/bin/minimap2'

samtools = '/root/miniconda3/envs/bme415/bin/samtools'

def wsl_path(path):
    path = os.path.abspath(path)
    drive, rest = os.path.splitdrive(path)
    return '/mnt/' + drive[0].lower() + rest.replace('\\', '/')

def wsl_command(program, arguments):
    return ['wsl.exe', '-d', 'Ubuntu', '-u', 'root', '--', program] + arguments

def run(command):
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    text = result.stdout.decode(errors='replace') + result.stderr.decode(errors='replace')
    return (result.returncode, text)

def load_samples():
    if not os.path.exists(sample_file):
        print('Run load_prepare.py first.')
        sys.exit(1)

    with open(sample_file, newline='') as f:
        samples = list(csv.DictReader(f, delimiter='\t'))

    for sample in samples:
        if sample['status'] == 'OK':
            _converted_1.append(sample)
    samples = _converted_1

    if not samples:
        print('No valid samples were found.')
        sys.exit(1)
    return samples

def check_tools():
    code, version = run(wsl_command(minimap2, ['--version']))

    if code != 0:
        print('minimap2 is not available in WSL.')
        sys.exit(1)
    code, samtools_version = run(wsl_command(samtools, ['--version']))

    if code != 0:
        print('samtools is not available in WSL.')
        sys.exit(1)
    print('Using minimap2', version.splitlines()[0])
    print('Using', samtools_version.splitlines()[0])

def alignment_stats(bam):
    code, output = run(wsl_command(samtools, ['flagstat', wsl_path(bam)]))

    if code != 0:
        return (0, '', 0, 0)
    mapped = 0
    mapped_percent = ''
    paired = 0
    secondary = 0

    for i in output.splitlines():
        if ' mapped (' in i and 'not mapped' not in i:
            mapped = int(i.split('+')[0])
            mapped_percent = i.split('(')[1].split('%')[0]

        elif 'properly paired' in i:
            paired = int(i.split('+')[0])

        elif 'secondary' in i:
            secondary = int(i.split('+')[0])
    return (mapped, mapped_percent, paired, secondary)

def align_sample(sample, dry_run):
    name = sample['sample_id']
    r1 = os.path.join(read_folder, name + '_R1.fastq.gz')
    r2 = os.path.join(read_folder, name + '_R2.fastq.gz')

    if not os.path.exists(r1) or not os.path.exists(r2):
        return (None, 'Preprocessed FASTQs are missing. Run preprocess_reads.py first.')
    folder = os.path.join(output_folder, name)
    bam = os.path.join(folder, name + '_aligned.bam')
    bai = bam + '.bai'

    if os.path.exists(bam) and os.path.exists(bai):
        return ((bam, bai), 'already exists')

    if dry_run:
        return ((bam, bai), 'would align')
    os.makedirs(folder, exist_ok=True)
    log_file = os.path.join(folder, 'alignment.log')
    align = wsl_command(minimap2, ['-ax', 'sr', '-t', str(threads), '-I', '2G', wsl_path(reference), wsl_path(r1), wsl_path(r2)])
    sort = wsl_command(samtools, ['sort', '-@', '2', '-o', wsl_path(bam), '-'])

    with open(log_file, 'w') as log:
        first = subprocess.Popen(align, stdout=subprocess.PIPE, stderr=log)
        second = subprocess.Popen(sort, stdin=first.stdout, stdout=subprocess.DEVNULL, stderr=log)
        first.stdout.close()
        second.wait()
        first.wait()

    if first.returncode != 0 or second.returncode != 0:
        return (None, 'Alignment failed. Read ' + log_file)
    code, _ = run(wsl_command(samtools, ['index', wsl_path(bam)]))

    if code != 0:
        return (None, 'samtools index failed')
    return ((bam, bai), '')

def main():
    dry_run = '--dry-run' in sys.argv
    print('Module 2: align preprocessed reads')

    if not os.path.exists(reference):
        print('Reference genome not found:', reference)
        sys.exit(1)
    check_tools()
    samples = load_samples()
    os.makedirs(output_folder, exist_ok=True)
    rows = []
    failures = []

    for sample in samples:
        print('Aligning', sample['sample_id'])
        files, note = align_sample(sample, dry_run)

        if files is None:
            print('  ERROR:', note)
            failures.append(sample['sample_id'])
            rows.append([sample['sample_id'], sample['condition'], sample['replicate'], 0, '', 0, 0, '', '', 'FAILED', note])
            continue
        bam, bai = files

        if dry_run:
            print('  ', note)
            continue
        mapped, percent, paired, secondary = alignment_stats(bam)
        rows.append([sample['sample_id'], sample['condition'], sample['replicate'], mapped, percent, paired, secondary, bam, bai, 'OK', note])

    if dry_run:
        print('Dry run finished. No files were changed.')
        return

    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['sample_id', 'condition', 'replicate', 'mapped_reads', 'mapped_pct', 'properly_paired_reads', 'multimapped_reads', 'bam_file', 'bai_file', 'status', 'notes'])
        writer.writerows(rows)
    print('Wrote', summary_file)

    if failures:
        print('Failed samples:', ', '.join(failures))
        sys.exit(1)

if __name__ == '__main__':
    main()
