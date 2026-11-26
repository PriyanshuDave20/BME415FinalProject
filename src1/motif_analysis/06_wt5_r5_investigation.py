import csv

import math

import os

import re

import statistics

import subprocess

import sys

from bisect import bisect_right

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

out_folder = 'results/wt5_r5_investigation'

motif_folder = os.path.join(out_folder, 'motif_files')

fimo_folder = os.path.join(out_folder, 'fimo')

figures_folder = os.path.join(out_folder, 'figures')

motif_folder_module5 = 'results/motifs'

fasta_folder = 'results/motifs/fasta'

wt_fasta = 'results/motifs/wt_sequences.fasta'

r168c_fasta = 'results/motifs/r168c_sequences.fasta'

ds_wt_fasta = 'results/motifs/fasta/ds_wt_sequences.fasta'

ds_r168c_fasta = 'results/motifs/fasta/ds_r168c_sequences.fasta'

overlap_file = 'results/comparison/region_overlap.tsv'

gtf_file = 'reference/gencode.v44.annotation.gtf'

wt5_consensus = 'CCGTCTCAA'

r5_consensus = 'GGTCTYGMWC'

fimo_pvalue_threshold = '1e-4'

streme_wt = 'wsl.exe -d Ubuntu -u root -- /root/miniconda3/envs/bme415/bin/streme'

tomtom_exe = '/root/miniconda3/envs/bme415/bin/tomtom'

fimo_exe = '/root/miniconda3/envs/bme415/bin/fimo'

samtools_exe = '/root/miniconda3/envs/bme415/bin/samtools'

wt_tag = 'streme_wt'

r_tag = 'streme_r168c'

def wslp(p):
    p = os.path.abspath(p)
    drive, rest = os.path.splitdrive(p)
    return '/mnt/' + drive[0].lower() + rest.replace('\\', '/')

def wsl_cmd(exe, args):
    return ['wsl.exe', '-d', 'Ubuntu', '-u', 'root', '--', exe] + args

def parse_streme_motifs(path):
    motifs = []
    cur = None
    f = open(path)

    for line in f:
        if line.startswith('MOTIF'):
            x = line.split()
            con = ''
            mid = x[-1] if len(x) > 1 else 'motif'

            if len(x) >= 3:
                con = x[1]
                mid = x[2]

            if '-' in con and con.split('-')[0].isdigit():
                con = con.split('-', 1)[1]
            cur = {'id': mid, 'con': con, 'w': 0, 'nsites': 0, 'e': None, 'pwm': []}
            motifs.append(cur)

        elif cur is not None and line.startswith('letter-probability'):
            import re
            mw = re.search('w=\\s*(\\d+)', line)
            mn = re.search('nsites=\\s*(\\d+)', line)
            me = re.search('E=\\s*(\\S+)', line)

            if mw:
                cur['w'] = int(mw.group(1))

            if mn:
                cur['nsites'] = int(mn.group(1))

            if me:
                cur['e'] = float(me.group(1))

        elif cur is not None:
            x = line.split()

            if len(x) == 4:
                try:
                    for v in x:
                        _converted_1.append(float(v))
                    cur['pwm'].append(_converted_1)
                except ValueError:
                    continue

                if len(cur['pwm']) == cur['w']:
                    cur = None
    f.close()
    return motifs

def write_meme_file(path, motifs):
    f = open(path, 'w')
    f.write('MEME version 4\n\n')
    f.write('ALPHABET= ACGT\n\n')
    f.write('strands: + -\n\n')
    f.write('Background letter frequencies (from uniform background):\n')
    f.write('A 0.25 C 0.25 G 0.25 T 0.25\n\n')

    for m in motifs:
        f.write('MOTIF %s %s\n' % (m['con'], m['id']))
        f.write('letter-probability matrix: alength= 4 w= %d nsites= %d E= %s\n' % (m['w'], m['nsites'], m['e']))

        for row in m['pwm']:
            f.write('  %.6f  %.6f  %.6f  %.6f\n' % tuple(row))
        f.write('\n')
    f.close()
    print('  wrote', path)

def find_motif(motifs, consensus):
    for m in motifs:
        if m['con'] == consensus:
            return m
    return None

def run_tomtom(query_file, target_file, outdir):
    if os.path.exists(os.path.join(outdir, 'tomtom.tsv')):
        print('  tomtom already done, skipping')
        return
    os.makedirs(outdir, exist_ok=True)
    cmd = wsl_cmd(tomtom_exe, ['--oc', wslp(outdir), '-thresh', '10', '-evalue', wslp(query_file), wslp(target_file)])

    with open(os.path.join(outdir, 'tomtom_stderr.log'), 'w') as errf:
        result = subprocess.run(cmd, stdout=errf, stderr=errf)

    if result.returncode != 0:
        print("Error: RuntimeError('tomtom failed: ' + open(os.path.join(outdir, 'tomtom_stderr.log')).read()[-400:])")
        sys.exit(1)

def parse_tomtom(outdir):
    rows = []
    path = os.path.join(outdir, 'tomtom.tsv')
    f = open(path)

    for line in f:
        if line.startswith('Query_ID') or line.startswith('#'):
            continue
        x = line.rstrip('\n').split('\t')

        if len(x) >= 9:
            rows.append({'query': x[0], 'target': x[1], 'offset': int(x[2]), 'p': float(x[3]), 'e': float(x[4]), 'q': float(x[5]), 'overlap': int(x[6]), 'orientation': x[8]})
    f.close()
    return rows

def build_shared_fasta():
    coords = {}

    for tag, fasta in (('WT', wt_fasta), ('R168C', r168c_fasta)):
        f = open(fasta)

        for line in f:
            if line.startswith('>'):
                rid = line[1:].strip()
                chrom, span = rid.split('|')[1].split(':')
                start, end = span.split('-')
                coords[tag, rid] = (chrom, int(start), int(end))
        f.close()
    pairs = []
    f = open(os.path.join('results/comparison', 'region_overlap.tsv'), newline='')

    for r in csv.DictReader(f, delimiter='\t'):
        if r['category'] != 'shared' or not r['partner_id']:
            continue

        if not r['region_id'].startswith('WT_'):
            continue
        wk = ('WT', r['region_id'])
        rk = ('R168C', r['partner_id'])

        if wk not in coords or rk not in coords:
            continue
        wc = coords[wk]
        rc = coords[rk]

        if wc[0] != rc[0]:
            continue
        ov_start = max(wc[1], rc[1])
        ov_end = min(wc[2], rc[2])

        if ov_end - ov_start < 20:
            continue
        pairs.append((r['region_id'], r['partner_id'], wc[0], ov_start, ov_end))
    f.close()
    import gzip
    by_chrom = {}

    for wt_id, r_id, chrom, start, end in pairs:
        by_chrom.setdefault(chrom, []).append((wt_id, r_id, start, end))
    seqs = []
    found = {}
    cur = None
    parts = []
    f = gzip.open(ref, 'rt')

    def flush(chrom, parts):
        if chrom is None or chrom not in by_chrom:
            return
        seq = ''.join(parts)

        for wt_id, r_id, start, end in by_chrom[chrom]:
            centre = start + (end - start) // 2
            half = 50
            s = seq[centre - half:centre - half + 100].upper()

            if len(s) < 100 or 'N' in s:
                continue
            key = (chrom, start, end)

            if key in found:
                continue
            found[key] = True
            seqs.append((wt_id, r_id, chrom, start, end, s))

    for line in f:
        if line.startswith('>'):
            flush(cur, parts)
            cur = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    flush(cur, parts)
    f.close()
    path = os.path.join(out_folder, 'shared_regions.fasta')
    f = open(path, 'w')

    for wt_id, r_id, chrom, start, end, s in seqs:
        f.write('>%s|%s|%s:%d-%d\n' % (wt_id, r_id, chrom, start, end))

        for i in range(0, len(s), 60):
            f.write(s[i:i + 60] + '\n')
    f.close()
    print('  wrote %s (%d shared windows)' % (path, len(seqs)))
    return (path, len(seqs))

def run_fimo(tag, motif_file, fasta_file):
    outdir = os.path.join(fimo_folder, tag)

    if os.path.exists(os.path.join(outdir, 'fimo.tsv')):
        print('  fimo %s: already done, skipping' % tag)
        return
    os.makedirs(outdir, exist_ok=True)
    cmd = wsl_cmd(fimo_exe, ['--oc', wslp(outdir), '--max-stored-matches', '1000000', '--thresh', fimo_pvalue_threshold, wslp(motif_file), wslp(fasta_file)])

    with open(os.path.join(outdir, 'fimo_stderr.log'), 'w') as errf:
        result = subprocess.run(cmd, stdout=errf, stderr=errf)

    if result.returncode != 0:
        print("Error: RuntimeError('fimo failed for %s: %s' % (tag, open(os.path.join(outdir, 'fimo_stderr.log')).read()[-300:]))")
        sys.exit(1)
    print('  fimo %s: done' % tag)

def parse_fimo(tag):
    path = os.path.join(fimo_folder, tag, 'fimo.tsv')
    n_matches = 0
    seq_names = set()
    strands = {'+': 0, '-': 0}
    f = open(path)

    for line in f:
        if line.startswith('#') or line.startswith('motif_id'):
            continue
        x = line.rstrip('\n').split('\t')

        if len(x) >= 8:
            n_matches += 1
            seq_names.add(x[2])

            if x[7] == '+':
                strands['+'] += 1
            else:
                strands['-'] += 1
    f.close()
    return (n_matches, seq_names, strands)

def count_fasta(path):
    return countfa(path)

def countfa(path):
    n = 0
    f = open(path)

    for line in f:
        if line.startswith('>'):
            n = n + 1
    f.close()
    return n

def fisher_and_or(with_a, total_a, with_b, total_b):
    from scipy.stats import fisher_exact
    a = with_a
    b = total_a - with_a
    c = with_b
    d = total_b - with_b
    odds, p = fisher_exact([[a, b], [c, d]])

    if with_b == 0 or total_b == 0:
        enrichment = float('inf') if with_a > 0 else float('nan')
    else:
        enrichment = with_a / total_a / (with_b / total_b)
    return (p, odds, enrichment)

def mcnemar_exact(both, only_wt5, only_r5, neither):
    from scipy.stats import binomtest
    discordant = only_wt5 + only_r5

    if discordant == 0:
        return 1.0
    return binomtest(only_wt5, discordant, 0.5).pvalue

def parse_gtf_lite(path):
    exons = {}
    genes = {}
    utr5 = {}
    utr3 = {}
    cds = {}
    f = open(path)

    for line in f:
        if line.startswith('#'):
            continue
        x = line.rstrip('\n').split('\t')

        if len(x) < 9:
            continue
        feature = x[2]

        if feature not in ('exon', 'transcript', 'five_prime_utr', 'three_prime_utr', 'CDS'):
            continue
        chrom = x[0]
        start = int(x[3]) - 1
        end = int(x[4])
        strand = x[6]

        if feature == 'exon':
            exons.setdefault(chrom, []).append((start, end, strand))

        elif feature == 'transcript':
            genes.setdefault(chrom, []).append((start, end, strand))

        elif feature == 'five_prime_utr':
            utr5.setdefault(chrom, []).append((start, end, strand))

        elif feature == 'three_prime_utr':
            utr3.setdefault(chrom, []).append((start, end, strand))

        elif feature == 'CDS':
            cds.setdefault(chrom, []).append((start, end, strand))
    f.close()

    def sortd(d):
        for chrom in d:
            d[chrom].sort()
        return d
    return (sortd(exons), sortd(genes), sortd(utr5), sortd(utr3), sortd(cds))

def build_index(intervals):
    for s, e, st in intervals:
        _converted_2.append(s)
    starts = _converted_2
    return starts

def overlaps(intervals, starts, chrom, start, end):
    if chrom not in intervals:
        return False
    idx = bisect_right(starts[chrom], start)
    candidates = intervals[chrom][max(0, idx - 1):idx + 1]

    for s, e, st in candidates:
        if s < end and start < e:
            return True
    i = idx
    ints = intervals[chrom]

    while i < len(ints) and ints[i][0] < end:
        s, e = (ints[i][0], ints[i][1])

        if s < end and start < e:
            return True
        i += 1
    return False

def classify_occurrence(chrom, start, end, strand, exons, genes, utr5, utr3, cds):
    if chrom not in genes:
        return ('intergenic', '')
    exon_hits = overlapping_features(exons, chrom, start, end)

    if exon_hits:
        if overlapping_features(utr5, chrom, start, end):
            label = "5' UTR"

        elif overlapping_features(utr3, chrom, start, end):
            label = "3' UTR"

        elif overlapping_features(cds, chrom, start, end):
            label = 'CDS'
        else:
            label = 'exon'
        gene_strand = exon_hits[0][2]
        sense = 'sense' if gene_strand == strand else 'antisense'
        return (label, sense)
    gene_hits = overlapping_features(genes, chrom, start, end)

    if gene_hits:
        gene_strand = gene_hits[0][2]
        sense = 'sense' if gene_strand == strand else 'antisense'
        return ('intron', sense)
    return ('intergenic', '')

def overlapping_features(feature_dict, chrom, start, end):
    if chrom not in feature_dict:
        return []
    ints = feature_dict[chrom]

    for s, e, st in ints:
        _converted_3.append(s)
    starts = _converted_3
    idx = bisect_right(starts, start)
    hits = []
    i = max(0, idx - 1)

    while i < len(ints) and ints[i][0] < end:
        s, e, st = ints[i]

        if s < end and start < e:
            hits.append(ints[i])
        i += 1
    return hits

def main():
    print('=' * 70)
    print('WT5 vs R5 motif investigation')
    print('=' * 70)

    for d in (out_folder, motif_folder, fimo_folder, figures_folder):
        os.makedirs(d, exist_ok=True)
    print('\nStep 1: locating WT5 and R5 PWMs ...')
    wt_motifs = parse_streme_motifs(os.path.join(motif_folder_module5, 'streme_wt', 'streme.txt'))
    r_motifs = parse_streme_motifs(os.path.join(motif_folder_module5, 'streme_r168c', 'streme.txt'))
    wt5 = find_motif(wt_motifs, wt5_consensus)
    r5 = find_motif(r_motifs, r5_consensus)

    if wt5 is None or r5 is None:
        print('ERROR: could not find WT5 (%s) or R5 (%s) in the STREME outputs' % (wt5_consensus, r5_consensus))
        sys.exit(1)
    print('  WT5: %s (id %s, width %d, E %.2e, %d sites)' % (wt5['con'], wt5['id'], wt5['w'], wt5['e'], wt5['nsites']))
    print('  R5:  %s (id %s, width %d, E %.2e, %d sites)' % (r5['con'], r5['id'], r5['w'], r5['e'], r5['nsites']))
    write_meme_file(os.path.join(motif_folder, 'wt5.meme'), [wt5])
    write_meme_file(os.path.join(motif_folder, 'r5.meme'), [r5])
    write_meme_file(os.path.join(motif_folder, 'wt5_and_r5.meme'), [wt5, r5])
    print('\nStep 2: TOMTOM WT5 vs R5 ...')
    tt_dir = os.path.join(out_folder, 'tomtom')
    run_tomtom(os.path.join(motif_folder, 'wt5_and_r5.meme'), os.path.join(motif_folder, 'wt5_and_r5.meme'), tt_dir)
    tomtom_rows = parse_tomtom(tt_dir)

    for r in tomtom_rows:
        if r['query'].startswith('WT5') and r['target'].startswith('R5') or (r['query'].startswith('R5') and r['target'].startswith('WT5')):
            _converted_4.append(r)
    tomtom_wt5_r5 = _converted_4
    print('  TOMTOM results:')

    for r in tomtom_wt5_r5:
        print('    %s vs %s: offset %d, overlap %d, orientation %s, p %.2e, E %.2e, q %.2e' % (r['query'], r['target'], r['offset'], r['overlap'], r['orientation'], r['p'], r['e'], r['q']))
    print('\nStep 3: building the shared-region sequence collection ...')
    shared_fa, n_shared = build_shared_fasta()
    print('\nStep 4: FIMO cross-scans ...')
    collections = [('wt', wt_fasta), ('r168c', r168c_fasta), ('shared', shared_fa), ('ds_wt', ds_wt_fasta), ('ds_r168c', ds_r168c_fasta)]
    brs = getbr()

    for nm in sorted(brs):
        collections.append((nm, os.path.join(fasta_folder, nm.replace('-', '_').lower() + '_sequences.fasta')))
    scans = {}

    for motif_tag, motif_file in (('WT5', os.path.join(motif_folder, 'wt5.meme')), ('R5', os.path.join(motif_folder, 'r5.meme'))):
        for coll_tag, fasta in collections:
            scan_tag = '%s_in_%s' % (motif_tag, coll_tag)
            run_fimo(scan_tag, motif_file, fasta)
            n_matches, seqs, strands = parse_fimo(scan_tag)
            n_total = count_fasta(fasta)
            scans[scan_tag] = {'n_matches': n_matches, 'seqs': seqs, 'n_total': n_total, 'strands': strands}
            print('  %-22s %5d matches in %5d of %6d sequences' % (scan_tag, n_matches, len(seqs), n_total))
    print('\nStep 5: occurrence summary and statistics ...')
    path = os.path.join(out_folder, 'occurrence_summary.tsv')
    f = open(path, 'w', newline='')
    w = csv.writer(f, delimiter='\t')
    w.writerow(['motif', 'collection', 'n_sequences', 'n_with_motif', 'pct_sequences_with_motif', 'n_matches', 'matches_per_1000_sequences', 'plus_strand_matches', 'minus_strand_matches'])

    for key in sorted(scans):
        motif_tag, coll_tag = key.split('_in_')
        s = scans[key]
        pct = 100.0 * len(s['seqs']) / s['n_total'] if s['n_total'] else 0.0
        per1000 = 1000.0 * s['n_matches'] / s['n_total'] if s['n_total'] else 0.0
        w.writerow([motif_tag, coll_tag, s['n_total'], len(s['seqs']), '%.2f' % pct, s['n_matches'], '%.2f' % per1000, s['strands']['+'], s['strands']['-']])
    f.close()
    print('  wrote', path)
    fisher_rows = []

    for motif_tag, con in (('WT5', wt5_consensus), ('R5', r5_consensus)):
        w_key = '%s_in_wt' % motif_tag
        r_key = '%s_in_r168c' % motif_tag
        with_w = len(scans[w_key]['seqs'])
        tot_w = scans[w_key]['n_total']
        with_r = len(scans[r_key]['seqs'])
        tot_r = scans[r_key]['n_total']
        p, odds, enr = fisher_and_or(with_w, tot_w, with_r, tot_r)
        fisher_rows.append([motif_tag, with_w, tot_w, with_r, tot_r, '%.3g' % p, '%.2f' % odds, '%.2f' % enr])
        print('  %s: WT %d/%d vs R168C %d/%d -> Fisher p %.3g, OR %.2f, enrichment %.2f' % (motif_tag, with_w, tot_w, with_r, tot_r, p, odds, enr))
    wt5_shared = scans['WT5_in_shared']['seqs']
    r5_shared = scans['R5_in_shared']['seqs']
    both = len(wt5_shared & r5_shared)
    only_wt5 = len(wt5_shared - r5_shared)
    only_r5 = len(r5_shared - wt5_shared)
    neither = n_shared - both - only_wt5 - only_r5
    p_mcnemar = mcnemar_exact(both, only_wt5, only_r5, neither)
    print('  shared collection (n = %d identical sequences):' % n_shared)
    print('    both motifs: %d, WT5 only: %d, R5 only: %d, neither: %d' % (both, only_wt5, only_r5, neither))
    print('    exact McNemar p = %.3g' % p_mcnemar)
    path = os.path.join(out_folder, 'fisher_tests.tsv')
    f = open(path, 'w', newline='')
    w = csv.writer(f, delimiter='\t')
    w.writerow(['test', 'motif', 'with_motif_condition1', 'total_condition1', 'with_motif_condition2', 'total_condition2', 'p_value', 'odds_ratio', 'enrichment_ratio_condition1_over_2'])

    for row in fisher_rows:
        w.writerow(['fisher_exact WT vs R168C'] + row)
    w.writerow(['mcnemar_exact shared collection', 'WT5 vs R5', both + only_wt5, n_shared, both + only_r5, n_shared, '%.3g' % p_mcnemar, '', ''])
    f.close()
    print('  wrote', path)
    print('\nStep 6: replicate consistency ...')
    br_rows = []

    for nm in sorted(brs):
        for motif_tag in ('WT5', 'R5'):
            key = '%s_in_%s' % (motif_tag, nm)
            s = scans[key]
            pct = 100.0 * len(s['seqs']) / s['n_total'] if s['n_total'] else 0.0
            br_rows.append([nm, motif_tag, s['n_total'], len(s['seqs']), '%.2f' % pct])
            print('  %-10s %-4s: %5d/%6d sequences (%.2f%%)' % (nm, motif_tag, len(s['seqs']), s['n_total'], pct))
    path = os.path.join(out_folder, 'per_br_occurrence.tsv')
    f = open(path, 'w', newline='')
    w = csv.writer(f, delimiter='\t')
    w.writerow(['replicate', 'motif', 'n_sequences', 'n_with_motif', 'pct_sequences_with_motif'])

    for row in br_rows:
        w.writerow(row)
    f.close()
    print('  wrote', path)
    print('\nStep 7: transcript context (this parses the 1.6 GB GENCODE GTF and takes several minutes) ...')
    exons, genes, utr5, utr3, cds = parse_gtf_lite(gtf_file)
    starts_by_chrom = {}

    for chrom in exons:
        starts_by_chrom[chrom] = build_index(exons[chrom])
    context_rows = []
    region_coords = {}

    for tag, fasta in (('WT', wt_fasta), ('R168C', r168c_fasta)):
        f = open(fasta)

        for line in f:
            if line.startswith('>'):
                rid = line[1:].strip()
                chrom, span = rid.split('|')[1].split(':')
                start, end = span.split('-')
                region_coords[tag, rid] = (chrom, int(start), int(end))
        f.close()

    for motif_tag in ('WT5', 'R5'):
        for coll_tag, source_tag in (('wt', 'WT'), ('r168c', 'R168C')):
            scan_tag = '%s_in_%s' % (motif_tag, coll_tag)
            fimo_path = os.path.join(fimo_folder, scan_tag, 'fimo.tsv')
            f = open(fimo_path)

            for line in f:
                if line.startswith('#') or line.startswith('motif_id'):
                    continue
                x = line.rstrip('\n').split('\t')

                if len(x) < 8:
                    continue
                seq_name = x[2]
                m_start = int(x[3]) - 1
                m_end = int(x[4])
                m_strand = x[7]
                key = (source_tag, seq_name)

                if key not in region_coords:
                    continue
                chrom, r_start, r_end = region_coords[key]
                g_start = r_start + m_start
                g_end = r_start + m_end
                label, sense = classify_occurrence(chrom, g_start, g_end, m_strand, exons, genes, utr5, utr3, cds)
                context_rows.append([motif_tag, coll_tag, seq_name, chrom, g_start, g_end, m_strand, label, sense])
    f.close()
    path = os.path.join(out_folder, 'transcript_context.tsv')
    f = open(path, 'w', newline='')
    w = csv.writer(f, delimiter='\t')
    w.writerow(['motif', 'condition_collection', 'sequence', 'chrom', 'genomic_start', 'genomic_end', 'motif_strand', 'location', 'sense'])

    for row in context_rows:
        w.writerow(row)
    f.close()
    print('  wrote', path, '(%d occurrences classified)' % len(context_rows))
    summary = {}

    for row in context_rows:
        key = (row[0], row[1], row[7])
        summary[key] = summary.get(key, 0) + 1
    path = os.path.join(out_folder, 'transcript_context_summary.tsv')
    f = open(path, 'w', newline='')
    w = csv.writer(f, delimiter='\t')
    w.writerow(['motif', 'collection', 'location', 'n_occurrences'])

    for key in sorted(summary):
        w.writerow([key[0], key[1], key[2], summary[key]])
    f.close()
    print('  wrote', path)
    print('\nMaking figures ...')
    wt_pwm = wt5['pwm']
    r_pwm = r5['pwm']
    plt.figure(figsize=(1.4 * max(wt5['w'], r5['w']) + 2, 3))
    render(plt.subplot(1, 2, 1), wt_pwm, 'WT5: %s' % wt5_consensus)
    render(plt.subplot(1, 2, 2), r_pwm, 'R5: %s' % r5_consensus)
    plt.suptitle('WT5 vs R5 PWMs')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_folder, 'wt5_r5_logos.png'), dpi=150)
    plt.close()
    labels = ['WT regions', 'R168C regions', 'shared regions']
    wt5_pct = [100.0 * len(scans['WT5_in_wt']['seqs']) / scans['WT5_in_wt']['n_total'], 100.0 * len(scans['WT5_in_r168c']['seqs']) / scans['WT5_in_r168c']['n_total'], 100.0 * len(scans['WT5_in_shared']['seqs']) / scans['WT5_in_shared']['n_total']]
    r5_pct = [100.0 * len(scans['R5_in_wt']['seqs']) / scans['R5_in_wt']['n_total'], 100.0 * len(scans['R5_in_r168c']['seqs']) / scans['R5_in_r168c']['n_total'], 100.0 * len(scans['R5_in_shared']['seqs']) / scans['R5_in_shared']['n_total']]
    plt.figure(figsize=(6, 4.5))
    x = [0, 1, 2]
    width = 0.35

    for i in x:
        _converted_5.append(i - width / 2)
    plt.bar(_converted_5, wt5_pct, width, label='WT5 (CCGTCTCAA)')

    for i in x:
        _converted_6.append(i + width / 2)
    plt.bar(_converted_6, r5_pct, width, label='R5 (GGTCTYGMWC)')
    plt.xticks(x, labels)
    plt.ylabel('% sequences with >= 1 motif match')
    plt.title('WT5 and R5 occurrence across sequence collections')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figures_folder, 'occurrence_bars.png'), dpi=150)
    plt.close()

    for row in br_rows:
        if row[1] == 'WT5':
            _converted_7.append(float(row[4]))
    br_wt5 = _converted_7

    for row in br_rows:
        if row[1] == 'R5':
            _converted_8.append(float(row[4]))
    br_r5 = _converted_8
    plt.figure(figsize=(7, 4.5))
    x = range(len(brs))
    width = 0.35

    for i in x:
        _converted_9.append(i - width / 2)
    plt.bar(_converted_9, br_wt5, width, label='WT5')

    for i in x:
        _converted_10.append(i + width / 2)
    plt.bar(_converted_10, br_r5, width, label='R5')
    plt.xticks(x, sorted(brs), fontsize=7)
    plt.ylabel('% sequences with >= 1 motif match')
    plt.title('WT5 and R5 occurrence per biological replicate')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figures_folder, 'per_br_occurrence.png'), dpi=150)
    plt.close()
    print('  wrote 4 figures to', figures_folder)
    print('\nFinished the WT5-R5 investigation')
    print('All outputs are in', out_folder)

import glob

import time

if __name__ == '__main__':
    main()
