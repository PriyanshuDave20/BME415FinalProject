import os

from collections import defaultdict

OUT_DIR = 'results/motif_target_analysis/07_final_reports'

ANALYSIS_ROOT = 'results/motif_target_analysis'

def read_tsv(path):
    rows = []

    if not os.path.exists(path):
        return rows

    with open(path) as f:
        header = f.readline().strip().split('\t')

        for i in f:
            cols = i.rstrip('\n').split('\t')
            row = dict(zip(header, cols))
            rows.append(row)
    return rows

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    comp = read_tsv(os.path.join(ANALYSIS_ROOT, '01_motif_validation/motif_comparison.tsv'))
    all_matches = read_tsv(os.path.join(ANALYSIS_ROOT, '02_motif_scanning/all_motif_matches.tsv'))
    diff_sites = read_tsv(os.path.join(ANALYSIS_ROOT, '02_motif_scanning/differential_sites.tsv'))
    annotated = read_tsv(os.path.join(ANALYSIS_ROOT, '03_transcript_annotation/annotated_motif_sites.tsv'))
    candidates = read_tsv(os.path.join(ANALYSIS_ROOT, '05_functional_analysis/candidate_sites.tsv'))

    for m in all_matches:
        if m.get('Motif_name') == 'WT5':
            _converted_1.append(m)
    wt5_matches = _converted_1

    for m in all_matches:
        if m.get('Motif_name') == 'R5':
            _converted_2.append(m)
    r5_matches = _converted_2

    for r in diff_sites:
        if r.get('Classification') == 'shared':
            _converted_3.append(1)
    n_shared = sum(_converted_3)

    for r in diff_sites:
        if r.get('Classification') == 'WT5_only':
            _converted_4.append(1)
    n_wt_only = sum(_converted_4)

    for r in diff_sites:
        if r.get('Classification') == 'R5_only':
            _converted_5.append(1)
    n_r_only = sum(_converted_5)
    feature_counts = defaultdict(int)

    for row in annotated:
        feature_counts[row.get('Feature', 'unknown')] += 1
    confidence_counts = defaultdict(int)

    for row in candidates:
        confidence_counts[row.get('Confidence', 'unknown')] += 1
    gene_summary = defaultdict(lambda: {'WT5': 0, 'R5': 0})

    for row in all_matches:
        gene = row.get('Gene', '')
        motif = row.get('Motif_name', '')

        if gene and motif in ('WT5', 'R5'):
            gene_summary[gene][motif] += 1
    f = open(os.path.join(OUT_DIR, 'FINAL_SUMMARY.md'), 'w')
    f.write('# Final summary: WT5 vs R5 motif-to-target analysis\n\n')
    f.write('## 1. Research question\n\n')
    f.write('How might the difference between WT5 (CCGTCTCAA) and R5 (GGTCTYGMWC) motifs affect potential recognition of AR-V7, NRXN1, NRXN2 and NRXN3 transcripts?\n\n')
    f.write('## 2. Motif comparison\n\n')
    f.write('- WT5 consensus: `CCGTCTCAA` (width 9)\n')
    f.write('- R5 consensus:  `GGTCTYGMWC` (width 10)\n')
    f.write('- PWM cosine similarity: ~0.76 (from earlier comparison)\n')
    f.write('- TOMTOM q-values: ~0.05 (borderline)\n')
    f.write('- The evidence that WT5 and R5 are definitively different binding specificities is **not consistent** across comparisons.\n\n')
    f.write('## 3. Target transcript scanning\n\n')
    f.write('Scanning method: IUPAC-aware consensus matching (not calibrated PWM scanning).\n\n')
    f.write('| Gene | WT5 matches | R5 matches |\n')
    f.write('|------|------------|-----------|\n')

    for gene in ('AR', 'NRXN1', 'NRXN2', 'NRXN3'):
        f.write('| %s | %d | %d |\n' % (gene, gene_summary[gene]['WT5'], gene_summary[gene]['R5']))
    f.write('\nTotal WT5: %d, Total R5: %d\n\n' % (len(wt5_matches), len(r5_matches)))
    f.write('## 4. Differential classification\n\n')
    f.write('- Shared positions: %d\n' % n_shared)
    f.write('- WT5-only positions: %d\n' % n_wt_only)
    f.write('- R5-only positions: %d\n\n' % n_r_only)
    f.write('## 5. Feature annotation\n\n')
    f.write('| Feature | Count |\n|---------|-------|\n')

    for i in sorted(feature_counts):
        f.write('| %s | %d |\n' % (i, feature_counts[i]))
    f.write('\n')
    f.write('## 6. Prioritization\n\n')
    f.write('| Confidence | N sites |\n|------------|---------|\n')

    for i in ('high_computational', 'moderate_computational', 'exploratory'):
        f.write('| %s | %d |\n' % (i, confidence_counts.get(i, 0)))
    f.write('\n')
    f.write('## 7. Integrated interpretation\n\n')
    f.write('Based on the current evidence, the difference between WT5 and R5 motifs may or may not reflect a genuine change in binding specificity. The TOMTOM comparison was borderline and the motifs share substantial similarity. The target transcript scanning identified consensus matches in all four target genes, but without calibrated PWM thresholds or experimental validation it is not possible to conclude that the R168C mutation would change which transcript regions are recognized.\n\n')
    f.write('## 8. Limitations\n\n')
    f.write('See LIMITATIONS.md for a detailed discussion.\n')
    f.close()
    print('wrote FINAL_SUMMARY.md')
    f = open(os.path.join(OUT_DIR, 'METHODS.md'), 'w')
    f.write('# Methods\n\n')
    f.write('## Motif source\n\n')
    f.write('WT5 and R5 PWMs were extracted from STREME outputs in the Module 5 analysis. MEME-format files are in `01_motif_validation/`.\n\n')
    f.write('## Transcript scanning\n\n')
    f.write('IUPAC-aware consensus matching was used against the transcript sequences in `reference/target_transcripts/`.\n\n')
    f.write('## Transcript annotation\n\n')
    f.write('Motif matches in transcript coordinates were mapped to genomic coordinates using GENCODE v44 exon structures, accounting for strand.\n\n')
    f.write('## Accessibility\n\n')
    f.write('A GC-content heuristic was used as a rough proxy. ViennaRNA was not available.\n\n')
    f.write('## Prioritization\n\n')
    f.write('A 0-4 score was assigned based on feature, junction proximity, isoform sharing, and heuristic accessibility.\n')
    f.close()
    print('wrote METHODS.md')
    f = open(os.path.join(OUT_DIR, 'LIMITATIONS.md'), 'w')
    f.write('# Limitations\n\n')
    f.write('1. **Motif sample size**: Both WT5 and R5 were discovered from a limited number of sequences.\n')
    f.write('2. **Borderline similarity statistics**: TOMTOM q-values were around 0.05 in both directions.\n')
    f.write('3. **PWM uncertainty**: Neither motif was consistently recovered across all biological replicates.\n')
    f.write('4. **Consensus matching only**: Without FIMO or calibrated PWM thresholds, the scan used IUPAC-aware consensus matching, which is less sensitive than a proper PWM scan.\n')
    f.write('5. **No RNA structure prediction**: ViennaRNA was not available; only a GC-content heuristic was used.\n')
    f.write('6. **Transcript annotation**: Not all transcript isoforms may be captured in the current target_transcripts FASTA files.\n')
    f.write('7. **No experimental validation**: None of the candidate sites have been experimentally tested.\n')
    f.write('8. **Motif occurrence != binding**: Finding a motif sequence in a transcript does not prove that KHDRBS2 binds there.\n')
    f.close()
    print('wrote LIMITATIONS.md')
    results_path = os.path.join(OUT_DIR, 'RESULTS_TABLE.tsv')

    with open(results_path, 'w') as f:
        f.write('Gene\tWT5_matches\tR5_matches\n')

        for gene in ('AR', 'NRXN1', 'NRXN2', 'NRXN3'):
            f.write('%s\t%d\t%d\n' % (gene, gene_summary[gene]['WT5'], gene_summary[gene]['R5']))
        f.write('\nClassification\tCount\n')
        f.write('shared\t%d\n' % n_shared)
        f.write('WT5_only\t%d\n' % n_wt_only)
        f.write('R5_only\t%d\n' % n_r_only)
    print('wrote RESULTS_TABLE.tsv')
    manifest_path = os.path.join(OUT_DIR, 'reproducibility_manifest.txt')

    with open(manifest_path, 'w') as f:
        f.write('Reproducibility manifest\n')
        f.write('========================\n\n')
        f.write('Input files:\n')
        f.write('  reference/gencode.v44.annotation.gtf\n')
        f.write('  reference/target_transcripts/*.fa\n')
        f.write('  results/wt5_r5_investigation/motif_files/wt5.meme\n')
        f.write('  results/wt5_r5_investigation/motif_files/r5.meme\n\n')
        f.write('Scripts executed (in order):\n')
        f.write('  scripts/01_validate_motifs.py\n')
        f.write('  scripts/02_scan_target_transcripts.py\n')
        f.write('  scripts/03_map_transcript_features.py\n')
        f.write('  scripts/04_compare_isoforms.py\n')
        f.write('  scripts/05_analyze_rna_accessibility.py\n')
        f.write('  scripts/06_prioritize_candidate_sites.py\n')
        f.write('  scripts/07_generate_visualizations.py\n')
        f.write('  scripts/08_generate_final_report.py\n\n')
        f.write('Python version: 3.12\n')
        f.write('External dependencies used: none beyond standard library\n')
    print('wrote reproducibility_manifest.txt')

if __name__ == '__main__':
    main()
