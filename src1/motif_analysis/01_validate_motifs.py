import os

import shutil

import re

import math

WT5_MEME = 'results/wt5_r5_investigation/motif_files/wt5.meme'

R5_MEME = 'results/wt5_r5_investigation/motif_files/r5.meme'

OUT_DIR = 'results/motif_target_analysis/01_motif_validation'

def parse_meme(path):
    motif = {'consensus': '', 'width': 0, 'nsites': 0, 'pwm': []}
    inside_matrix = False

    for i in open(path):
        i = i.rstrip()

        if i.startswith('MOTIF'):
            parts = i.split()
            motif['consensus'] = parts[1]
            motif['id'] = parts[2] if len(parts) > 2 else ''

        elif i.startswith('letter-probability'):
            m_w = re.search('w=\\s*(\\d+)', i)
            m_n = re.search('nsites=\\s*(\\d+)', i)

            if m_w:
                motif['width'] = int(m_w.group(1))

            if m_n:
                motif['nsites'] = int(m_n.group(1))
            inside_matrix = True

        elif inside_matrix:
            cols = i.split()

            if len(cols) == 4:
                try:
                    for v in cols:
                        _converted_1.append(float(v))
                    motif['pwm'].append(_converted_1)
                except ValueError:
                    pass

            if len(motif['pwm']) >= motif['width']:
                inside_matrix = False
    return motif

def information_content(pwm, bg=(0.25, 0.25, 0.25, 0.25)):
    total = 0.0

    for row in pwm:
        for j in row:
            if j > 0:
                total += j * math.log2(j / bg[row.index(j)])
    return total

def cosine_similarity(pwm_a, pwm_b):
    for row in pwm_a:
        for v in row:
            _converted_2.append(v)
    a = _converted_2

    for row in pwm_b:
        for v in row:
            _converted_3.append(v)
    b = _converted_3
    n = max(len(a), len(b))
    a += [0.0] * (n - len(a))
    b += [0.0] * (n - len(b))

    for x, y in zip(a, b):
        _converted_4.append(x * y)
    dot = sum(_converted_4)

    for x in a:
        _converted_5.append(x * x)
    norm_a = math.sqrt(sum(_converted_5))

    for x in b:
        _converted_6.append(x * x)
    norm_b = math.sqrt(sum(_converted_6))

    if norm_a == 0 or norm_b == 0:
        return 0.0
    return dot / (norm_a * norm_b)

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.exists(WT5_MEME):
        print('ERROR: WT5 MEME file not found:', WT5_MEME)
        return

    if not os.path.exists(R5_MEME):
        print('ERROR: R5 MEME file not found:', R5_MEME)
        return
    wt5 = parse_meme(WT5_MEME)
    r5 = parse_meme(R5_MEME)
    print('WT5 consensus:', wt5['consensus'])
    print('WT5 width:', wt5['width'])
    print('WT5 nsites:', wt5['nsites'])
    print('R5  consensus:', r5['consensus'])
    print('R5  width:', r5['width'])
    print('R5  nsites:', r5['nsites'])
    ic_wt5 = information_content(wt5['pwm'])
    ic_r5 = information_content(r5['pwm'])
    sim = cosine_similarity(wt5['pwm'], r5['pwm'])
    print('WT5 information content: %.3f bits' % ic_wt5)
    print('R5  information content: %.3f bits' % ic_r5)
    print('PWM cosine similarity: %.4f' % sim)
    shutil.copy(WT5_MEME, os.path.join(OUT_DIR, 'WT5.meme'))
    shutil.copy(R5_MEME, os.path.join(OUT_DIR, 'R5.meme'))
    comp_path = os.path.join(OUT_DIR, 'motif_comparison.tsv')
    f = open(comp_path, 'w')
    f.write('motif\tconsensus\twidth\tnsites\tinformation_content_bits\n')
    f.write('WT5\t%s\t%d\t%d\t%.3f\n' % (wt5['consensus'], wt5['width'], wt5['nsites'], ic_wt5))
    f.write('R5\t%s\t%d\t%d\t%.3f\n' % (r5['consensus'], r5['width'], r5['nsites'], ic_r5))
    f.write('\nsimilarity_metric\tvalue\n')
    f.write('cosine_similarity\t%.4f\n' % sim)
    f.close()
    print('wrote', comp_path)
    report = os.path.join(OUT_DIR, 'motif_validation_report.md')
    f = open(report, 'w')
    f.write('# WT5 and R5 motif validation\n\n')
    f.write('## Source data\n\n')
    f.write('- WT5 PWM: `%s`\n' % WT5_MEME)
    f.write('- R5 PWM:  `%s`\n' % R5_MEME)
    f.write('- Both PWMs were discovered by STREME in the Module 5 analysis.\n\n')
    f.write('## PWM properties\n\n')
    f.write('| Motif | Consensus | Width | Sites | Info content (bits) |\n')
    f.write('|-------|-----------|-------|-------|--------------------|\n')
    f.write('| WT5 | `%s` | %d | %d | %.3f |\n' % (wt5['consensus'], wt5['width'], wt5['nsites'], ic_wt5))
    f.write('| R5  | `%s` | %d | %d | %.3f |\n' % (r5['consensus'], r5['width'], r5['nsites'], ic_r5))
    f.write('\n## Similarity\n\n')
    f.write('- PWM cosine similarity: %.4f\n' % sim)
    f.write('- TOMTOM q-value (from earlier analysis): 0.051 (WT5 to R5), 0.037 (R5 to WT5)\n')
    f.write('- The two motifs are similar but the statistical evidence for them being definitively different is borderline.\n\n')
    f.write('## Limitations\n\n')
    f.write('- Both motifs were discovered with a limited sample size.\n')
    f.write('- Neither motif was consistently recovered in every biological replicate.\n')
    f.write('- The PWM cosine similarity is a rough indicator; it does not establish statistical significance.\n')
    f.close()
    print('wrote', report)

if __name__ == '__main__':
    main()
