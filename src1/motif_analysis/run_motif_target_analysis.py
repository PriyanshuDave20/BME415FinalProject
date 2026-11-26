import subprocess

import sys

import os

SCRIPTS_DIR = 'scripts'

STAGES = [('validate', '01_validate_motifs.py'), ('scan', '02_scan_target_transcripts.py'), ('annotate', '03_map_transcript_features.py'), ('isoforms', '04_compare_isoforms.py'), ('accessibility', '05_analyze_rna_accessibility.py'), ('prioritize', '06_prioritize_candidate_sites.py'), ('visualize', '07_generate_visualizations.py'), ('report', '08_generate_final_report.py')]

def run_stage(stage_name, script_name):
    script_path = os.path.join(SCRIPTS_DIR, script_name)

    if not os.path.exists(script_path):
        print('ERROR: script not found:', script_path)
        return False
    print('=' * 60)
    print('Running stage: %s (%s)' % (stage_name, script_name))
    print('=' * 60)
    result = subprocess.run([sys.executable, script_path])

    if result.returncode != 0:
        print('ERROR: stage %s failed with exit code %d' % (stage_name, result.returncode))
        return False
    print()
    return True

def main():
    args = sys.argv[1:]

    if not args or args[0] not in ('--stage',):
        print('Usage: python run_motif_target_analysis.py --stage <stage>')

        for s in STAGES:
            _converted_1.append(s[0])
        print('Stages:', ', '.join(_converted_1), 'or all')
        sys.exit(1)
    stage = args[1] if len(args) > 1 else ''

    if stage == 'all':
        for stage_name, script_name in STAGES:
            ok = run_stage(stage_name, script_name)

            if not ok:
                sys.exit(1)
        print('All stages completed.')
    else:
        found = False

        for stage_name, script_name in STAGES:
            if stage_name == stage:
                ok = run_stage(stage_name, script_name)

                if not ok:
                    sys.exit(1)
                found = True
                break

        if not found:
            print('Unknown stage:', stage)

            for s in STAGES:
                _converted_2.append(s[0])
            print('Available stages:', ', '.join(_converted_2), 'or all')
            sys.exit(1)

if __name__ == '__main__':
    main()
