import subprocess

import sys

import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

STAGES = [('validate', '01_validate_references.py'), ('extract', '02_extract_transcripts.py'), ('select', '03_select_isoforms.py'), ('scan', '04_scan_motifs.py'), ('map', '05_map_coordinates.py'), ('extract_regions', '06_extract_docking_regions.py'), ('rna', '07_prepare_rna_fasta.py'), ('report', '08_generate_reports.py')]

def run_stage(stage, script):
    path = os.path.join(SCRIPT_DIR, script)

    if not os.path.exists(path):
        print('ERROR: script not found:', path)
        return False
    print('=' * 60)
    print('Stage: %s (%s)' % (stage, script))
    print('=' * 60)
    result = subprocess.run([sys.executable, path])

    if result.returncode != 0:
        print('ERROR: stage %s failed' % stage)
        return False
    print()
    return True

def main():
    args = sys.argv[1:]

    if not args or args[0] != '--stage':
        print('Usage: run_pipeline.py --stage <stage>')

        for s in STAGES:
            _converted_1.append(s[0])
        print('Stages:', ', '.join(_converted_1), 'or all')
        sys.exit(1)
    stage = args[1] if len(args) > 1 else ''

    if stage == 'all':
        for s, script in STAGES:
            ok = run_stage(s, script)

            if not ok:
                sys.exit(1)
        print('All stages completed.')
    else:
        for s, script in STAGES:
            if s == stage:
                ok = run_stage(s, script)

                if not ok:
                    sys.exit(1)
                return
        print('Unknown stage:', stage)
        sys.exit(1)

if __name__ == '__main__':
    main()
