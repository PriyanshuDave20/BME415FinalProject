import os

import sys

import urllib.request

import json

wt_sequence = 'MEEEKYLPELMAEKDSLDPSFVHASRLLAEEIEKFQGSDGKKEDEEKKYLDVISNKNIKLSERVLIPVKQYPKFNFVGKLLGPRGNSLKRLQEETGAKMSILGKGSMRDKAKEEELRKSGEAKYAHLSDELHVLIEVFAPPGEAYSRMSHALEEIKKFLVPDYNDEIRQEQLRELSYLNGSEDSGRGRGIRGRGIRIAPTAPSRGRGGAIPPPPPPGRGVLTPRGSTVTRGALPVPPVARGVPTPRARGAPTVPGYRAPPPPAHEAYEEYGYDDGYGGEYDDQTYETYDNSYATQTQSVPEYYDYGHGVSEDAYDSYAPEEWATTRSSLKAPPQRSARGGYREHPYGRY'

mutation_position = 168

wt_residue = 'R'

mutant_residue = 'C'

uniprot_accession = 'Q5VWX1'

out_folder = 'results/structural_analysis'

failures = []

def check(name, condition, detail=''):
    status = 'PASS' if condition else 'FAIL'
    print('  [%s] %s %s' % (status, name, detail))

    if not condition:
        failures.append(name + ' ' + detail)
    return condition

def fetch_uniprot(accession):
    url = 'https://rest.uniprot.org/uniprotkb/%s.json' % accession

    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            return json.loads(response.read().decode())
    except Exception as exc:
        print('  WARNING: could not reach UniProt (%s)' % exc)
        return None

def main():
    os.makedirs(out_folder, exist_ok=True)
    print('Step 1: sequence validation')
    print('  WT sequence length: %d' % len(wt_sequence))
    check('WT residue at position 168 is R', wt_sequence[mutation_position - 1] == wt_residue, '(found %s)' % wt_sequence[mutation_position - 1])
    mutant_sequence = list(wt_sequence)
    mutant_sequence[mutation_position - 1] = mutant_residue
    mutant_sequence = ''.join(mutant_sequence)
    check('mutant residue at position 168 is C', mutant_sequence[mutation_position - 1] == mutant_residue, '(found %s)' % mutant_sequence[mutation_position - 1])
    check('WT and mutant have identical lengths', len(wt_sequence) == len(mutant_sequence), '(%d vs %d)' % (len(wt_sequence), len(mutant_sequence)))
    differences = []

    for i in range(len(wt_sequence)):
        if wt_sequence[i] != mutant_sequence[i]:
            differences.append(i + 1)
    check('exactly one substitution', len(differences) == 1 and differences == [mutation_position], '(changed positions: %s)' % differences)
    uniprot_seq = None
    domains = []
    entry = fetch_uniprot(uniprot_accession)

    if entry is not None:
        uniprot_seq = entry['sequence']['value']
        print('  UniProt %s: %s (%s), length %d' % (uniprot_accession, entry['proteinDescription']['recommendedName']['fullName']['value'], entry['organism']['scientificName'], len(uniprot_seq)))
        check('UniProt sequence matches the WT sequence used here', uniprot_seq == wt_sequence)

        for i in entry.get('features', []):
            if i['type'] in ('Domain', 'Region', 'Motif'):
                domains.append((i['type'], i['location']['start']['value'], i['location']['end']['value'], i.get('description', '')))

        for ftype, start, end, description in domains:
            where = 'contains' if start <= mutation_position <= end else 'outside'
            print('  UniProt %s %d-%d: %s (residue 168 is %s this feature)' % (ftype, start, end, description or 'n/a', where))
    fasta_path = os.path.join(out_folder, 'r168c_sequences.fasta')

    with open(fasta_path, 'w') as f:
        f.write('>KHDRBS2_WT\n')

        for i in range(0, len(wt_sequence), 60):
            f.write(wt_sequence[i:i + 60] + '\n')
        f.write('>KHDRBS2_R168C\n')

        for i in range(0, len(mutant_sequence), 60):
            f.write(mutant_sequence[i:i + 60] + '\n')
    print('  wrote', fasta_path)
    report = os.path.join(out_folder, 'sequence_validation_report.md')

    with open(report, 'w') as f:
        f.write('# Sequence validation report\n\n')
        f.write('- WT sequence length: %d residues\n' % len(wt_sequence))
        f.write('- Residue 168 (1-based) in WT: %s\n' % wt_sequence[mutation_position - 1])
        f.write('- Residue 168 (1-based) in mutant: %s\n' % mutant_sequence[mutation_position - 1])
        f.write('- Number of substitutions: %d (R168C only)\n' % len(differences))
        f.write('- WT and mutant lengths identical: %s\n' % (len(wt_sequence) == len(mutant_sequence)))

        if uniprot_seq is not None:
            f.write('- UniProt %s sequence identical to the WT sequence used here: %s\n' % (uniprot_accession, uniprot_seq == wt_sequence))
        else:
            f.write('- UniProt cross-check: not performed (offline)\n')

        if domains:
            f.write('- UniProt annotated features (domain/region/motif):\n')

            for ftype, start, end, description in domains:
                where = 'contains' if start <= mutation_position <= end else 'outside'
                f.write('  - %s %d-%d %s (residue 168 is %s)\n' % (ftype, start, end, description or 'n/a', where))
        f.write('- FASTA written to: %s\n' % fasta_path)
        f.write('- Validation failures: %s\n' % (', '.join(failures) if failures else 'none'))
    print('  wrote', report)

    if failures:
        print('VALIDATION FAILED')
        sys.exit(1)
    print('Sequence validation passed.')

if __name__ == '__main__':
    main()
