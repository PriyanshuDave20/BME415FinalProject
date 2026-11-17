# 04_analyze_local_environment.py
# Analyzes the local structural environment of residue 168 in the WT and
# (when available) R168C structures: solvent accessibility, nearby charged,
# polar and hydrophobic residues, and potential interaction changes caused
# by the R -> C substitution.

import csv
import glob
import os
import statistics
import sys

from Bio.PDB import MMCIFParser, PDBParser, ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

out_folder = "results/structural_analysis"
position_of_interest = 168
local_start, local_end = 158, 178
contact_distance = 6.0      # angstrom, for 'nearby residue' contacts
hbond_distance = 3.5        # angstrom, generous donor-acceptor distance

charged_positive = ["ARG", "LYS", "HIS"]
charged_negative = ["ASP", "GLU"]
polar = ["ASN", "GLN", "SER", "THR", "TYR", "CYS", "TRP", "HIS"]
hydrophobic = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "PRO", "GLY"]


def load_structures():
    """Loads every structure file found (WT and R168C)."""
    structures = []
    parser = MMCIFParser(QUIET=True)
    pdb_parser = PDBParser(QUIET=True)
    for path in (glob.glob(os.path.join(out_folder, "*.cif"))
                 + glob.glob(os.path.join(out_folder, "*.pdb"))
                 + glob.glob(os.path.join(out_folder, "WT", "*.cif"))
                 + glob.glob(os.path.join(out_folder, "WT", "*.pdb"))
                 + glob.glob(os.path.join(out_folder, "R168C", "*.cif"))
                 + glob.glob(os.path.join(out_folder, "R168C", "*.pdb"))):
        name = os.path.basename(path).upper()
        if ("WT" not in name and "R168C" not in name):
            continue
        if path.endswith(".cif"):
            structure = parser.get_structure("s", path)
        else:
            structure = pdb_parser.get_structure("s", path)
        structures.append((name, structure))
    return structures


def analyze_environment(structure, name):
    """SASA (Shrake-Rupley, Bio.PDB) and contact analysis for the chain
    containing residue 168. Returns a result dictionary."""
    chain = None
    target = None
    for model in structure:
        for c in model:
            residues = [r for r in c if r.id[0] == " "]
            numbers = [r.id[1] for r in residues]
            if position_of_interest in numbers:
                chain = c
                target = [r for r in residues
                          if r.id[1] == position_of_interest][0]
                break
        if target is not None:
            break
    if target is None:
        return None

    # solvent accessibility: Shrake-Rupley on the whole chain, values per
    # residue are the fraction of the max sphere area kept accessible
    sasa = ShrakeRupley()
    sasa.compute(chain, level="R")

    residues = [r for r in chain if r.id[0] == " "]
    by_number = {r.id[1]: r for r in residues}

    # nearby residues within the contact distance of any atom of residue 168
    nearby = []
    for r in residues:
        if r.id[1] == position_of_interest:
            continue
        min_dist = None
        for atom_a in target:
            for atom_b in r:
                dist = atom_a - atom_b
                if min_dist is None or dist < min_dist:
                    min_dist = dist
        if min_dist is not None and min_dist <= contact_distance:
            nearby.append((r.id[1], r.resname, min_dist))

    result = {
        "name": name,
        "resname_168": target.resname,
        "sasa_168": round(target.sasa, 1),
        "sasa_local": {},
        "nearby_positive": [(n, rn, round(d, 1)) for (n, rn, d) in nearby
                            if rn in charged_positive],
        "nearby_negative": [(n, rn, round(d, 1)) for (n, rn, d) in nearby
                            if rn in charged_negative],
        "nearby_polar": [(n, rn, round(d, 1)) for (n, rn, d) in nearby
                         if rn in polar],
        "nearby_hydrophobic": [(n, rn, round(d, 1)) for (n, rn, d) in nearby
                               if rn in hydrophobic],
        "n_contacts": len(nearby),
    }
    for n in range(local_start, local_end + 1):
        if n in by_number and "sasa" in by_number[n].xtra:
            result["sasa_local"][n] = round(by_number[n].sasa, 1)
    return result


def main():
    structures = load_structures()
    wt = [(n, s) for (n, s) in structures if "WT" in n.upper()]
    mut = [(n, s) for (n, s) in structures if "R168C" in n.upper()]

    if not wt:
        print("ERROR: no WT structure found in", out_folder)
        sys.exit(1)

    results = []
    for name, structure in wt + mut:
        print("Analyzing", name, "...")
        result = analyze_environment(structure, name)
        if result is None:
            print("  ERROR: residue %d not found in this structure"
                  % position_of_interest)
            sys.exit(1)
        results.append(result)
        print("  residue %d is %s, SASA %.1f A^2, %d residues within %.1f A"
              % (position_of_interest, result["resname_168"],
                 result["sasa_168"], result["n_contacts"], contact_distance))
        for label in ("nearby_positive", "nearby_negative", "nearby_polar",
                      "nearby_hydrophobic"):
            print("    %s: %s" % (label, result[label]))

    # comparison when both variants are available
    comparison_pending = len(mut) == 0
    if comparison_pending:
        print("\nWARNING: no R168C structure available yet - the WT-vs-mutant")
        print("environment comparison will be done after ColabFold predictions")
        print("are placed in", out_folder)
        print("Submit results/structural_analysis/r168c_sequences.fasta to")
        print("ColabFold (see colabfold_instructions.md in the same folder)")

    # write the environment table
    path = os.path.join(out_folder, "local_environment_comparison.tsv")
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["structure", "resname_168", "sasa_168",
                         "n_residues_within_6A", "nearby_positive",
                         "nearby_negative", "nearby_polar",
                         "nearby_hydrophobic"])
        for r in results:
            writer.writerow([r["name"], r["resname_168"], r["sasa_168"],
                             r["n_contacts"],
                             "; ".join("%d %s %.1fA" % c
                                       for c in r["nearby_positive"]),
                             "; ".join("%d %s %.1fA" % c
                                       for c in r["nearby_negative"]),
                             "; ".join("%d %s %.1fA" % c
                                       for c in r["nearby_polar"]),
                             "; ".join("%d %s %.1fA" % c
                                       for c in r["nearby_hydrophobic"])])
    print("\n  wrote", path)

    if comparison_pending:
        print("\nSTATUS: WT analysis complete, R168C comparison pending")
        print("ColabFold predictions")


if __name__ == "__main__":
    main()
