# 02_inspect_structures.py
# Inspects the AlphaFold structure files in the project: chain identifiers,
# residue numbering, missing residues, and pLDDT confidence values.

import csv
import glob
import os
import statistics
import sys

from Bio.PDB import MMCIFParser, PDBParser, ShrakeRupley

out_folder = "results/structural_analysis"
local_start = 158           # local region of interest
local_end = 178
position_of_interest = 168

# pLDDT values reported by the user (for verification against the model)
reported_plddt = {156: 95.12, 157: 94.19, 158: 94.81, 159: 94.69, 160: 92.44,
                  161: 89.50, 162: 86.50, 163: 76.12, 164: 83.69,
                  165: 88.25, 166: 91.44, 167: 89.81, 168: 90.25,
                  169: 92.12, 170: 92.94}


def find_structure_files():
    """Finds AlphaFold model files (cif or pdb) in the project."""
    found = []
    for pattern in ("*.cif", "*.pdb"):
        for path in glob.glob(pattern) + glob.glob(os.path.join("results",
                                            "structural_analysis", pattern)):
            found.append(path)
    return found


def inspect_structure(path):
    """Reads one structure file and returns its metadata:
    chains, residue numbering, pLDDT values, missing residues."""
    if path.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", path)
    model = structure[0]

    chains = {}
    for chain in model:
        residues = [r for r in chain if r.id[0] == " "]
        numbers = sorted(r.id[1] for r in residues)
        plddt = {}
        for r in residues:
            if "CA" in r:
                plddt[r.id[1]] = r["CA"].bfactor
        # missing residues = gaps in the numbering
        missing = []
        for a, b in zip(numbers, numbers[1:]):
            if b - a > 1:
                missing.append("%d-%d" % (a, b - 1))
        chains[chain.id] = {"n_residues": len(residues),
                            "first": numbers[0], "last": numbers[-1],
                            "missing": missing, "plddt": plddt,
                            "residues": residues}
    return {"path": path, "chains": chains}


# UniProt annotations for Q5VWX1 (verified 01_validate_sequences.py)
UNIPROT_FEATURES = [
    ("KH domain", 65, 135),
    ("Disordered region", 182, 284),
    ("Disordered region", 319, 349),
]

# focus sub-regions requested for the structural analysis
FOCUS_REGIONS = [
    ("KH domain (65-135)", 65, 135),
    ("KH C-terminal part (120-145)", 120, 145),
    ("Linker before disordered region (136-181)", 136, 181),
    ("Linker around R168 (150-181)", 150, 181),
    ("Disordered region start (182-200)", 182, 200),
    ("Annotated disordered region (182-284)", 182, 284),
    ("Annotated disordered region (319-349)", 319, 349),
]


def region_context(chain):
    """Region-level analysis for one chain: pLDDT statistics, per-residue
    SASA (Shrake-Rupley) and the fraction of exposed residues per region,
    plus the contacts between the linker and the KH domain."""
    residues = [r for r in chain if r.id[0] == " "]
    by_number = {r.id[1]: r for r in residues}

    ShrakeRupley().compute(chain, level="R")

    region_rows = []
    for name, start, end in FOCUS_REGIONS + UNIPROT_FEATURES:
        vals = [by_number[n]["CA"].bfactor for n in range(start, end + 1)
                if n in by_number]
        sasa_vals = [by_number[n].sasa for n in range(start, end + 1)
                     if n in by_number]
        if not vals:
            continue
        # a residue is counted as exposed when its SASA is above 30 A^2
        exposed = sum(1 for s in sasa_vals if s > 30)
        region_rows.append([name, start, end, len(vals),
                            round(statistics.mean(vals), 2),
                            round(min(vals), 2),
                            round(statistics.mean(sasa_vals), 1),
                            round(100.0 * exposed / len(sasa_vals), 1)])
    return region_rows, by_number, residues


def linker_kh_contacts(chain):
    """Counts linker residues (136-181, containing R168) whose atoms come
    within 6 A of the KH domain (65-135), and lists the closest KH-domain
    contacts of residue 168."""
    residues = [r for r in chain if r.id[0] == " "]
    by_number = {r.id[1]: r for r in residues}
    linker = [r for r in residues if 136 <= r.id[1] <= 181]
    kh = [r for r in residues if 65 <= r.id[1] <= 135]

    contacting = []
    for r in linker:
        touching = set()
        for atom_a in r:
            for k in kh:
                for atom_b in k:
                    if atom_a - atom_b <= 6.0:
                        touching.add(k.id[1])
        if touching:
            contacting.append((r.id[1], r.resname,
                               min(touching), max(touching), len(touching)))

    # closest KH-domain residue to residue 168
    r168 = by_number[position_of_interest]
    closest = None
    for k in kh:
        dist = min(atom_a - atom_b for atom_a in r168 for atom_b in k)
        if closest is None or dist < closest[2]:
            closest = (k.id[1], k.resname, dist)

    # pLDDT of residue 168 and its neighbours relative to the two regions
    return {"linker_n": len(linker),
            "linker_contacting_kh": contacting,
            "closest_kh_to_168": closest,
            "plddt_168": by_number[position_of_interest]["CA"].bfactor}


def main():
    os.makedirs(out_folder, exist_ok=True)
    files = find_structure_files()
    if not files:
        print("ERROR: no structure files (.cif or .pdb) found in the project")
        sys.exit(1)

    print("Found %d structure file(s)" % len(files))
    summary_rows = []
    context_rows = []
    for path in files:
        variant = "R168C" if "R168C" in os.path.basename(path).upper() \
            else "WT"
        print("\nStructure:", path, "(treated as %s)" % variant)
        meta = inspect_structure(path)
        for chain_id, c in meta["chains"].items():
            print("  chain %s: %d residues, numbering %d-%d"
                  % (chain_id, c["n_residues"], c["first"], c["last"]))
            if c["missing"]:
                print("  missing residues: %s" % ", ".join(c["missing"]))
            else:
                print("  missing residues: none")

            plddt = c["plddt"]
            local_vals = [plddt[n] for n in range(local_start, local_end + 1)
                          if n in plddt]
            print("  mean pLDDT (all residues): %.2f"
                  % statistics.mean(plddt.values()))
            print("  mean pLDDT (residues %d-%d): %.2f"
                  % (local_start, local_end, statistics.mean(local_vals)))
            print("  pLDDT at residue %d: %.2f"
                  % (position_of_interest,
                     plddt.get(position_of_interest, float("nan"))))

            # verify the reported pLDDT values against the model
            mismatches = []
            for n in sorted(reported_plddt):
                model_value = round(plddt.get(n, -1), 2)
                if abs(model_value - reported_plddt[n]) > 0.01:
                    mismatches.append("%d: model %.2f vs reported %.2f"
                                      % (n, model_value, reported_plddt[n]))
            if mismatches:
                print("  WARNING: reported pLDDT values do not match: %s"
                      % "; ".join(mismatches))
            else:
                print("  reported pLDDT values match the model exactly")

            summary_rows.append([variant, os.path.basename(path), chain_id,
                                 c["n_residues"], c["first"], c["last"],
                                 "%.2f" % statistics.mean(plddt.values()),
                                 "%.2f" % statistics.mean(local_vals),
                                 "%.2f" % plddt.get(position_of_interest,
                                                    float("nan")),
                                 "; ".join(c["missing"]) or "none"])

            # region-level context (UniProt annotations + focus regions)
            # re-parse the file to get the Bio.PDB chain object
            if path.endswith(".cif"):
                bio_parser = MMCIFParser(QUIET=True)
            else:
                bio_parser = PDBParser(QUIET=True)
            bio_chain = bio_parser.get_structure("s", path)[0][chain_id]
            region_rows, by_number, residues = region_context(bio_chain)
            for row in region_rows:
                context_rows.append([variant, os.path.basename(path)] + row)

            # where does residue 168 sit structurally?
            r168 = by_number[position_of_interest]
            linker_rows = [row for row in region_rows
                           if row[0].startswith("Linker before")]
            if linker_rows:
                linker_mean_plddt = linker_rows[0][4]
                print("  linker 136-181: mean pLDDT %.2f (R168 pLDDT %.2f)"
                      % (linker_mean_plddt, plddt[position_of_interest]))

            contacts = linker_kh_contacts(bio_chain)
            print("  R168 is %s at position %d, %d residues before the"
                  " annotated disordered region (182)"
                  % (r168.resname, position_of_interest,
                     182 - position_of_interest))
            print("  linker residues in contact with the KH domain: %d of %d"
                  % (len(contacts["linker_contacting_kh"]),
                     contacts["linker_n"]))
            if contacts["linker_contacting_kh"]:
                spans = ["%d (%s, KH %d-%d)" % (n, rn, kmin, kmax)
                         for (n, rn, kmin, kmax, _c)
                         in contacts["linker_contacting_kh"]]
                print("    contacting linker residues: %s" % ", ".join(spans))
            ck = contacts["closest_kh_to_168"]
            print("  closest KH-domain residue to R168: %d (%s) at %.1f A"
                  % (ck[0], ck[1], ck[2]))
            # the region-level table includes the linker and KH rows

    # write the metadata summary table
    path = os.path.join(out_folder, "structural_metadata_summary.tsv")
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["variant", "file", "chain", "n_residues",
                         "first_residue", "last_residue", "mean_plddt_all",
                         "mean_plddt_158_178", "plddt_168",
                         "missing_residues"])
        for row in summary_rows:
            writer.writerow(row)
    print("\n  wrote", path)

    # write the region context table
    path = os.path.join(out_folder, "region_context.tsv")
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["variant", "file", "region", "start", "end",
                         "n_residues", "mean_plddt", "min_plddt",
                         "mean_sasa_A2", "pct_exposed_gt_30A2"])
        for row in context_rows:
            writer.writerow(row)
    print("  wrote", path)


if __name__ == "__main__":
    main()
