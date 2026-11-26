# KHDRBS2 WT vs R168C project

This project compares wild-type KHDRBS2 with the R168C mutation using easyCLIP-seq,
motif-to-transcript analysis, predicted protein structures, and protein-RNA docking.
The repo is organized as five connected stages:

1. Process easyCLIP-seq reads and identify enriched binding regions.
2. Compare the WT5 and R5 motifs against target transcripts.
3. Compare WT and R168C protein structures.
4. Prepare selected RNA target isoforms and 3D structures for docking.
5. Generate docking models and analyze WT/R168C protein-RNA contacts.

## Repository layout

```text
raw_data/                     Original easyCLIP-seq FASTQ inputs (ignored)
reference/                    GRCh38 genome, GENCODE annotation, target transcripts
src/                          Organized analysis and pipeline scripts
results/                      Final reports, tables, figures, selected structures, target 3D inputs
data_cache/                   Regenerable intermediates and archived/overlapping outputs (ignored)
WT/                           Top-10 WT protein-RNA docking models by target
R168C/                        Top-10 R168C protein-RNA docking models by target
3D Analysis/                  HDOCK .out records in Target_structs/
docking_contact_results/      Docking contact tables and mutation-impact summary
summaries/                    Convenience copies of cross-stage summaries
HDOCKlite/                    Local HDOCKlite tool (ignored)
.venv/                        Local Python environment (ignored)
```

### Source-code organization

```text
src/read_processing/          FASTQ validation, trimming/UMI handling, alignment, regions, comparison
src/motif_analysis/           STREME motifs, WT5/R5 target scanning, annotation, reports
src/structural_analysis/      Sequence/structure validation, RMSD and residue-168 analysis
src/target_isoform_docking/   Isoform selection, motif scanning, RNA-region/FASTA preparation
src/docking/                  Docking-model generation, contact analysis, impact summary
src/transcript_processing/    Target transcript extraction and curation helpers
```

## Data and cache policy

`results/` should contain curated outputs that support conclusions. Heavy intermediates
belong in `data_cache/`, grouped by topic. Current cache folders include:

- `data_cache/preprocessed/` and `data_cache/alignments*/`
- `data_cache/dedup/`
- `data_cache/comparison/`
- `data_cache/motif_target_analysis/`
- `data_cache/wt5_r5_investigation/`
- `data_cache/run_logs/`

`raw_data/`, `data_cache/`, `HDOCKlite/`, and `.venv/` are ignored. Some legacy
ignore rules remain for older `results/` paths in case an earlier script is rerun
before every script path has been modernized.

## Stage A: easyCLIP-seq read processing

Run from the repository root:

```powershell
py -3 src/read_processing/load_prepare.py
py -3 src/read_processing/preprocess_reads.py
py -3 src/read_processing/align_reads.py
py -3 src/read_processing/find_binding_regions.py
py -3 src/read_processing/compare_conditions.py
py -3 src/motif_analysis/motifs.py
```

Main curated outputs:

- `results/sample_summary.tsv`
- `results/preprocessing_summary.tsv`
- `results/binding_regions/`
- `results/comparison/comparison_summary.md`
- `results/downsampling/downsampling_summary.tsv`
- `results/motifs/motif_summary.md`

The workflow identifies WT5 and R5 as the motifs used by later stages. The motif
scan uses WT5 `CCGTCTCAA` and R5 `GGTCTYGMWC`.

## Stage B: WT5/R5 motif-to-target analysis

The numbered scripts are under `src/motif_analysis/`. The final report is the
source of truth for motif-to-target counts:

```text
results/motif_target_analysis/07_final_reports/FINAL_SUMMARY.md
```

The consolidated mutation-facing motif summary is:

```text
results/motif_target_analysis/07_final_reports/R5_WT5_motif_impact_summary.md
```

Key target-scan totals are:

| Gene | WT5 | R5 |
|---|---:|---:|
| AR | 1 | 2 |
| NRXN1 | 2 | 2 |
| NRXN2 | 1 | 2 |
| NRXN3 | 1 | 2 |
| **Total** | **5** | **8** |

The 13 positions are classified as 3 shared, 2 WT5-only, and 5 R5-only. All are
annotated as exonic transcript positions.

## Stage C: WT and R168C structural comparison

Run from the repository root:

```powershell
py -3 src/structural_analysis/01_validate_sequences.py
py -3 src/structural_analysis/02_inspect_structures.py
py -3 src/structural_analysis/03_compare_structures.py
py -3 src/structural_analysis/04_analyze_local_environment.py
py -3 src/structural_analysis/05_generate_report.py
```

Main outputs:

- `results/structural_analysis/WT/`
- `results/structural_analysis/R168C/`
- `results/structural_analysis/rmsd_comparison.tsv`
- `results/structural_analysis/local_environment_comparison.tsv`
- `results/structural_analysis/structural_analysis_report.md`

## Stage D: target isoform and RNA docking preparation

Run the target-preparation pipeline:

```powershell
py -3 src/target_isoform_docking/run_pipeline.py --stage all
```

Then run the standalone validation step:

```powershell
py -3 src/target_isoform_docking/09_final_validation.py
```

Main outputs:

- `results/target_isoform_docking/selected_isoforms/`
- `results/target_isoform_docking/motif_scanning/`
- `results/target_isoform_docking/docking_targets/`
- `results/target_isoform_docking/rna_sequences/`
- `results/target_isoform_docking/target_3d_structs/`
- `results/target_isoform_docking/reports/`

The target 3D inputs used for docking are `NRXN1.pdb`, `NRXN2.pdb`, `NRXN3.pdb`,
and `AR-V7.pdb`.

## Stage E: docking and contact analysis

Generate the top-10 rigid-body, HDOCKlite-compatible complex models:

```powershell
py -3 src/docking/run_docking_all_targets.py
```

The models are written to:

```text
WT/{AR_V7,NRXN1,NRXN2,NRXN3}/model_1.pdb ... model_10.pdb
R168C/{AR_V7,NRXN1,NRXN2,NRXN3}/model_1.pdb ... model_10.pdb
```

Run the contact analysis and write the integrated report:

```powershell
py -3 src/docking/run_contact_analysis_all_targets.py
py -3 src/docking/generate_mutation_impact_summary.py
```

Main outputs:

- `docking_contact_results/{target}/contact_summary.csv`
- `docking_contact_results/{target}/protein_residue_contacts.csv`
- `docking_contact_results/viable_contact_summary_all_targets.csv`
- `docking_contact_results/mutation_impact_summary_all_targets.md`

HDOCK `.out` records for the eight WT/R168C-by-target runs are retained in
`3D Analysis/Target_structs/`.

## Software requirements

- Python 3 with `matplotlib`, `numpy`, and `cutadapt`
- WSL Ubuntu with `minimap2`, `samtools`, `umi_tools`, and MEME Suite/STREME
- HDOCKlite for the retained HDOCK `.out` records

## Main results

The motif-to-target analysis reports 8 R5 consensus matches versus 5 WT5 matches
across AR, NRXN1, NRXN2, and NRXN3. Three positions are shared, five are R5-only,
and two are WT5-only. The heuristic accessibility mean is 0.563 for R5 versus
0.491 for WT5. See `FINAL_SUMMARY.md` for the authoritative counts.

Structurally, R168C does not introduce a clear local backbone rearrangement near
residue 168, but it replaces an exposed positive arginine with cysteine in an
acidic surface pocket.

At the 5 A docking-contact cutoff, R168C shows the largest contact increase for
NRXN1 (median viable residue pairs 13.0 for WT versus 52.5 for R168C), with a
smaller increase for NRXN3 and similar values for NRXN2 and AR-V7. Contact
results are computational predictions, not experimental binding measurements.
