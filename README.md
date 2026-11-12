KHDRBS2 easyCLIP-seq pipeline
BME415 final project
WT versus R168C comparison using dataset GSE154168 (12 FASTQ runs, human)

What this project does

This pipeline takes raw easyCLIP-seq FASTQ files for the RNA binding protein
KHDRBS2 (wild type and the R168C cancer mutant), processes and aligns them,
finds candidate binding regions, compares the two conditions, and discovers
enriched sequence motifs. It answers whether the R168C mutation changes
where KHDRBS2 binds or the sequence composition of its binding regions.


How to run

Run the scripts in this order with Python 3 (each is independently
rerunnable and skips work that is already done):

python load_prepare.py
python preprocess_reads.py
python align_reads.py
python find_binding_regions.py
python compare_conditions.py
python motifs.py

Requirements: Python 3 with matplotlib and cutadapt, plus a WSL Ubuntu
installation with minimap2, samtools, STAR, the MEME Suite (streme) and
umi_tools in a conda environment called bme415. scripts/setup_wsl_alignment.sh
installs those tools. The human reference genome (GRCh38, GENCODE v44) must
be downloaded into the reference folder (about 1 GB compressed, addresses in
the script comments). About 16 GB of free space is needed for intermediate
files.


Modules

load_prepare.py
- Validates all 12 FASTQ runs (read structure, mate pairing, barcodes)
- Assigns each run to WT or R168C with its biological replicate
- Output: results/sample_summary.tsv

preprocess_reads.py
- Trims the easyCLIP 3 prime adapters with cutadapt (the exact adapter
  sequences from the original workflow)
- Puts the UMI into the read name so duplicates can be removed later
- Discards pairs shorter than 14 nt and pairs with a malformed UMI
- Output: results/preprocessed/ and results/preprocessing_summary.tsv

align_reads.py
- Aligns the preprocessed reads to GRCh38 with minimap2 (STAR also
  supported) and sorts and indexes the BAM files
- Output: results/alignments_preprocessed/ and
  results/alignment_summary_preprocessed.tsv

find_binding_regions.py
- Counts fragments (first in pair, MAPQ 10 or higher) in 4 kb bins
- Calls enriched bins with a per-sample Poisson tail threshold and merges
  neighbouring bins into candidate regions
- Output: results/binding_regions/ (per-sample and combined region tables)

compare_conditions.py
- Verifies the replicate structure (4 biological replicates per condition,
  MiSeq runs are technical repeats of the HiSeq libraries)
- Removes PCR duplicates with umi_tools after merging technical runs
- Compares region overlap and normalized signal (FPM) between conditions
- Counts k-mers (k = 4, 5, 6) in the region sequences
- Runs a matched-replicate downsampling sensitivity analysis
- Compares the 4 kb method against a 250 bp crosslink-site method
- Output: results/comparison/, results/downsampling/,
  results/binding_regions/method_comparison.md

motifs.py
- Extracts 100 nt sequence windows around the condition-level regions
- Runs STREME separately on WT, R168C, each biological replicate and the
  downsampled regions
- Compares motifs between conditions by PWM similarity
- Output: results/motifs/ (FASTA files, STREME outputs, comparison tables,
  sequence logos and figures, motif_summary.md)


Main results

- WT and R168C share the large majority of candidate binding regions
- The apparent excess of R168C regions was mostly a sequencing depth and
  PCR duplication artifact, shown by deduplication and downsampling
- Both conditions yield the same GC-rich motif families with no evidence
  of an R168C specific motif
- Candidate regions are inferred from read coverage and are not confirmed
  binding sites; no statistical significance testing was applied


Folder layout

results             all pipeline outputs, figures and reports
src/
    load_prepare.py     module 1
    preprocess_reads.py module 2a
    align_reads.py      module 2b
    find_binding_regions.py  module 3
    compare_conditions.py    module 4
    motifs.py           module 5


Reports

results/alignment_comparison.md
  raw read versus preprocessed read alignment comparison
results/binding_regions/method_comparison.md
  4 kb bin method versus 250 bp crosslink-site method
results/comparison/comparison_summary.md
  full condition comparison report

