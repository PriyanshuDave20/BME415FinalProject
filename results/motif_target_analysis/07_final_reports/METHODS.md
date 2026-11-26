# Methods

## Motif source

WT5 and R5 PWMs were extracted from STREME outputs in the Module 5 analysis. MEME-format files are in `01_motif_validation/`.

## Transcript scanning

IUPAC-aware consensus matching was used against the transcript sequences in `reference/target_transcripts/`.

## Transcript annotation

Motif matches in transcript coordinates were mapped to genomic coordinates using GENCODE v44 exon structures, accounting for strand.

## Accessibility

A GC-content heuristic was used as a rough proxy. ViennaRNA was not available.

## Prioritization

A 0-4 score was assigned based on feature, junction proximity, isoform sharing, and heuristic accessibility.
