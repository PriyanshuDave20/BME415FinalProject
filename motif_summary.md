# Module 5: motif discovery summary (STREME)

STREME (MEME Suite 5.5.9) was run on 100 nt windows centred on
the condition-level candidate regions from the Module 4
biological-replicate analysis (deduplicated BAMs, deeper
HiSeq run per library). Without control sequences STREME
shuffles the primary sequences preserving dinucleotide
frequencies, so the background comes from each condition's
own regions. Settings: --minw 4 --maxw 12 --nmotifs 5
(condition and downsampled runs) / 1 (replicate runs),
--seed 154168.

## Sequences analysed

- DS-R168C: 12467 sequences kept (1 excluded with N, 0 duplicates)
- DS-WT: 12335 sequences kept (2 excluded with N, 0 duplicates)
- R168C: 14634 sequences kept (2 excluded with N, 0 duplicates)
- R168C-BR1: 10402 sequences kept (2 excluded with N, 0 duplicates)
- R168C-BR2: 13907 sequences kept (1 excluded with N, 0 duplicates)
- R168C-BR3: 4010 sequences kept (1 excluded with N, 0 duplicates)
- R168C-BR4: 2440 sequences kept (1 excluded with N, 0 duplicates)
- WT: 12320 sequences kept (2 excluded with N, 0 duplicates)
- WT-BR1: 8037 sequences kept (2 excluded with N, 0 duplicates)
- WT-BR2: 11679 sequences kept (2 excluded with N, 0 duplicates)
- WT-BR3: 1180 sequences kept (1 excluded with N, 0 duplicates)
- WT-BR4: 2205 sequences kept (1 excluded with N, 0 duplicates)

Windows are region-centred (Module 3 stored no per-base
summits), so motif signal is diluted by the 4 kb bin width.
Sequences are genomic DNA; candidate regions are not
reliably stranded, so no per-strand orientation was used.

## Motifs found

### WT (5 motifs)

- GCTGGGATTACA (width 12, 842 sites, 845 sequences, E = 9.50e-23)
- CTCRGCCTCCCR (width 12, 972 sites, 985 sequences, E = 6.80e-22)
- RGYGYRGTGGC (width 11, 889 sites, 906 sequences, E = 1.80e-18)
- CAGCCTGGSCAA (width 12, 664 sites, 677 sequences, E = 6.40e-14)
- CCGTCTCAA (width 9, 634 sites, 685 sequences, E = 1.50e-10)

### R168C (5 motifs)

- GCCACYRCRCYC (width 12, 1213 sites, 1267 sequences, E = 7.00e-28)
- GGGAGGCYGAGG (width 12, 976 sites, 988 sequences, E = 1.10e-26)
- GCTGGGATTACA (width 12, 962 sites, 982 sequences, E = 4.70e-26)
- GCCTGGSCAACA (width 12, 849 sites, 863 sequences, E = 3.70e-23)
- GGTCTYGMWC (width 10, 720 sites, 778 sequences, E = 4.30e-11)

## WT vs R168C comparison

Motifs were matched by best-offset PWM cosine similarity in
both orientations; similarity >= 0.85 counts as the same
motif (see motif_comparison.tsv).

- shared motifs: 4
  - CAGCCTGGSCAA ~ GCCTGGSCAACA (similarity 1.00)
  - GCTGGGATTACA ~ GCTGGGATTACA (similarity 1.00)
  - RGYGYRGTGGC ~ GCCACYRCRCYC (similarity 1.00)
  - CTCRGCCTCCCR ~ GGGAGGCYGAGG (similarity 1.00)
- CCGTCTCAA is listed as WT-associated, but it is a close
  relative of the shared motif family
- GGTCTYGMWC is listed as R168C-associated, but it is a close
  relative of the shared motif family

## Motif recurrence across biological replicates

STREME was run separately on each of the 4 biological
replicates per condition (1 motif each). A condition motif
counts as detected in a replicate when the PWM similarity is
>= 0.85 (see replicate_motif_summary.tsv):

- WT GCTGGGATTACA: detected in 4 of 4 replicates (1.00, 1.00, 0.99, 1.00)
- WT CTCRGCCTCCCR: detected in 0 of 4 replicates (0.56, 0.56, 0.63, 0.51)
- WT RGYGYRGTGGC: detected in 0 of 4 replicates (0.60, 0.61, 0.53, 0.54)
- WT CAGCCTGGSCAA: detected in 0 of 4 replicates (0.57, 0.58, 0.59, 0.41)
- WT CCGTCTCAA: detected in 0 of 4 replicates (0.51, 0.52, 0.58, 0.57)
- R168C GCCACYRCRCYC: detected in 0 of 4 replicates (0.50, 0.56, 0.51, 0.57)
- R168C GGGAGGCYGAGG: detected in 0 of 4 replicates (0.50, 0.53, 0.50, 0.52)
- R168C GCTGGGATTACA: detected in 4 of 4 replicates (1.00, 1.00, 1.00, 1.00)
- R168C GCCTGGSCAACA: detected in 0 of 4 replicates (0.64, 0.57, 0.65, 0.57)
- R168C GGTCTYGMWC: detected in 0 of 4 replicates (0.57, 0.51, 0.57, 0.51)

## Downsampled sensitivity check

STREME was re-run on the Module 4 downsampled region windows
(equalized sequencing depth within matched replicate pairs):

- WT GCTGGGATTACA: similarity 1.00 (recovered)
- WT CTCRGCCTCCCR: similarity 1.00 (recovered)
- WT RGYGYRGTGGC: similarity 1.00 (recovered)
- WT CAGCCTGGSCAA: similarity 1.00 (recovered)
- WT CCGTCTCAA: similarity 0.60 (NOT recovered)
- R168C GCCACYRCRCYC: similarity 1.00 (recovered)
- R168C GGGAGGCYGAGG: similarity 0.99 (recovered)
- R168C GCTGGGATTACA: similarity 1.00 (recovered)
- R168C GCCTGGSCAACA: similarity 1.00 (recovered)
- R168C GGTCTYGMWC: similarity 0.57 (NOT recovered)

## Interpretation limits

- An enriched motif in candidate regions is not a confirmed
  KHDRBS2 binding motif.
- A motif seen in only one condition can reflect candidate-
  region selection and sequencing depth, not binding
  specificity.
- Module 4 showed the two conditions' region sets overlap
  strongly and their overall k-mer composition is nearly
  identical, so large motif differences would be surprising.
- STREME finds enriched patterns; it does not establish
  biological function or a direct biochemical effect of the
  R168C mutation.

wt_motif_id	r168c_motif_id	wt_consensus	r168c_consensus	width_wt	width_r168c	similarity	wt_evalue	r168c_evalue	wt_seqs_with_motif	r168c_seqs_with_motif	relationship
# STREME-4	STREME-4	CAGCCTGGSCAA	GCCTGGSCAACA	12	12	1.000	6.40e-14	3.70e-23	677	863	shared
# STREME-1	STREME-3	GCTGGGATTACA	GCTGGGATTACA	12	12	1.000	9.50e-23	4.70e-26	845	982	shared
# STREME-3	STREME-1	RGYGYRGTGGC	GCCACYRCRCYC	11	12	0.998	1.80e-18	7.00e-28	906	1267	shared
# STREME-2	STREME-2	CTCRGCCTCCCR	GGGAGGCYGAGG	12	12	0.996	6.80e-22	1.10e-26	985	988	shared
# STREME-5		CCGTCTCAA		9		1.50e-10		685		WT-associated (best sim to R168C: 0.76)
# STREME-5		GGTCTYGMWC		10		4.30e-11		778	R168C-associated (best sim to WT: 0.76)
