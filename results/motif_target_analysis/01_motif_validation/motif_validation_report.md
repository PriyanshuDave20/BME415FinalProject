# WT5 and R5 motif validation

## Source data

- WT5 PWM: `results/wt5_r5_investigation/motif_files/wt5.meme`
- R5 PWM:  `results/wt5_r5_investigation/motif_files/r5.meme`
- Both PWMs were discovered by STREME in the Module 5 analysis.

## PWM properties

| Motif | Consensus | Width | Sites | Info content (bits) |
|-------|-----------|-------|-------|--------------------|
| WT5 | `CCGTCTCAA` | 9 | 634 | 14.660 |
| R5  | `GGTCTYGMWC` | 10 | 720 | 14.555 |

## Similarity

- PWM cosine similarity: 0.1762
- TOMTOM q-value (from earlier analysis): 0.051 (WT5 to R5), 0.037 (R5 to WT5)
- The two motifs are similar but the statistical evidence for them being definitively different is borderline.

## Limitations

- Both motifs were discovered with a limited sample size.
- Neither motif was consistently recovered in every biological replicate.
- The PWM cosine similarity is a rough indicator; it does not establish statistical significance.
