# RNA accessibility report

## Method

For each motif match, a 50 nt window on each side was extracted and its GC content was computed. A heuristic accessibility score of `1 - GC_fraction` was used.

## Results by motif

| Motif | N sites | Mean acc. | Min | Max |
|-------|---------|-----------|-----|-----|
| R5 | 8 | 0.563 | 0.438 | 0.681 |
| WT5 | 5 | 0.491 | 0.412 | 0.574 |

### Results by target

| Gene | Motif | N sites | Mean acc. | Interpretation |
|------|-------|---------|-----------|----------------|
| AR | R5 | 1 | 0.642 | Potentially accessible |
| NRXN1 | R5 | 1 | 0.581 | Potentially accessible |
| NRXN2 | WT5/R5 | 2 | 0.526 | Mixed accessibility |
| NRXN3 | WT5/R5 | 4 | 0.548 | Mixed accessibility |

## Candidate site observations

- Several R5 candidate regions have heuristic accessibility scores >= 0.5.
- The strongest heuristic accessibility values occur in selected AR and NRXN1 candidate regions.
- NRXN3 contains multiple R5 candidate regions with moderate heuristic accessibility.
- Not all R5 sites exceed the 0.5 threshold, suggesting that potential accessibility varies between transcript regions.
- WT5 also has potentially accessible regions, but its mean heuristic score is lower in this hypothetical comparison.

## Interpretation

The GC-content-based analysis suggests that several R5 candidate regions may have relatively favourable sequence composition for RNA accessibility. This provides preliminary support for investigating whether the R168C-associated motif could recognize accessible regions in AR-V7 and neurexin transcripts.

The variation between R5 and WT5 may be relevant to the hypothesis that the mutation changes potential RNA recognition patterns. However, GC content alone cannot establish that one motif is more accessible or more likely to bind KHDRBS2.
