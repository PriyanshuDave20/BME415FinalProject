# Final summary: WT5 vs R5 motif-to-target analysis

## 1. Research question

How might the difference between WT5 (CCGTCTCAA) and R5 (GGTCTYGMWC) motifs affect potential recognition of AR-V7, NRXN1, NRXN2 and NRXN3 transcripts?

## 2. Motif comparison

- WT5 consensus: `CCGTCTCAA` (width 9)
- R5 consensus: `GGTCTYGMWC` (width 10)
- PWM cosine similarity: ~0.76 (from earlier comparison)
- TOMTOM q-values: ~0.05 (borderline)
- The motifs demonstrate some sequence differences, but the evidence for definitively distinct binding specificities remains inconclusive.

## 3. Target transcript scanning

Scanning method: IUPAC-aware consensus matching (not calibrated PWM scanning).

| Gene | WT5 matches | R5 matches |
|------|------------|-----------|
| AR | 1 | 2 |
| NRXN1 | 2 | 2 |
| NRXN2 | 1 | 2 |
| NRXN3 | 1 | 2 |

Total WT5: 5, Total R5: 8

## 4. Differential classification

- Shared positions: 3
- WT5-only positions: 2
- R5-only positions: 5

R5 produced more consensus matches than WT5, including five R5-only positions. However, both motifs also shared some candidate recognition positions.

This suggests a possible difference in potential target-recognition patterns rather than a completely distinct recognition preference.

## 5. Feature annotation

| Feature | Count |
|---------|-------|
| exon | 13 |

The candidate matches were identified in exonic transcript regions. Some occur in regions potentially relevant to transcript processing, although their precise relationship to alternative splice events requires further verification.

## 6. Prioritization

| Confidence | N sites |
|------------|---------|
| high_computational | 2 |
| moderate_computational | 3 |
| exploratory | 2 |

### Top candidates

| Gene | Motif | Transcript | Score | Feature | Reasons |
|------|-------|------------|-------|---------|---------|
| NRXN3 | WT5 | ENST00000281127.11 | 3 | exon | near_splice_junction;shared_across_2_isoforms;potentially_accessible |
| NRXN3 | WT5/R5 | ENST00000335750.7 | 3 | exon | near_splice_junction;shared_across_2_isoforms;potentially_accessible |
| AR | R5 | ENST00000504326.5 | 3 | cryptic_exon_associated_region | near_splice_junction;potentially_accessible;AR-V7_specific_region |
| NRXN2 | WT5/R5 | ENST00000496291.2 | 2 | exon | shared_across_2_isoforms;potentially_accessible |
| NRXN1 | R5 | verified_transcript_id_required | 2 | exon | near_splice_junction;potentially_accessible |
| NRXN3 | WT5/R5 | ENST00000281127.11 | 1 | exon | shared_across_2_isoforms |
| NRXN2 | WT5 | ENST00000496291.2 | 1 | exon | shared_across_2_isoforms |

## 7. RNA accessibility

Important caveat: ViennaRNA was not installed during this analysis. Accessibility values were calculated using a GC-content-based heuristic, not actual RNA secondary structure predictions.

For each motif match, a 50 nt window on each side was extracted and its GC content was computed. A heuristic accessibility score of `1 - GC_fraction` was used.

### Results by motif

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

### Candidate site observations

- Several R5 candidate regions have heuristic accessibility scores >= 0.5.
- The strongest heuristic accessibility values occur in selected AR and NRXN1 candidate regions.
- NRXN3 contains multiple R5 candidate regions with moderate heuristic accessibility.
- Not all R5 sites exceed the 0.5 threshold, suggesting variation in sequence composition between candidate regions.
- WT5 also has potentially accessible regions, but its mean heuristic score is lower in this hypothetical comparison.

## 8. Integrated interpretation

The hypothetical analysis identifies a greater number of R5 consensus matches than WT5 matches across the four target transcript groups. Five R5-only positions and three shared positions suggest that the R168C-associated motif may have a somewhat different potential recognition pattern.

The prioritization analysis identifies two high-computational and three moderate-computational candidates. NRXN3 represents an interesting target because multiple transcript isoforms contain R5-associated candidate regions. AR-V7 is also potentially relevant because one hypothetical candidate is located near its cryptic exon-associated region.

The GC-content-based accessibility analysis indicates that several R5 candidate regions have relatively favourable sequence composition, with a mean heuristic accessibility score of 0.563 compared with 0.491 for WT5.

Together, these findings provide preliminary computational support for the hypothesis that the R168C-associated motif may recognize a somewhat different distribution of RNA regions across AR-V7 and neurexin transcripts.

However, the evidence does not establish a definitive change in binding specificity, actual KHDRBS2 binding, or an effect on alternative splicing.
