
# Functional interpretation of candidate sites

## Prioritization method

Each site was scored on a 0-4 scale:

1. +1 if the site is in a CDS or 3' UTR
2. +1 if the site is near a splice junction (within 50 bp)
3. +1 if the site is shared across multiple transcript isoforms
4. +1 if the heuristic accessibility score >= 0.5

A score of 3-4 is high computational confidence; 2 is moderate; 0-1 is exploratory.

## Results

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

## WT5 versus R5 comparison

R5 produced several candidate sites in transcript regions potentially associated with alternative RNA processing. Some R5 candidates were located near splice junctions and were shared across multiple transcript isoforms.

WT5 produced fewer comparable high-scoring candidate matches in these particular regions. This suggests a possible difference in the potential target-recognition patterns of the WT and R168C-associated motifs.

However, the difference requires confirmation using matched PWM scoring thresholds and sequence-background controls.

## Interpretation

The results suggest that the R168C-associated R5 motif may have a different potential recognition pattern across selected AR and neurexin transcripts.

NRXN3 represents the strongest candidate because multiple transcript isoforms contain potentially accessible R5 sites near alternative splicing-associated regions. The AR-V7 candidate is also interesting because of its potential proximity to the cryptic exon-associated splice event.

These findings provide preliminary computational support for investigating whether altered KHDRBS2 RNA recognition could influence transcript processing.

The results do not establish actual KHDRBS2 binding or demonstrate changes in alternative splicing.
