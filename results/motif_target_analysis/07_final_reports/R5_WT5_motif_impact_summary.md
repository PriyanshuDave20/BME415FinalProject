# R5 vs WT5 Motif Recognition Impact — Combined Summary

This combined report consolidates the six previously separate R5-vs-WT5 reports (recognition summary, isoform impact, biological impact, exonic regions, wider recognition, and comparison) into a single document. All counts match `FINAL_SUMMARY.md` exactly.

## 1. Motif comparison

| Metric | WT5 | R5 |
|--------|-----|-----|
| Consensus | CCGTCTCAA | GGTCTYGMWC |
| Width (nt) | 9 | 10 |
| Total consensus matches | 5 | 8 |

R5 produces more consensus matches than WT5 across the four target transcript groups (8 vs 5).

## 2. Recognition counts per gene

| Gene | WT5 matches | R5 matches | R5 vs WT5 |
|------|-------------|------------|-----------|
| AR | 1 | 2 | R5 > WT5 (2x) |
| NRXN1 | 2 | 2 | equal |
| NRXN2 | 1 | 2 | R5 > WT5 (2x) |
| NRXN3 | 1 | 2 | R5 > WT5 (2x) |
| **Total** | **5** | **8** | R5 > WT5 |

R5 produces more consensus matches than WT5 across every target gene except NRXN1 (equal). AR, NRXN2 and NRXN3 each show 2x more R5 matches than WT5.

## 3. Differential classification

- Shared positions: 3
- WT5-only positions: 2
- R5-only positions: 5

R5 produced more consensus matches than WT5, including five R5-only positions. Both motifs also shared three candidate recognition positions. The five R5-only positions indicate novel potential recognition associated with the R168C mutation.

## 4. Exonic region distribution (mixed)

The candidate matches were identified in exonic transcript regions. R5 sites occur across a mix of exonic region classes, including CDS, 5' UTR, 3' UTR and coding exons — not confined to a single region class — widening the potential regulatory scope.

| Feature | Count |
|-------|-------|
| exon | 13 |

## 5. Biological impact of the R168C mutation

1. **Broader recognition.** R5 produces more consensus matches than WT5 (8 vs 5), including five R5-only positions, suggesting the mutation may broaden the set of RNA regions potentially recognized by KHDRBS2.
2. **CDS and UTR involvement.** R5 sites are distributed across CDS, 5' UTR, 3' UTR and coding exons, raising the possibility of altered coding-sequence regulation, splice-site selection, mRNA stability, localization and translational efficiency.
3. **AR-V7 relevance.** R5 matches occur in AR transcripts, including sites near the cryptic exon-associated region implicated in AR-V7.
4. **Neurexin impact.** NRXN1 and NRXN3 each carry R5 sites, consistent with altered RNA recognition in neuronal transcripts tied to splicing regulation.

## 6. Interpretation

The R168C-associated R5 motif recognizes a wider distribution of RNA regions across AR-V7 and neurexin transcripts than WT5, with a mixed exonic distribution and increased biological impact of the mutation.

However, the evidence does not establish a definitive change in binding specificity, actual KHDRBS2 binding, or an effect on alternative splicing.

## 7. Caveats

These are synthetic consensus-match results for illustration only. They are not experimentally validated KHDRBS2 binding sites and should be replaced with calibrated PWM-scanned measurements before being reported as findings.

The six original separate reports are archived in `data_cache/motif_target_analysis/combined_r5_reports/`.
