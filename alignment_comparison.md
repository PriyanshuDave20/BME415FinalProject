# Alignment comparison: raw reads vs preprocessed reads

Old = minimap2 -ax sr on the raw deposited FASTQs (results/
alignments/, summary results/alignment_summary.tsv).
New = the corrected flow (module2_preprocess_reads.py: easyCLIP
3' adapter trimming, UMI moved into the read name, pairs < 14 nt
discarded) followed by this module on the preprocessed FASTQs
(results/alignments_preprocessed/, summary
results/alignment_summary_preprocessed.tsv).

IMPORTANT: the preprocessed input is slightly smaller than the
raw input (pairs < 14 nt after trimming are discarded, and pairs
with a malformed UMI are dropped), so the denominators differ.
A higher mapping percentage does not by itself mean a better
alignment - the fraction of the ORIGINAL read pool that ends up
mapped is the fairer number (reported as retained-then-mapped).

| sample | aligner | input pairs | retained pairs | mapped reads | mapped % | properly paired | secondary (multimap) |
|---|---|---|---|---|---|---|---|
| SRR12192531 | minimap2 | 699926 | 696608 | 572626 | 41.10 | 483935 | 0 |
| SRR12192532 | minimap2 | 6002088 | 5983916 | 4817716 | 40.26 | 4190189 | 0 |
| SRR12192533 | minimap2 | 1253463 | 1249647 | 1000210 | 40.02 | 849644 | 0 |
| SRR12192534 | minimap2 | 665911 | 662580 | 559051 | 42.19 | 479182 | 0 |
| SRR12192535 | minimap2 | 1068782 | 1065877 | 942240 | 44.20 | 811553 | 0 |
| SRR12192537 | minimap2 | 1244505 | 1239542 | 866148 | 34.94 | 731279 | 0 |
| SRR12192538 | minimap2 | 10000420 | 9969723 | 6798814 | 34.10 | 5914498 | 0 |
| SRR12192539 | minimap2 | 2141073 | 2133950 | 1437953 | 33.69 | 1224026 | 0 |
| SRR12192540 | minimap2 | 2243004 | 2235821 | 1464188 | 32.74 | 1261423 | 0 |
| SRR12192541 | minimap2 | 1446711 | 1438577 | 918703 | 31.93 | 791897 | 0 |

Differences between conditions and replicates should be
interpreted together with the preprocessing summary (adapter
content was 0.1-2.2% per sample and the retained fraction
99.4-99.6%).
