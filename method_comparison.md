# Region-identification method comparison

Method A (baseline): 4 kb bins counted at the fragment
leftmost aligned base, per-sample Poisson-tail threshold,
consecutive enriched bins merged (Module 3).

Method B (comparison): 250 bp bins counted at the read 5' end
(the easyCLIP crosslink-site proxy), same Poisson-tail
threshold and merging, run on the UMI-deduplicated
biological-replicate BAMs.

The original easyCLIP workflow calls peaks with CLIPper on
deduplicated BAMs. CLIPper could not be reproduced here (its
python2-era toolchain and dependency pins are impractical in
this environment), so Method B implements the crosslink-site
resolution principle with the same threshold logic as the
baseline. Neither method is a calibrated significance test:
the Poisson background ignores mappability and local
read-density structure.

| metric | WT | R168C |
|---|---|---|
| 4 kb regions (Method A) | 12322 | 14636 |
| 250 bp regions (Method B) | 14276 | 21385 |
| 250 bp regions overlapping a 4 kb region | 13817 (96.8%) | 20586 (96.3%) |
| 4 kb regions containing a 250 bp region | 5820 (47.2%) | 7933 (54.2%) |

Interpretation: the 4 kb regions are an order of magnitude
coarser. Most 250 bp crosslink-site regions fall inside 4 kb
regions, but only a fraction of 4 kb regions contain a
detected crosslink-site cluster - the 4 kb method calls
broader, weaker signal that the narrow method does not
resolve. For downstream motif and transcript work the
narrower regions are preferable; the 4 kb regions remain a
conservative whole-region screen. Per-replicate narrow
regions are in results/binding_regions/narrow_bins/.
