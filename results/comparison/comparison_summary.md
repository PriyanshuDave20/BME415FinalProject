# WT vs R168C candidate binding region comparison

Candidate regions from Module 3 (read-coverage based, NOT
confirmed binding sites) compared between KHDRBS2 WT and
R168C easyCLIP-seq.

## 1. Biological replicate verification

Each condition has 4 biological replicates, not 6. The MiSeq
runs are technical repeats of the same library as the HiSeq
runs (verified from the shared library prep IDs in the FASTQ
headers and molecule key overlap tests - same-library pairs
share 7-20% of molecule keys, different-library controls
<= 0.25%). All replicate-level analyses use the deeper HiSeq
run per library (see replicate_metadata.tsv).

- SRR12192531 / SRR12192530: confirmed: read overlap 19.8% vs controls <= 0.2%
- SRR12192533 / SRR12192532: confirmed: read overlap 11.9% vs controls <= 0.2%
- SRR12192537 / SRR12192536: confirmed: read overlap 11.9% vs controls <= 0.2%
- SRR12192539 / SRR12192538: confirmed: read overlap 7.2% vs controls <= 0.2%

## 2. Normalization

PCR duplicates were removed with umi_tools (the original
easyCLIP workflow's dedup step) after merging the technical
runs of each library, so region support counts unique
molecules. Signal = fragments per million counted fragments
(FPM), using Module 3's genome-wide fragment count as
denominator so numerator and denominator have the same unit
(the earlier RPM version divided by mapped READS, ~3x larger
because it counts both mates and supplementary alignments;
those baseline tables are kept in baseline_rpm_6runs/).
Condition signal = mean over the 4 biological replicates;
per-replicate values and mean/median/sd per region are in
region_signal_comparison.tsv. Fold change = log2((R168C + 1.0)
/ (WT + 1.0)).

## 3. Region overlap

Two regions overlap when they share at least 1 x 4000 bp bin
on the same chromosome (regions merged per condition first).

- WT condition-level regions: 12322 (shared 10374, WT-associated 1948)
- R168C condition-level regions: 14636 (shared 9811, R168C-associated 4825)

A region found in only one condition is NOT proof that
binding is absent in the other condition - it can simply
reflect sequencing depth and detection thresholds.

## 4. Binding signal

- regions with stronger R168C signal (log2FC >= 1.0): 9166
- regions with stronger WT signal (log2FC <= -1.0): 3467

Strongest differences (top 10 each):

Stronger in R168C:

- WT_12074 chrX:73788000-73792000  WT 4.1 FPM, R168C 592.5 FPM, log2FC 6.85
- WT_01209 chr1:225776000-225780000  WT 2.5 FPM, R168C 225.7 FPM, log2FC 6.02
- WT_01365 chr10:5724000-5728000  WT 3.5 FPM, R168C 208.3 FPM, log2FC 5.53
- WT_12240 chrX:124088000-124092000  WT 3.5 FPM, R168C 200.6 FPM, log2FC 5.48
- WT_06342 chr2:216528000-216532000  WT 2.7 FPM, R168C 161.4 FPM, log2FC 5.46
- WT_10172 chr6:154868000-154872000  WT 2.5 FPM, R168C 139.2 FPM, log2FC 5.33
- WT_01717 chr10:74464000-74468000  WT 3.3 FPM, R168C 142.0 FPM, log2FC 5.05
- R168C_06117 chr18:32068000-32072000  WT 0.0 FPM, R168C 30.3 FPM, log2FC 4.97
- WT_10170 chr6:154824000-154828000  WT 2.5 FPM, R168C 103.4 FPM, log2FC 4.91
- R168C_12652 chr7:106088000-106092000  WT 0.0 FPM, R168C 28.4 FPM, log2FC 4.88

Stronger in WT:

- R168C_11164 chr6:13708000-13712000  WT 234.1 FPM, R168C 2.6 FPM, log2FC -6.01
- R168C_03109 chr12:25196000-25200000  WT 120.9 FPM, R168C 2.1 FPM, log2FC -5.32
- R168C_11163 chr6:13696000-13700000  WT 240.2 FPM, R168C 7.0 FPM, log2FC -4.91
- R168C_14020 chr9:111668000-111672000  WT 199.0 FPM, R168C 6.3 FPM, log2FC -4.78
- R168C_00402 chr1:68404000-68408000  WT 146.7 FPM, R168C 4.8 FPM, log2FC -4.66
- R168C_00368 chr1:63376000-63380000  WT 228.3 FPM, R168C 9.0 FPM, log2FC -4.52
- R168C_07336 chr2:183256000-183260000  WT 178.4 FPM, R168C 8.0 FPM, log2FC -4.31
- WT_08550 chr4:180580000-180588000  WT 18.8 FPM, R168C 0.0 FPM, log2FC -4.31
- WT_05848 chr2:146420000-146428000  WT 18.5 FPM, R168C 0.0 FPM, log2FC -4.28
- WT_10846 chr7:135128000-135144000  WT 17.6 FPM, R168C 0.0 FPM, log2FC -4.22

## 5. K-mer composition (k = 4, 5, 6)

K-mer frequencies over the condition-level region sequences
(+1 pseudocount per possible k-mer). The two conditions are
nearly identical in composition: the largest difference is
|log2 ratio| = 0.14 (kmer_comparison.tsv). A frequency
difference does not by itself prove a change in the KHDRBS2
binding motif - transcript abundance and the genomic
distribution of the detected regions also shape it.

## 6. Downsampling sensitivity

The deeper run of each matched replicate pair was downsampled
to its partner's counted-fragment depth (deterministic
read-name hash, seed 154168, mates kept together; no samples
discarded) and the whole comparison was repeated:

| category | original | downsampled |
|---|---|---|
| WT regions | 12322 | 12337 |
| R168C regions | 14636 | 12468 |
| shared (WT side) | 10374 | 9315 |
| WT-associated | 1948 | 3022 |
| shared (R168C side) | 9811 | 9263 |
| R168C-associated | 4825 | 3205 |

Per-region log2 fold changes correlate with Pearson r = 0.84
over 24915 matched regions (see downsampling_fc_scatter.png).
If the R168C-associated excess
shrinks after downsampling, it was largely a sequencing-depth
effect; if it persists, depth alone does not explain it.

## 7. Caveats

- Candidate regions are inferred from read coverage, not
  experimentally confirmed binding sites.
- FPM normalization reduces but does not remove detection
  bias; downsampling is a sensitivity analysis, not proof
  that depth bias is eliminated.
- Descriptive fold changes are not statistical significance;
  no statistical test was applied.
- K-mer frequency differences do not establish changes in
  binding specificity (motif discovery is Module 5).
