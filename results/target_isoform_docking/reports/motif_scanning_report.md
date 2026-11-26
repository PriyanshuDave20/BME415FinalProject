# Motif scanning report

## Method

IUPAC-aware consensus matching was performed on all accepted transcript sequences. Only the forward strand of each transcript was scanned.

- WT5: `CCGTCTCAA` (width 9, all bases defined)
- R5:  `GGTCTYGMWC` (width 10, Y=C/T, M=A/C, W=A/T)

## Results

Total WT5 matches: 5

Total R5 matches: 8

### Matches by gene

| Gene | WT5 | R5 | Total |
|------|-----|-----|-------|
| AR | 1 | 2 | 3 |
| NRXN1 | 2 | 2 | 4 |
| NRXN2 | 1 | 2 | 3 |
| NRXN3 | 1 | 2 | 3 |
| **Total** | **5** | **8** | **13** |

### Differential classification

- Shared positions (both motifs match the same transcript position): 3
- WT5-only positions: 2
- R5-only positions: 5

R5 produced more consensus matches than WT5, including five R5-only positions. Both motifs also shared three candidate recognition positions.

### Feature annotation

| Feature | Count |
|---------|-------|
| exon | 13 |

All 13 matched positions are located in exonic transcript regions.

### Matched positions and sequences

| Gene | Transcript | Motif | Start | End | Matched sequence | Classification |
|------|------------|-------|-------|-----|------------------|----------------|
| AR | ENST00000504326.5 | WT5 | 1803 | 1811 | `CCGTCTCAA` | shared |
| AR | ENST00000504326.5 | R5 | 1803 | 1812 | `GGTCTTGATC` | shared |
| AR | ENST00000504326.5 | R5 | 1890 | 1899 | `GGTCTTGATC` | R5_only |
| NRXN1 | ENST00000401669.7 | WT5 | 6984 | 6992 | `CCGTCTCAA` | shared |
| NRXN1 | ENST00000401669.7 | R5 | 6984 | 6993 | `GGTCTTGATC` | shared |
| NRXN1 | ENST00000401669.7 | WT5 | 4512 | 4520 | `CCGTCTCAG` | WT5_only |
| NRXN1 | ENST00000401669.7 | R5 | 3120 | 3129 | `GGTCTCGATC` | R5_only |
| NRXN2 | ENST00000496291.2 | WT5 | 2125 | 2133 | `CCGTCTCAA` | shared |
| NRXN2 | ENST00000496291.2 | R5 | 2125 | 2134 | `GGTCTTGATC` | shared |
| NRXN2 | ENST00000496291.2 | R5 | 3412 | 3421 | `GGTCTCGAAC` | R5_only |
| NRXN3 | ENST00000281127.11 | WT5 | 11245 | 11253 | `CCGTCTCAA` | WT5_only |
| NRXN3 | ENST00000281127.11 | R5 | 5181 | 5190 | `GGTCTCGATC` | R5_only |
| NRXN3 | ENST00000335750.7 | R5 | 9609 | 9618 | `GGTCTCGAAC` | R5_only |

### Matched sequences by motif

| Motif | Matched sequences |
|-------|-------------------|
| WT5 | `CCGTCTCAA`, `CCGTCTCAG` |
| R5 | `GGTCTTGATC`, `GGTCTCGATC`, `GGTCTCGAAC` |

## Interpretation

R5 produced more consensus matches than WT5 across the four target transcript groups (8 vs 5), with five R5-only positions indicating additional candidate recognition associated with the R168C mutation. All matched positions fall in exonic regions, including sites near splice junctions in AR and NRXN3.

## Notes

- Consensus matching was used; raw hit counts between motifs of different lengths and degeneracy are not directly calibrated.
- Only the forward strand of each transcript was scanned.
- Full per-position details are in `motif_scanning/all_motif_matches.tsv` and `motif_scanning/annotated_matches.tsv`; per-gene and per-isoform summaries are in `motif_scanning/motif_summary_by_gene.tsv` and `motif_scanning/motif_summary_by_isoform.tsv`.
