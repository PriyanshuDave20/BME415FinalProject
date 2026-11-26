# Limitations

1. **Motif sample size**: Both WT5 and R5 were discovered from a limited number of sequences.
2. **Borderline similarity statistics**: TOMTOM q-values were around 0.05 in both directions.
3. **PWM uncertainty**: Neither motif was consistently recovered across all biological replicates.
4. **Consensus matching only**: Without FIMO or calibrated PWM thresholds, the scan used IUPAC-aware consensus matching, which is less sensitive than a proper PWM scan.
5. **No RNA structure prediction**: ViennaRNA was not available; only a GC-content heuristic was used.
6. **Transcript annotation**: Not all transcript isoforms may be captured in the current target_transcripts FASTA files.
7. **No experimental validation**: None of the candidate sites have been experimentally tested.
8. **Motif occurrence != binding**: Finding a motif sequence in a transcript does not prove that KHDRBS2 binds there.
