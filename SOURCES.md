# Sources and documentation used for implementation

This file maps the likely and authoritative sources behind the project workflow. It is not a claim that every source was used for every line of code, but it identifies the main documentation and methods that support the implementation.

## General language, runtime, and workflow

- **Python 3** documentation for file handling, command-line interfaces, string processing, and standard-library usage.
- **Git** and **GitHub** documentation for version control, branching, and repository organization.
- **Windows Subsystem for Linux (WSL)** documentation for running command-line bioinformatics tools from Windows.
- **Conda / pip** documentation for managing Python environments and dependencies.

## Read preprocessing and alignment

- **cutadapt** documentation and paper for adapter trimming, quality filtering, and read preprocessing.
- **UMI-tools** documentation and paper for UMI extraction and deduplication.
- **minimap2** documentation and paper for short-read alignment to the GRCh38 genome.
- **SAMtools** documentation and paper for BAM manipulation, sorting, indexing, and alignment summaries.
- **pysam** documentation for programmatic access to BAM and SAM files.

## Binding-region discovery and motif analysis

- **BEDTools** documentation for genomic interval arithmetic and overlap-based region analysis.
- **pybedtools** documentation when used as a Python interface to BEDTools.
- **MEME Suite** documentation and papers for motif discovery and motif scanning.
- **STREME** documentation for enriched motif discovery.
- **FIMO** documentation for motif scanning against transcript sequences.
- **IUPAC nucleotide ambiguity codes** documentation for interpreting motif sequences.

## Transcript and annotation resources

- **GENCODE** documentation and release notes for transcript annotations, exon coordinates, and transcript identifiers.
- **Ensembl** documentation for genome assembly and transcript metadata.
- **NCBI RefSeq** documentation for transcript and gene identifiers.
- **GTF/GFF format documentation** for parsing exon, transcript, and gene feature files.

## Protein structure and docking

- **Biopython** documentation and paper for parsing and manipulating PDB/mmCIF structures.
- **RCSB PDB** documentation for PDB and mmCIF file formats and structural metadata.
- **AlphaFold2** documentation and paper for predicted protein structure interpretation.
- **AlphaFold Protein Structure Database** documentation for model metadata and confidence measures.
- **UniProt** documentation for protein sequence identifiers and functional annotations.
- **HDOCK** documentation and paper for protein-RNA docking and HDOCKlite usage.

## Post-docking analysis and visualization

- **Biopython** structural-analysis modules for residue parsing and distance calculations.
- **SciPy** documentation for spatial distance calculations and statistical summaries.
- **NumPy** documentation and paper for numerical operations.
- **pandas** documentation for tabular result assembly and reporting.
- **Matplotlib** documentation and paper for plotting and figure generation.

## Reproducibility and infrastructure

- Official documentation for **Snakemake** or **Nextflow** when pipeline orchestration is used.
- Official documentation for **Docker** or **Singularity** when containerized reproducibility is used.
- Official documentation for **GitHub Actions** when automated testing or workflow validation is used.

## Community troubleshooting and code-level references

- **Stack Overflow**
- **Bioinformatics Stack Exchange**
- **Biostars**
- **GitHub Issues**

These are useful for tool-specific edge cases, installation problems, file-format parsing, and command-line behavior. They are best treated as implementation aids rather than primary methodological references.

## References

1. Bailey, T. L., Boden, M., Buske, F. A., Frith, M., Grant, C. E., Clementi, L., Ren, J., Li, W. W., & Noble, W. S. (2009). MEME SUITE: tools for motif discovery and searching. *Nucleic Acids Research*, 37(W1), W202-W208. https://doi.org/10.1093/nar/gkp335
2. Berman, H. M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T. N., Weissig, H., Shindyalov, I. N., & Bourne, P. E. (2000). The Protein Data Bank. *Nucleic Acids Research*, 28(1), 235-242. https://doi.org/10.1093/nar/28.1.235
3. Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. L. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422-1423. https://doi.org/10.1093/bioinformatics/btp163
4. Frankish, A., et al. (2021). GENCODE 2021. *Nucleic Acids Research*, 49(D1), D916-D923. https://doi.org/10.1093/nar/gkaa1087
5. Harris, C. R., Millman, K. J., van der Walt, S. J., et al. (2020). Array programming with NumPy. *Nature*, 585(7825), 357-362. https://doi.org/10.1038/s41586-020-2649-2
6. Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9(3), 90-95. https://doi.org/10.1109/MCSE.2007.55
7. Jumper, J., Evans, R., Pritzel, A., et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature*, 596(7873), 583-589. https://doi.org/10.1038/s41586-021-03819-2
8. Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094-3100. https://doi.org/10.1093/bioinformatics/bty191
9. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078-2079. https://doi.org/10.1093/bioinformatics/btp352
10. Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*, 17(1), 10-12. https://doi.org/10.14806/ej.17.1.200
11. McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 51-56. https://doi.org/10.25080/Majora-92bf1922-00a
12. Smith, T., Heger, A., & Sudbery, I. (2017). UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. *Genome Research*, 27(3), 491-499. https://doi.org/10.1101/gr.209601.116
13. Virtanen, P., Gommers, R., Oliphant, T. E., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17(3), 261-272. https://doi.org/10.1038/s41592-019-0686-2
14. Yan, Y., Zhang, D., Zhou, P., Li, B., & Huang, S.-Y. (2017). HDOCK: a web server for protein-protein and protein-DNA/RNA docking based on a hybrid strategy. *Nucleic Acids Research*, 45(W1), W365-W373. https://doi.org/10.1093/nar/gkx407
