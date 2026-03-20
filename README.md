# **An AI-ready phosphorylation meta-analysis for _Saccharomyces cerevisiae_**


This analysis is based on a reanalysis of yeast phosphorylation data sets (PRIDE accession: PXD071918) which used [alanine as a decoy](https://github.com/PGB-LIV/mzidFLR) to control for the false localisation rate.

This repo contains a sub-folder for each element of the downstream analysis:

[01_motif_and_pathway_analysis/scripts](https://github.com/EllenBoswell/yeast_phospho/tree/main/01_motif_and_pathway_analysis/scripts) - generation of foreground and background 15mers, motif analysis (rmotifx) and gene ontology pathway enrichment analysis (clusterProfiler).

[03_disorder_analysis/scripts](https://github.com/EllenBoswell/yeast_phospho/tree/main/03_disorder_analysis/scripts) - processing of metapredict (V3) disorder predictions. The human analysis is based on data from [Kalyuzhnyy et al., 2025](https://pubs.acs.org/doi/10.1021/acs.jproteome.5c00278).

[04_AF3/scripts](https://github.com/EllenBoswell/yeast_phospho/tree/main/04_AF3/scripts) - generation of JSON input files for AlphaFold 3 and processing of AlphaFold models after DSSP.

[05_summary/scripts](https://github.com/EllenBoswell/yeast_phospho/tree/main/05_summary/scripts) - Summary figure for the yeast phosphobuild (Figure 1 in manuscript).
