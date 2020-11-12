# Stored Dairy Manure Microbiome Workflow

In this small pilot study I examine drivers of microbial community and antimicrobial resistance gene abundances in stored dairy manure. GROOT and AMRPlusPlus were used for read assignment to Antimicrobial Resistance (AMR) gene sequences and Phyloflash was used to assign reads to 16S rRNA gene sequences. The files listed below run different statistical analysis to and perform some preliminary plotting.

# Files
2020-09_Procrustes.ipynb: Here I evaluate how patters of organism and AMR gene abundances compare across farms using the Procrustes analysis. I am beginning to build code to run different user-specified dimensional reduction techniques and plot the results. 

2020-09_CONCOR_AMR_annotated.ipynb: This python notebook explores different data-preprocessing prior to running Canonical Correlation Analysis (CONCOR). The purpose of CONCOR is evaluate how metadata parameters correlate with AMR gene abundance patterns. 

ccabiplot_topn.py: This python script takes the output from CONCOR analysis and generates a plot. This script was modified from code from Kerby Shedden to enable the user to set the number of arrows plotted with the highest magnitudes.

summarize_taxa.py: This python script was made to make it easier to summarize data at different taxanomic levels to assess how robust between-sample differences are to the resolution of taxonomic assignment
