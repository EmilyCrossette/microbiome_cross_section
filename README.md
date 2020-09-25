# Stored Dairy Manure Microbiome Workflow

In this small pilot study I examine drivers of microbial community and antimicrobial resistance gene abundances in stored dairy manure. After using GROOT and AMRPlusPlus for read assignment to Antimicrobial Resistance (AMR) gene sequences and Phyloflash to assign reads to 

# Files
2020-09_CONCOR_AMR_annotated.ipynb: This python notebook explores different data-preprocessing prior to running Canonical Correlation Analysis (CONCOR). The purpose of CONCOR is evaluate how metadata parameters correlate with AMR gene abundance patterns. 

ccabiplot_topn.py: This python script takes the output from CONCOR analysis and generates a plot. This script was modified from code from Kerby Shedden to enable the user to set the number of arrows plotted with the highest magnitudes.

2020-09_Procrustes.ipynb: Here I evaluate how patters of organism and AMR gene abundances compare across farms using the Procrustes analysis. I am beginning to build code to run different user-specified dimensional reduction techniques and plot the results. 
