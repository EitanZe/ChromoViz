Description: The contents of this repo are the work of the LOY on HSPCs group from the 2023 Weizmann Institute Bessie F. Lawrence ISSI  research program. It includes a script for processing data from Anndata files and two scripts for visualization.

Purpose: The scripts can be used in conjuction to analyze the expression of the genes associated with a certain chromosome using 10x single cell RNA sequencing data and plot it against other features such as age, sex, red blood cell count, etc. It can also be used to identify chromosomal aberrations and their effects. The unique benefits of single cell RNA sequencing is in its ability to differentiate between cell types and to show the unique molecular identifier levels of genes for individual cells.

Instructions: From the data processing file, users can export csvs with scores per cell and per persion about the gene expression levels given a certain gene list. From that, the visualizer files can be used to generate different graphs analyzing the data.

visualizer.py is meant to be used as a library to be used in conjuction with scripts like data processing.ipynb to eliminate the need to export to csv then load csv, but can also be run as a standalone program.



Credits:
visualizer.py - Eitan Zemel
data processor.ipynb - Shane Weinberger
austin
maia