# Colorectal Cancer Study of Austria - PhIP-Seq analysis

This repository contains the analysis code used to process and visualize PhIP-Seq data for the colorectal cancer project.

The analysis is based on the `Phipseq_Sampleset` Python package, which is used to load the PhIP-Seq datasets, handle sample metadata, generate presence/absence and summary tables, and organize the data for downstream analyses.

Paper figures and exploratory analyses are produced in the accompanying Jupyter notebook. The notebook contains the plotting functions used to generate the manuscript figures, including repertoire-level analyses, prevalence-based filtering, antigen-level summaries, machine-learning outputs, and selected candidate antigen visualizations.

## Basic workflow

1. Load PhIP-Seq data and metadata using `Phipseq_Sampleset`.
2. Generate antigen presence/absence and annotation tables.
3. Run downstream analyses in the Jupyter notebook.
4. Export plots and tables used in the manuscript.

## Data analysis
Main data analysis were made through the Jupyter Notebook `Analysis_CORSA_Figures_and_tables.ipynb`, so the inspection of that is enough for observing all the results, including tables and figures. This loads all the relevant tables, and has all the functions generating the outputs.
Since we stored the data in a centralized database with all the generated data from our lab, we cannot provide access to the full DB, or the libraries full data (available upon request). 
We offer then the presence/absence and foldchange files which show for all the patients if a given oligopeptide was enriched in `data` folder and the metadata can be loaded in `metadata` folder.

In depth analysis were performed with 2 additional packages:
1. phiper v 0.4.1 (https://github.com/Polymerase3/phiper) - all relevant scripts for the data produced are in `phiper_analysis`
2. phipml (https://github.com/csReynaB/phipml) - all the relevant scripts found in `predictions/` folder (in the notebook the path to the joblib file may need to be adjusted)

## Notes

Input raw data are not included in this repository unless explicitly stated, but they will be uploaded to Dataverse server with access for reviewers. Since the analysis was performed locally, paths to local data files should be adjusted in the notebook before rerunning the analysis.
