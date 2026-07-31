# Colorectal Cancer Study of Austria - PhIP-Seq analysis

This repository contains the analysis code used to process and visualize PhIP-Seq data for the colorectal cancer project.

The analysis is based on the `Phipseq_Sampleset` Python package, which is used to load the PhIP-Seq datasets, handle sample metadata, generate presence/absence and summary tables, and organize the data for downstream analyses.

Paper figures and exploratory analyses are produced in the accompanying Jupyter notebook. The notebook contains the plotting functions used to generate the manuscript figures, including repertoire-level analyses, prevalence-based filtering, antigen-level summaries, machine-learning outputs, and selected candidate antigen visualizations.

## Basic workflow

1. Load PhIP-Seq data and metadata using `Phipseq_Sampleset`.
2. Generate antigen presence/absence and annotation tables.
3. Run downstream analyses in the Jupyter notebook.
4. Export plots and tables used in the manuscript.

## Notes

Input data are not included in this repository unless explicitly stated. Paths to local data files should be adjusted in the notebook before running the analysis.
