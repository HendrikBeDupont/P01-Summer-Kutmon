# Functional Bioinformatics Analysis Pipeline

This repository contains an R-based pipeline for performing functional bioinformatics analysis on differential gene expression datasets. 
The pipeline is designed to process transcriptomic data, identify differentially expressed genes (DEGs), perform pathway enrichment analysis, 
and construct protein-protein interaction (PPI) networks.

This workflow is particularly suited for analyzing datasets like GSE239914 (vEDS Fibroblasts) but can be adapted to other transcriptomic datasets.

# Authors

Originally developed by Dr. Summer-Kutmon and later refined by E. Retounioti, L. Cunningham, Q. Rosina, S. van Weeren, M. Reinartz, R. van den Brand, and H. Dupont.

# Overview

This pipeline follows a structured analysis workflow for transcriptomic data, incorporating key steps:

- Preprocessing & DEG Identification
  
Import structured gene expression data.
Filter and format data for analysis.
Identify DEGs based on log2 fold change and adjusted p-values.

- Visualization
  
Generate volcano plots to visualize DEGs at different cutoffs (log2FC = 0.26, 0.58, 1).
Create tree plots to show pathway relationships.

- Pathway Enrichment Analysis
  
Perform Over-Representation Analysis (ORA) using either:
Gene Ontology (GO)
KEGG Pathways
WikiPathways
MSigDB
Reactome
Evaluate functional relevance of significantly altered pathways.

- Protein-Protein Interaction (PPI) Network

Construct PPI networks using STRING and visualize them in Cytoscape.
Perform network topology analysis to identify hub genes.
Cluster genes into functional communities.

- Final Report
  
Interpret results based on enrichment analysis and network insights.

# Installation & Dependencies
To run the pipeline, install R (≥4.0) and the following R packages:

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

Install required Bioconductor packages
BiocManager::install(c("org.Hs.eg.db", "clusterProfiler", "enrichplot", "msigdbr", "RCy3"))

# Install CRAN packages
install.packages(c("ggplot2", "dplyr", "readxl", "EnhancedVolcano", "RColorBrewer"))
Additionally, install Cytoscape and the STRING plugin for network analysis.

# Usage
Run the script FunctionalAnalysis.R in RStudio or via command line:
source("FunctionalAnalysis.R")

- This will:

Process the dataset.
Generate volcano plots for differential expression.
Perform pathway enrichment analysis.
Construct PPI networks in Cytoscape.
Generate a final report.

# Example: Adjusting Cutoff Values
Modify the log2fc.cutoff parameter to filter DEGs with different fold-change thresholds:

log2fc.cutoff <- 1 # Change to 0.58 or 0.26 if needed

# Output
- The pipeline generates:
  
DEG tables (degs.tsv)
Volcano plots (volcano-plot.png)
Pathway enrichment tables (WP-Enrichment.txt)
Tree plots (WP_Treeplot.png)
PPI network files (for Cytoscape)
Interpretation Guide
