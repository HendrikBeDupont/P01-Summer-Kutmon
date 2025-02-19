# Functional Bioinformatics Analysis Pipeline

This repository contains an R-based pipeline for performing functional bioinformatics analysis on differential gene expression datasets. 
The pipeline is designed to process transcriptomic data, identify differentially expressed genes (DEGs), perform pathway enrichment analysis, 
and construct protein-protein interaction (PPI) networks.

The repository also contains a Python script rendering a diagraph of the workflow. This workflow is particularly suited for analyzing datasets like GSE239914 (vEDS Fibroblasts) but can be adapted to other transcriptomic datasets.

# Authors

Originally developed by [Dr. Summer-Kutmon](https://github.com/mkutmon) and later refined by E. Retounioti, L. Cunningham, Q. Rosina, S. van Weeren, M. Reinartz, R. van den Brand, and [H. Dupont](https://github.com/HendrikBeDupont).

# Overview

This project aims to analyze **differentially expressed genes (DEGs)** by performing:  
✔ **Preprocessing & Filtering** of RNA-seq data  
✔ **Volcano Plot Visualization** for DEG selection  
✔ **Pathway Enrichment Analysis** using various gene set databases  
✔ **Protein-Protein Interaction (PPI) Network Analysis**  
✔ **Clustering & Topology Analysis** to identify key regulatory genes  

The pipeline is implemented in **R** and integrates **Cytoscape** for network visualization.
