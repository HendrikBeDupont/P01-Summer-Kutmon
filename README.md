# Functional Bioinformatics Analysis Pipeline

This repository contains an R-based pipeline for performing functional bioinformatics analysis on differential gene expression datasets. 
The pipeline is designed to process transcriptomic data, identify differentially expressed genes (DEGs), perform pathway enrichment analysis, 
and construct protein-protein interaction (PPI) networks.

The repository also contains a Python script rendering a diagraph of the workflow. This workflow is particularly suited for analyzing datasets like GSE239914 (vEDS Fibroblasts) but can be adapted to other transcriptomic datasets.

# Authors

Originally developed by [Dr. Summer-Kutmon](https://github.com/mkutmon) and later refined by [E. Retounioti](https://github.com/ErsiRet?tab=overview&from=2025-02-01&to=2025-02-20), [L. Cunningham](https://github.com/Lilabelle), [Q. Rosina](https://github.com/RabanvdBrand), [S. van Weeren](https://github.com/Spitvuurtje), [M. Reinartz](https://github.com/marvinreinartz), [R. van den Brand](https://github.com/RabanvdBrand), and [H. Dupont](https://github.com/HendrikBeDupont).

# Overview

This project aims to analyze **differentially expressed genes (DEGs)** by performing:  
✔ **Data Import and Preprocessing** of RNA-seq data  
✔ **Visualizing Results with a Volcano Plot** for DEG selection  
✔ **Over-representation analysis** using a case-specific gene set collection
✔ **Pathway visualization** 
✔ **Protein-Protein Interaction (PPI)** Network visualization and Analysis
✔ **Clustering and Drug-target extension** to identify key regulatory genes  

The pipeline is implemented in **R** and integrates **Cytoscape** for network visualization.
