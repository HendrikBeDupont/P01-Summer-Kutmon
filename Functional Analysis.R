# Script:       FuncationalAnalysis.R
# Description:  In this script, we will explore a differential gene 
#               expression dataset comparing disease vs. control
#               samples, perform pathway enrichment analysis and 
#               build a PPI network
# Version: 1.0
# Last updated: 2025-02-03
# Author: mkutmon


# ##################################################################
# R INSTRUCTIONS
# ##################################################################

# * Lines that start with a # are comments
# * You can run a code line by placing the cursor in the line and clicking 
#   CTRL/Command + Enter

# ##################################################################
# R SETUP
# ##################################################################

# Here we install and load all required packages. 
if (!("BiocManager" %in% installed.packages())) { install.packages("BiocManager", update=FALSE) }
if (!("rstudioapi" %in% installed.packages())) { BiocManager::install("rstudioapi", update=FALSE) }
if (!("org.Hs.eg.db" %in% installed.packages())) { BiocManager::install("org.Hs.eg.db", update=FALSE) }
if (!("dplyr" %in% installed.packages())) { BiocManager::install("dplyr", update=FALSE) }
if (!("EnhancedVolcano" %in% installed.packages())) { BiocManager::install("EnhancedVolcano", update=FALSE) }
if (!("readxl" %in% installed.packages())) { BiocManager::install("readxl", update=FALSE) }
if (!("clusterProfiler" %in% installed.packages())) { BiocManager::install("clusterProfiler", update=FALSE) }
if (!("enrichplot" %in% installed.packages())) { BiocManager::install("enrichplot", update=FALSE) }
if (!("Rgraphviz" %in% installed.packages())) { BiocManager::install("Rgraphviz", update=FALSE) }
if (!("RCy3" %in% installed.packages())) { BiocManager::install("RCy3", update=FALSE) }
if (!("msigdbr" %in% installed.packages())) { BiocManager::install("msigdbr",update=FALSE) }
if (!("RColorBrewer" %in% installed.packages())) { BiocManager::install("RColorBrewer",update=FALSE) }
if (!("readr" %in% installed.packages())) { BiocManager::install("readr",update=FALSE) }

library(rstudioapi)
library(org.Hs.eg.db)
library(dplyr)
library(EnhancedVolcano)
library(readxl)
library(clusterProfiler)
library(enrichplot)
library(Rgraphviz)
library(RCy3)
library(msigdbr)
library(RColorBrewer)
library(readr)
library(ggplot2)

# We will set the working directory to the location where the current 
# script is located. This way, we can use relative file path locations. 
setwd(getwd())

# We will create an output folder where all figures and files will be stored
out.folder <- "output/"
dir.create(out.folder)

# ##################################################################
# 1.1 Data Import and Preprocessing
# ##################################################################

# Load dataset
data <- read_excel("GSE239914-differential-analysis.xlsx")

# Select specific columns (gene information, logFC, p-value)
data <- data[,c(8,1,6,2,3,10,11)]

# Filter for protein-coding genes and remove missing values
data.pc <- data %>%
  filter(GeneType == "protein-coding") %>%
  filter(!is.na(GeneID))

# Dropping a column 6 from the dataset 
data.pc <- data.pc[,-6]

# ##################################################################
# 1.2 Visualizing Results with a Volcano Plot
# ##################################################################

# Assign categories based on log2FC cutoffs
data.pc$Category <- "Not Significant"
data.pc$Category[data.pc$padj < 0.05 & abs(data.pc$log2FoldChange) > 0.26] <- "log2FC > 0.26"
data.pc$Category[data.pc$padj < 0.05 & abs(data.pc$log2FoldChange) > 0.58] <- "log2FC > 0.58"
data.pc$Category[data.pc$padj < 0.05 & abs(data.pc$log2FoldChange) > 1] <- "log2FC > 1"

# Convert Category to a factor
data.pc$Category <- factor(data.pc$Category, levels = c("Not Significant", "log2FC > 0.26", "log2FC > 0.58", "log2FC > 1"))

# Create a color vector that matches each gene row in data.pc
category_colors_map <- c(
  "Not Significant" = "grey",
  "log2FC > 0.26" = "blue",
  "log2FC > 0.58" = "purple",
  "log2FC > 1" = "red"
)

# Assign colors to each row based on Category
data.pc$Color <- category_colors_map[data.pc$Category]

# Ensure color vector is the same length as the data
stopifnot(length(data.pc$Color) == nrow(data.pc))

# Create enhanced volcano plot
p <- EnhancedVolcano(
  data.pc,
  title = "Volcano Plot with Multiple Log2FC Cutoffs (5384) DEGs",
  lab = data.pc$Symbol,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  labSize = 3,
  xlim = c(min(data.pc$log2FoldChange, na.rm = TRUE), max(data.pc$log2FoldChange, na.rm = TRUE)),
  ylim = c(0, max(-log10(data.pc$padj), na.rm = TRUE)),
  colCustom = data.pc$Color,  # Assign colors to genes
  legendPosition = "top",
  vline = c(-1, 1)  # Add both vertical dashed lines at log2FC = ±1
)

# Save the plot
filename <- paste0(out.folder, "volcano-plot-multi-cutoff.png")
png(filename, width = 2000, height = 1500, res = 150)
print(p)
dev.off()

# ##################################################################
# Decision: Choose log2FC Cutoff, Set confidence Cutoff
# ##################################################################

# Define thresholds for differential expression
log2fc.cutoff <- 1
pvalue.cutoff <- 0.05

# ##################################################################
# 1.3 Filtering for Differentially Expressed Genes (DEGs)
# ##################################################################

# Select differently expressed genes (DEGs) (genes meeting log2FC and p-value criteria)
degs <- data.pc[abs(data.pc$log2FoldChange) > log2fc.cutoff & data.pc$padj < pvalue.cutoff,]

# Save the DEGs to a file
write.table(degs, file=paste0(out.folder,"degs.tsv"), row.names = FALSE, sep="\t", quote = FALSE)

# ##################################################################
# Optional: Run analysis on up- or down-regulated genes
# ##################################################################

# Selects genes where log2FoldChange is greater than the threshold (1 by default)
genes.up <- nrow(degs[degs$log2FoldChange > log2fc.cutoff,])

# Selects genes where log2FoldChange is less than -1
genes.down <- nrow(degs[degs$log2FoldChange < -log2fc.cutoff,])

# Print results
cat("Number of Upregulated Genes: ", genes.up, "\n")
cat("Number of Downregulated Genes: ", genes.down, "\n")

# ##################################################################
# 2.1 Over-representation analysis (ORA)
# ##################################################################

# Fetch WikiPathways gene sets
genesets.wp <- msigdbr(species = "Homo sapiens", subcategory = "CP:WIKIPATHWAYS") %>% dplyr::select(gs_name, entrez_gene)

# Perform pathway enrichment analysis
res.wp <- clusterProfiler::enricher(degs$GeneID, TERM2GENE = genesets.wp, pAdjustMethod = "fdr", pvalueCutoff = 0.05, minGSSize = 5, maxGSSize = 400)
res.wp.df <- as.data.frame(res.wp)

# Save results 
write.table(res.wp.df, file=paste0(out.folder,"WP-Enrichment.txt"), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ##################################################################
# 2.2 Tree Plot Visualization 
# ##################################################################

# Tree plot visualization ->> pathways that are close together share a lot of genes
# you can do res.wp.up /down based on this code to see up and downregulated genes
res.wp.sim <- enrichplot::pairwise_termsim(res.wp)
treeplot(res.wp.sim, label_format = 0.3, showCategory = nrow(res.wp.df), cluster.params = list(label_words_n = 0))

# Save Plot
filename <- paste0(out.folder,"WP_Treeplot.png")
png(filename , width = 3000, height = 4000, res = 150)
plot(treeplot(res.wp.sim, label_format = 0.3, showCategory = nrow(res.wp.df), cluster.params = list(label_words_n = 0)))
dev.off()

# ##################################################################
# Optional: Pathway visualization 
# ##################################################################

# Check if Cytoscape is running - keep it open and check what is happening while you are 
# running the code
cytoscapePing()

# Check if WikiPathways app is installed
if(!"name: WikiPathways, version: 3.3.10, status: Installed" %in% RCy3::getInstalledApps()) {
  RCy3::installApp("WikiPathways")
}

# Open Pathway of interest - based on the res.wp.df, you can
# select pathways of interest
# Find the associated pathway identifier
# https://www.wikipathways.org/browse/table.html
# Make sure you select the ID of the human pathway
# Example: MAPK Cascade

pw.id <- "WP422"
RCy3::commandsRun(paste0('wikipathways import-as-pathway id=',pw.id)) 

toggleGraphicsDetails()

# load the data into Cytoscape (columns get added at the bottom)
loadTableData(data.pc, data.key.column = "EnsemblGeneID", table.key.column = "Ensembl")

# visualize the log2FC as a node fill color gradient
RCy3::setNodeColorMapping(table.column = 'log2FoldChange', mapping.type = 'c', table.column.values = c(-2,0,2), colors = paletteColorBrewerRdBu, default.color = '#FFFFFF', style.name = 'WikiPathways')

# Select significant genes and change border color
x <- RCy3::createColumnFilter('padj', 'padj', 0.05, "LESS_THAN")
RCy3::setNodeBorderColorBypass(x$nodes, new.colors = "#009900")
RCy3::setNodeBorderWidthBypass(x$nodes, new.sizes = 7)
RCy3::clearSelection()


# ##################################################################
# 3.1 Creation of the PPI Network
# ##################################################################

# make sure Cytoscape is running
RCy3::cytoscapePing()

if(!"name: stringApp, version: 2.0.3, status: Installed" %in% RCy3::getInstalledApps()) {
  RCy3::installApp("stringApp")
}

# Select strict DEGs (log2FC > 2)
degs.strict <- degs[abs(degs$log2FoldChange) > 2,]

# Convert gene list for STRING query
query <- format_csv(as.data.frame(degs.strict$Symbol), col_names=F, escape = "double", eol =",")

# Create STRING PPI network
commandsRun(paste0('string protein query cutoff=0.7 newNetName="PPI network" query="',query,'" limit=0 species="Homo sapiens"'))

# =======================
# if you run into this error: Error:  reason: URI Too Long
# uncomment (remove #) and run the following line - create the network "manually" in Cytoscape, then run the rest of the code normally (check video!!)

# windows:
# writeClipboard(query) 

# apple: 
install.packages("clipr")
library(clipr)

# Copy the gene list to clipboard
query <- paste(degs.strict$Symbol, collapse=",")
write_clip(query)
# =======================

# ##################################################################
# 3.2 Visualization and Analysis
# ##################################################################

# Analyze network topology
RCy3::analyzeNetwork()
hist(RCy3::getTableColumns(columns = "Degree")$Degree, breaks=100)

# network topology
RCy3::createVisualStyle("centrality")
RCy3::setNodeLabelMapping("display name", style.name = "centrality")
colors <-  c ('#FFFFFF', '#DD8855')
setNodeColorMapping("Degree", c(0,60), colors, style.name = "centrality", default.color = "#C0C0C0")
setNodeSizeMapping("Degree", table.column.values = c(0,60), sizes = c(30,100), mapping.type = "c", style.name = "centrality", default.size = 10)
RCy3::setVisualStyle("centrality")
RCy3::toggleGraphicsDetails()

# ##################################################################
# 3.3 Clustering and Drug-target Extension
# ##################################################################

# data vTRUEisualization clutering
RCy3::loadTableData(data=data.pc, data.key.column = "Symbol", table = "node", table.key.column = "query term")
RCy3::createVisualStyle("log2FC vis")
RCy3::setNodeLabelMapping("display name", style.name = "log2FC vis")
control.points <- c (-5.0, 0.0, 5.0)
colors <-  c ('#5588DD', '#FFFFFF', '#DD8855')
setNodeColorMapping("log2FoldChange", control.points, colors, style.name = "log2FC vis", default.color = "#C0C0C0")
RCy3::setVisualStyle("log2FC vis")
RCy3::lockNodeDimensions("TRUE", "log2FC vis")

###############DOESNT WORK################
# Drug traget extension, define the path to your LinkSet file (adjust this path to your actual file location)
linkset_path <- normalizePath("C:/Users/hendrikdupont/Documents/drugbank4-2.xgmml")

# Extend the network using CyTargetLinker
RCy3::commandsRun(paste('cytargetlinker extend',
                        'idAttribute="query term"',  
                        paste0('linkSetFiles="', linkset_path, '"'),
                        'network=current'))

# Optional: Visualize drug-target interactions with different node shapes
RCy3::setNodeShapeMapping("Type", c("Drug", "Protein"), c("diamond", "ellipse"), style.name = "log2FC vis")

# Optional: Highlight drug-target interactions using edge color
RCy3::setEdgeColorMapping("interaction", c("drugs-targets"), c("#FF5555"), style.name = "log2FC vis")
