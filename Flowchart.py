from graphviz import Digraph

# Create a Digraph object for the workflow
workflow = Digraph(format='png')

# Define nodes with labels and colors
workflow.attr("node", shape="box", style="filled")

# Dataset Node
workflow.node("Dataset", "Dataset: Transcriptome analysis comparing disease vs. control\nExample: GSE239914 (vEDS Fibroblasts)", fillcolor="white", shape="parallelogram")

# Preprocessing Step
workflow.node("Preprocessing", "1.1 Data Import and Preprocessing", fillcolor="lightblue")

# Volcano Plot (Helps decide cutoff)
workflow.node("Volcano", "1.2 Visualizing Results with a Volcano Plot", fillcolor="lightblue")

# Filtering for Differentially Expressed Genes (DEGs)
workflow.node("DEGs", "1.3 Filtering for Differentially Expressed Genes (DEGs)", fillcolor="lightblue")

# Run analysis on up- or down-regulated genes
workflow.node("UpDown", "Optional: Run analysis on up- or down-regulated genes", fillcolor="lightblue")

# First Decision Point: Cutoff
workflow.node("Cutoff", "Decision: Choose Log2FC Cutoff and Set p adjusted Cutoff", shape="diamond", fillcolor="lightyellow")

# Second Decision Point: Confidence
workflow.node("Confidence", "Decision: Set Confidence Cutoff", shape="diamond", fillcolor="lightyellow")


# Path A: Gene Sets Enrichment Analysis
workflow.node("GeneSets", "Decision: Choose Gene Set Collection", shape="diamond", fillcolor="lightyellow")
workflow.node("ORA", "2.1 Over-Representation Analysis (ORA)", fillcolor="crimson")
workflow.node("TreePlot", "2.2 Treeplot Visualization", fillcolor="crimson")
workflow.node("PathwayVisualization", "Optional: Pathway Visualization", fillcolor="crimson")


# Path B: Network-Based Analysis
workflow.node("PPI", "3.1 Creation of the PPI Network", fillcolor="lightgreen")
workflow.node("Clustering", "3.2 Visualization and Analysis", fillcolor="lightgreen")
workflow.node("Topology", "3.3 Clustering and Drug-target extension", fillcolor="lightgreen")

# Final Report
workflow.node("Report", " 1.4 Final Report: Functional Interpretation", fillcolor="lightblue")

# Define workflow edges
workflow.edge("Dataset", "Preprocessing")
workflow.edge("Preprocessing", "Volcano")
workflow.edge("Volcano", "Cutoff")
workflow.edge("Cutoff", "DEGs")
workflow.edge("DEGs", "UpDown")

# Path A: Gene Sets Enrichment Analysis
workflow.edge("DEGs", "GeneSets")
workflow.edge("GeneSets", "ORA")
workflow.edge("ORA", "PathwayVisualization")
workflow.edge("ORA", "TreePlot")
workflow.edge("TreePlot", "PathwayVisualization")
workflow.edge("TreePlot", "Report")

# Path B: Network-Based Analysis
workflow.edge("DEGs", "Confidence")
workflow.edge("Confidence", "PPI")
workflow.edge("PPI", "Clustering")
workflow.edge("PPI", "Topology")

# Clustering leads to database choice or final report
workflow.edge("Clustering", "GeneSets")
workflow.edge("Clustering", "Report")

# Topology leads directly to report
workflow.edge("Topology", "Report")

# Create a vertical legend in the top right corner using a subgraph cluster
with workflow.subgraph(name="cluster_legend") as legend:
    legend.attr(label="Legend", style="filled", color="white", rank="sink")
    legend.node("Legend_Blue", "Data Exploration & Results", fillcolor="lightblue", shape="box")
    legend.node("Legend_Yellow", "Decision Points (User Choice Required)", fillcolor="lightyellow", shape="box")
    legend.node("Legend_Red", "Pathway Enrichment Analysis", fillcolor="crimson", shape="box")
    legend.node("Legend_Green", "Network Analysis and Visualization", fillcolor="lightgreen", shape="box")

    # Arrange legend items vertically
    legend.edge("Legend_Blue", "Legend_Yellow", style="invis")
    legend.edge("Legend_Yellow", "Legend_Red", style="invis")
    legend.edge("Legend_Red", "Legend_Green", style="invis")

# Render and save flowchart with dataset node and vertical legend
workflow.render("functional_bioinformatics_workflow_final_v2")

print("Flowchart with legend saved as functional_bioinformatics_workflow_final_v3.png")
