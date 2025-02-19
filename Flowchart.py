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

# First Decision Point: Cutoff
workflow.node("Cutoff", "Decision: Choose Log2FC Cutoff", shape="diamond", fillcolor="lightyellow")

# Path A: Gene Sets Enrichment Analysis
workflow.node("GeneSets", "Decision: Choose Gene Set Collection", shape="diamond", fillcolor="lightyellow")
workflow.node("ORA", "2.1 Over-Representation Analysis (ORA)", fillcolor="lightgray")
workflow.node("TreePlot", "2.2 Tree Plot (Visualization of ORA)", fillcolor="lightgray")
workflow.node("PathwayVisualization", "2.3 Optional: Pathway Visualization", fillcolor="lightgray")


# Path B: Network-Based Analysis
workflow.node("Confidence", "Decision: Set p-value Cutoff", shape="diamond", fillcolor="lightyellow")
workflow.node("PPI", "3.1 Protein-Protein Interaction (PPI) Network", fillcolor="lightgreen")
workflow.node("Clustering", "3.2.1 Clustering: Finding Communities", fillcolor="lightgreen")
workflow.node("Topology", "3.2.2 Network Topology (Hubs & Bottlenecks)", fillcolor="lightgreen")

# Final Report
workflow.node("Report", " 4 Final Report: Functional Interpretation", fillcolor="lightblue")

# Define workflow edges
workflow.edge("Dataset", "Preprocessing")
workflow.edge("Preprocessing", "Volcano")
workflow.edge("Volcano", "DEGs")
workflow.edge("DEGs", "Cutoff")

# Path A: Gene Sets Enrichment Analysis
workflow.edge("Cutoff", "GeneSets")
workflow.edge("GeneSets", "ORA")
workflow.edge("ORA", "PathwayVisualization")
workflow.edge("ORA", "TreePlot")
workflow.edge("TreePlot", "PathwayVisualization")
workflow.edge("TreePlot", "Report")

# Path B: Network-Based Analysis
workflow.edge("Cutoff", "Confidence")
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
    legend.node("Legend_Blue", "General Data Processing & Results", fillcolor="lightblue", shape="box")
    legend.node("Legend_Yellow", "Decision Points (User Choice Required)", fillcolor="lightyellow", shape="box")
    legend.node("Legend_Gray", "Pathway Enrichment Analysis", fillcolor="lightgray", shape="box")
    legend.node("Legend_Green", "Network Analysis (Clustering & Topology)", fillcolor="lightgreen", shape="box")

    # Arrange legend items vertically
    legend.edge("Legend_Blue", "Legend_Yellow", style="invis")
    legend.edge("Legend_Yellow", "Legend_Gray", style="invis")
    legend.edge("Legend_Gray", "Legend_Green", style="invis")

# Render and save flowchart with dataset node and vertical legend
workflow.render("functional_bioinformatics_workflow_final_v2")

print("Flowchart with legend saved as functional_bioinformatics_workflow_final_v2.png")
