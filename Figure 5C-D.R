#                           .....ILC NK HIPEC.....
#                                  Figure 5C-D

# This script creates figures present in Figure 5C-D of the manuscript
# Marchalot et al. "Tumor-infiltrating immature innate lymphoid cells in colorectal 
# cancer are biased towards tissue-resident NK cell/ILC1 differentiation"
# Prerequisites: Prerequisites: Load 00_SETUP.R script
# Input: SingleCellExperiment object created from FCS files
# Output: plots as in Figure 5C-D in Marchalot et al.


# Figure 5C UMAP of ILC and NK cells clusters in flow cytometry ----------------------------------------------------
# Load SingleCellExperiment object created from FCS files
sce <- readRDS("path/to/sce.rds") 
# UMAP plot colored by clusters
plotDR(sce, "UMAP", color_by = "merging1")
ggsave("Figure_5C.jpeg", width = 8, height = 6, units = "in", dpi = 400)

# Figure 5D Heatmap of cluster markers ----------------------------------------------------
plotExprHeatmap(sce, features = "state",
                by = "cluster_id", k = "merging1", bars = TRUE, perc = TRUE)
ggsave("Figure_5D.jpeg", width = 10, height = 8, units = "in", dpi = 400)
