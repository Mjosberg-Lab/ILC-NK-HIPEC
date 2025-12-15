#                           .....ILC NK HIPEC.....
#                                  Figure 2

# This script creates figures present in Figure 2A-D of the manuscript and prepares data for figure 2E-F
# Marchalot et al. "Tumor-infiltrating immature innate lymphoid cells in colorectal 
# cancer are biased towards tissue-resident NK cell/ILC1 differentiation"
# Prerequisites: Load 00_SETUP.R script
# Input: ILC_NK Seurat object available in GSE302045
# Output: plots as in Figure 2 in Marchalot et al.

ILC_NK@active.ident <- ILC_NK$clusters
# Figure 2A Stacked Violins Immature subsets ----------------------------------------------------
Common_top10 <- c('TCF7','SELL','GPR183','CD55','IFITM3','XCL1','XCL2')
jpeg(filename = "Figure_2A.jpeg", width = 10, height = 8, units = "in",res=400)
Stacked_VlnPlot(ILC_NK, features = Common_top10, x_lab_rotate = TRUE, colors_use = hue_pal()(14))
dev.off()

# Figure 2B Stacked Violins early NK ----------------------------------------------------
eNK_top <- c('IL12RB2','CRTAM','KLRC1','IRAK3','GNLY','COTL1','TNFSF9','RNF130')
jpeg(filename = "Figure_2B.jpeg", width = 10, height = 8, units = "in",res=400)
Stacked_VlnPlot(ILC_NK, features = eNK_top, x_lab_rotate = TRUE, colors_use = hue_pal()(14))
dev.off()

# Figure 2C Stacked Violins naive ILC ----------------------------------------------------
nILC_top <- c('TNFRSF4','MFGE8','TNFRSF18','CDKN1A','KIT','LTB','ICOS','TCF7','TIMP1','SCN1B','FXYD7','SSBP2')
jpeg(filename = "Figure_2C.jpeg", width = 10, height = 8, units = "in",res=400)
Stacked_VlnPlot(ILC_NK, features = nILC_top, x_lab_rotate = TRUE, colors_use = hue_pal()(14))
dev.off()

# Figure 2D Monocle-3 trajectory analysis and plot ----------------------------------------------------
# Transform Seurat object in cell_data_set format for Monocle
Mono_ILC_NK <- ILC_NK
Mono_ILC_NK$seurat_clusters <- Mono_ILC_NK@active.ident
Mono_ILC_NK <- as.cell_data_set(ILC_NK)
# Get feature/gene metadata
fData(Mono_ILC_NK)
fData(Mono_ILC_NK)$gene_short_name <- rownames(fData(Mono_ILC_NK))
head(fData(Mono_ILC_NK))
# Retrieve clustering information from Seurat object
recreate.partitions <- c(rep(1, length(Mono_ILC_NK@colData@rownames)))
names(recreate.partitions) <- Mono_ILC_NK@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions
Mono_ILC_NK@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
# Assign cluster info
list.cluster <- ILC_NK@active.ident
Mono_ILC_NK@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
# Assign UMAP coordinates
Mono_ILC_NK@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- ILC_NK@reductions$umap@cell.embeddings
# Learn trajectory
Mono_ILC_NK <- learn_graph(Mono_ILC_NK, use_partition = F)
# Order cells in Pseudotime
Mono_ILC_NK <- order_cells(Mono_ILC_NK, reduction_method = "UMAP", root_cells = colnames(Mono_ILC_NK[, clusters(Mono_ILC_NK) == 'early NK']))
# Order cells by Monocle 3 pseudotime
Mono_ILC_NK$monocle3_pseudotime <- pseudotime(Mono_ILC_NK)
data.pseudo <- as.data.frame(colData(Mono_ILC_NK))
jpeg(filename = "Figure_2D.jpeg", width = 6, height = 5, units = "in",res=400)
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot() + NoLegend()
dev.off()

# Figure 2E-F Export Seurat object to H5ad for Velocity analysis in Python ----------------------------------------------------
ILC_NK$cell_type <- ILC_NK@active.ident
ILC_NK$cell_type <- as.character(ILC_NK$cell_type)
path <- "/path/to/seurat_object_ILC_NK"
SaveH5Seurat(ILC_NK, filename = path, overwrite = T)
Convert(source =  paste(path,'.h5seurat',sep = ""), assay ="RNA", dest = "h5ad", overwrite = T)