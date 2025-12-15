#                           .....ILC NK HIPEC.....
#                                  Figure 1

# This script creates figures present in Figure 1 of the manuscript
# Marchalot et al. "Tumor-infiltrating immature innate lymphoid cells in colorectal 
# cancer are biased towards tissue-resident NK cell/ILC1 differentiation"
# Prerequisites: Load 00_SETUP.R script
# Input: ILC_NK Seurat object available in GSE302045
# Output: plots as in Figure 1 in Marchalot et al.

# Figure 1B Clustering Tree Plot ----------------------------------------------------
ILC_NK@active.ident <- ILC_NK$seurat.clusters
ILC_NK <- BuildClusterTree(ILC_NK, reduction = "harmony")
jpeg(filename = 'Figure_1B.jpeg', width = 8, height = 8, units = "in",res=400)
PlotClusterTree(ILC_NK, font = 1, type = 'phylogram', srt = 90, adj = 0.5, 
                edge.color = 'black', tip.color	= 'black', node.color	= 'black')
dev.off()

# Figure 1C UMAP with clusters ----------------------------------------------------
ILC_NK@active.ident <- ILC_NK$seurat.clusters
jpeg(filename = 'Figure_1C.jpeg', width = 8, height = 8, units = "in",res=400)
DimPlot(ILC_NK, reduction = "umap", label = TRUE, label.size = 5) + NoLegend()
dev.off()

# Figure 1D scCostumize DotPlot  ----------------------------------------------------
all_markers <- FindAllMarkers(ILC_NK)
top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 5, named_vector = FALSE, make_unique = TRUE)
jpeg(filename = 'Figure_1D.jpeg', width = 10, height = 8, units = "in",res=400)
Clustered_DotPlot(seurat_object = ILC_NK, features = top_markers,
                  ggplot_default_colors = TRUE,
                  x_lab_rotate = TRUE,
                  row_label_size = 10,
                  k = 10,
                  show_parent_dend_line = TRUE
                  )
dev.off()

# Figure 1E Stacked Violin Plots Transcription Factors ----------------------------------------------------
TF <- c('TCF7','ID2','ZBTB16','RORC','AHR','IKZF2','BCL6','GATA3','MAF','EOMES','TBX21','ZNF683','IKZF3','PRDM1')
jpeg(filename = 'Figure_1E.jpeg', width = 10, height = 8, units = "in",res=400)
Stacked_VlnPlot(ILC_NK, features = TF, x_lab_rotate = TRUE, colors_use = hue_pal()(14))
dev.off()

# Figure 1F UMAP with cluster names ----------------------------------------------------
ILC_NK@active.ident <- ILC_NK$clusters
jpeg(filename = 'Figure_1F.jpeg', width = 8, height = 8, units = "in",res=400)
DimPlot(ILC_NK, reduction = "umap", label = TRUE, label.size = 5) + NoLegend()
dev.off()

# Figure 1G Violin Plot Module scoring tissue resident NK ----------------------------------------------------
trNK <- c('COTL1','CAPG','KLRC1','XCL2','ITGA1','XCL1','GSTP1','RGS1','IGFBP2','PIK3R1','KLRC2','ZFHX3',
          'ITM2A','EPAS1','TNFRSF18','PDCD4','ITGAX','B3GNT7','CD63','SRGAP3','KRT86','TRDC','DAPK2',
          'NR4A2','TNFAIP3','AFAP1L2','IL2RB','AREG','MAFF','KRT81') #Garcia-Alonso et al Nat Gen 2021
ILC_NK <- AddModuleScore(ILC_NK, features = list(trNK), name = "trNK_")

jpeg(filename = 'Figure_1G.jpeg', width = 10, height = 8, units = "in",res=400)
VlnPlot(ILC_NK, features = "trNK_1", pt.size = 0)
dev.off()

# Figure 1H Violin Plot Module scoring circulating NK ----------------------------------------------------
cNK <- c('FCGR3A','FGFBP2','CX3CR1','S1PR5','SPON2','AKR1C3','PLAC8','CD226','MYOM2','ADGRG1','SELL',
         'TGFBR3','KLF2','IL18RAP','COTL1','PLEKHG3','RASGRP2','ASCL2','RAP1GAP2','ARL4C','THEMIS2',
         'SBK1','TTC38','TBX21','CYBA','TTC16','CAPG','CD3G','KLRC1','EMP3') #Garcia-Alonso et al Nat Gen 2021
ILC_NK <- AddModuleScore(ILC_NK, features = list(cNK), name = "cNK_")
jpeg(filename = 'Figure_1H.jpeg', width = 10, height = 8, units = "in",res=400)
VlnPlot(ILC_NK, features = "cNK_1", pt.size = 0)
dev.off()

