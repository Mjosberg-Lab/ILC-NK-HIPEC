#                           .....ILC NK HIPEC.....
#                                  Figure 3A-D

# This script creates figures present in Figure 3A-B of the manuscript and prepares data for figure 3D
# Marchalot et al. "Tumor-infiltrating immature innate lymphoid cells in colorectal 
# cancer are biased towards tissue-resident NK cell/ILC1 differentiation"
# Prerequisites: Load ILC_NK Seurat object available in GSE302045
# Input: Seurat object
# Output: plots as in Figure 3A-B in Marchalot et al.


# Figure 3A Barplots tissue composition ----------------------------------------------------
# Downsampling to the smallest tissue size (n=3185 cells in Colon)
x <- ILC_NK
length(x$group[x$group %in% 'Colon']) 
length(x$group[x$group %in% 'CRC'])  
length(x$group[x$group %in% 'PC'])  
set.seed(6)
Colon <- subset(x, subset = group == 'Colon')
CRC <- subset(x, subset = group == 'CRC')
CRC <- CRC[, sample(colnames(CRC), size = length(x$group[x$group %in% 'Colon']), replace = F)]
PC <- subset(x, subset = group == 'PC')
PC <- PC[, sample(colnames(PC), size = length(x$group[x$group %in% 'Colon']), replace = F)]
data.downsG <- merge(Colon, list(CRC,PC))

# Composition of clusters by group - Downsampled
Counts <- table(Idents(data.downsG), data.downsG$group)
Counts <- as.data.frame(Counts)
plot_file1 <- paste("Figure_3A.jpeg",sep = "")
jpeg(filename = plot_file1, width = 12, height = 10, units = "in",res=400)
ggplot(data=Counts, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar(stat = "Identity", position = "fill") +
  labs(x = "Clusters", y = "Frequency", fill = "Tissue") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1.2),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))
dev.off()

# Figure 3B UMAP of cells by tissues ----------------------------------------------------
jpeg(filename = "Figure_3B.jpeg", width = 10, height = 8, units = "in",res=400)
DimPlot(ILC_NK, reduction = "umap", label.size = 5, group.by = "group")
dev.off()

# Figure 3D Export Seurat objects split by tissues to H5ad for Velocity analysis in Python ----------------------------------------------------
Colon <- subset(x, subset = group == 'Colon')
SaveH5Seurat(Colon, filename = "/path/to/seurat_object_Colon", overwrite = T)
Convert(source =  paste(path,'.h5seurat',sep = ""), assay ="RNA", dest = "h5ad", overwrite = T)
CRC <- subset(x, subset = group == 'CRC')
SaveH5Seurat(CRC, filename = "/path/to/seurat_object_CRC", overwrite = T)
Convert(source =  paste(path,'.h5seurat',sep = ""), assay ="RNA", dest = "h5ad", overwrite = T)
PC <- subset(x, subset = group == 'PC')
SaveH5Seurat(PC, filename = "/path/to/seurat_object_PC", overwrite = T)
Convert(source =  paste(path,'.h5seurat',sep = ""), assay ="RNA", dest = "h5ad", overwrite = T)

