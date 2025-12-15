#                           .....ILC NK HIPEC.....
#                                  Figure 4

# This script creates figures present in Figure 4 of the manuscript
# Marchalot et al. "Tumor-infiltrating immature innate lymphoid cells in colorectal 
# cancer are biased towards tissue-resident NK cell/ILC1 differentiation"
# Prerequisites: Load 00_SETUP.R script
# Input: ILC_NK Seurat object available in GSE302045
# Output: plots as in Figure 4 in Marchalot et al.


# Figure 4A Volcano plot of nILC DE genes ----------------------------------------------------
# Differential expression analysis nILC Colon vs CRC
nILC <- subset(ILC_NK, subset = clusters == 'nILC')
DEGs <- FindMarkers(nILC, ident.1 = 'Colon', ident.2 = 'CRC', group.by = 'group', only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'wilcox')
DEGs$gene <- rownames(DEGs)
# Clean RP and MT genes
gene.list <- rownames(DEGs)  
new.gene.list <- gsub('^MT-(.*)', NA, gene.list)
new.gene.list <- gsub('^RP(.*)', NA, new.gene.list)
DEGs <- DEGs[rownames(DEGs) == new.gene.list, ]
DEGs <- na.omit(DEGs)

#Volcano plots DE
selection <- c('HSPA1A','GNLY','RNF130','HSPA1B','KLRC1','RGS2','MFGE8','DNAJB1','CARD19','PDE7B',
               'CEBPB','HLA-DRA','NKG7','IGFLR1','CTSW','HLA-DPA1','CTSC','KLRD1','APOBEC3G','ITM2C','ITGAX',
               'CSRNP1','NR4A3','IL4I1','IL1R1') # nILC Colon vs CRC selection for gene highlighting

DEGs_neg <- DEGs
DEGs_neg$avg_log2FC = DEGs_neg$avg_log2FC*(-1) #Invert FC values for figures

jpeg(filename = "Figure_4A.jpeg", width = 10, height = 10, units = "in",res=400)
EnhancedVolcano(DEGs_neg,
                lab = rownames(DEGs_neg),
                selectLab = selection,
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff= 0.01,
                FCcutoff = 0.5,
                drawConnectors = TRUE
                )
dev.off()

# Figure 4B-C Stacked violin plots of selected DE genes ----------------------------------------------------
ILC_NK@active.ident <- as.factor(ILC_NK$group)
ILC_NK_sub <- subset(x = ILC_NK, idents = c('Colon','CRC'))
ILC_NK_sub@active.ident <- ILC_NK$clusters #Subset only Colon and CRC to display

list1 <- c('GNLY','KLRC1','RGS2','CEBPB','CTSW') # All neg DEG nILC Colon vs CRC
list2 <- c('LTB','PTGER4','TIPARP','IL1R1','IL4I1') #All pos DEG nILC Colon vs CRC

jpeg(filename = "Figure_4B.jpeg", width = 12, height = 8, units = "in",res=400)
Stacked_VlnPlot(ILC_NK_sub, features = lis1, split.by = 'group', idents = 'early NK',
                x_lab_rotate = TRUE, colors_use = hue_pal()(3))
dev.off()

jpeg(filename = "Figure_4C.jpeg", width = 12, height = 8, units = "in",res=400)
Stacked_VlnPlot(ILC_NK_sub, features = list2, split.by = 'group', idents = 'nILC',
                x_lab_rotate = TRUE, colors_use = hue_pal()(3))
dev.off() 

# Figure 4D Volcano plot of eNK DE genes ----------------------------------------------------
# Differential expression analysis eNK Colon vs CRC
eNK <- subset(ILC_NK, subset = clusters == 'early NK')
DEGs <- FindMarkers(eNK, ident.1 = 'Colon', ident.2 = 'CRC', group.by = 'group', only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'wilcox')
DEGs$gene <- rownames(DEGs)
# Clean RP and MT genes
gene.list <- rownames(DEGs)  
new.gene.list <- gsub('^MT-(.*)', NA, gene.list)
new.gene.list <- gsub('^RP(.*)', NA, new.gene.list)
DEGs <- DEGs[rownames(DEGs) == new.gene.list, ]
DEGs <- na.omit(DEGs)

#Volcano plots DE
selection <- c('GNLY','SELL','HSPA1A','HSPA1B','DNAJB1','LINC00996','HAVCR2','HSPB1','BIN1',
               'LYST','IGFBP2','CRTAM','HLA-DRA','IL32','HIF1A','IL10RA','TCF7',
               'RGS1','RGS2','TNFRSF18','KLRC1','GZMK','EOMES','AHR',
               'S100A4','CD3E','ZBTB16','GZMA','LTB','IL12RB2','FOSB','KLRB1','NR4A2','MAP3K8',
               'TIPARP','AFF3','PTGER4','FOSB','ZFP36','EIF1','IER2','VPS37B',
               'H2AFX','TOX','ICAM1','NR4A3','IRF8','RORA','TEX14','SRGN','CCL4','CCL3','CCL4L2') #eNK Colon vs CRC selection for gene highlighting

DEGs_neg <- DEGs
DEGs_neg$avg_log2FC = DEGs_neg$avg_log2FC*(-1) #Invert FC values for figures

jpeg(filename = "Figure_4D.jpeg", width = 10, height = 10, units = "in",res=400)
EnhancedVolcano(DEGs_neg,
                lab = rownames(DEGs_neg),
                selectLab = selection,
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff= 0.01,
                FCcutoff = 0.5,
                drawConnectors = TRUE
                )
dev.off()

# Figure 4E-F Stacked violin plots of selected DE genes ----------------------------------------------------
list1 <- c('GNLY','CRTAM','SELL','TCF7','IL10RA','HIF1A') # All neg DEG eNK Colon vs CRC
list2 <- c('CD3E','ZBTB16','EOMES','IRF8','LTB','AHR') #All pos DEG eNK Colon vs CRC

jpeg(filename = "Figure_4E.jpeg", width = 12, height = 8, units = "in",res=400)
Stacked_VlnPlot(ILC_NK_sub, features = list1, split.by = 'group', idents = 'early NK',
                x_lab_rotate = TRUE, colors_use = hue_pal()(3))
dev.off()

jpeg(filename = "Figure_4F.jpeg", width = 12, height = 8, units = "in",res=400)
Stacked_VlnPlot(ILC_NK_sub, features = list2, split.by = 'group', idents = 'early NK',
                x_lab_rotate = TRUE, colors_use = hue_pal()(3))
dev.off() 