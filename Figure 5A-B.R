#                           .....ILC NK HIPEC.....
#                                  Figure 5A-B

# This script creates figures present in Figure 5A-B of the manuscript
# Marchalot et al. "Tumor-infiltrating immature innate lymphoid cells in colorectal 
# cancer are biased towards tissue-resident NK cell/ILC1 differentiation"
# Prerequisites: Load 00_SETUP.R script
# Input: ILC_NK Seurat object available in GSE302045
# Output: plots as in Figure 5A-B in Marchalot et al.


# Figure 5A Dotplots of selected genes expression for flow cytometry panel ----------------------------------------------------
jpeg(filename = 'Figure_5A.jpeg', width = 10, height = 8, units = "in",res=400)
plot <- DotPlot(ILC_NK, 
                features = c('IL7R','PTGDR2','KIT',"NCR2","HLA-DRA",'HLA-DRB1',
                             'ITGA1','ITGAE','CD7','KLRD1','KLRF1','NCAM1','FCGR3A','B3GAT1'),
                cols = c("yellow","blue3"),
                cluster.idents = FALSE
)
plot + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

# Figure 5B Dotplots of selected ADT expression for flow cytometry panel ----------------------------------------------------
DefaultAssay(ILC_NK) <- "Protein"
jpeg(filename = 'Figure_5B.jpeg', width = 10, height = 8, units = "in",res=400)
plot <- DotPlot(ILC_NK, 
                features = c('CD127-TotalSeqC','CD94-TotalSeqC','CD336-TotalSeqC',"CD45RA-TotalSeqC",
                             "CD117-TotalSeqC"),
                cols = c("yellow","blue3"),
                cluster.idents = FALSE
)
plot + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()
DefaultAssay(ILC_NK) <- "RNA"

# Figure 5C UMAP of ILC and NK cells selected on FCS files ----------------------------------------------------
