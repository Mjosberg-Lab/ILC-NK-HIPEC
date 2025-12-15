#                           .....ILC NK HIPEC.....
#                       00 - Setup & Package Loading

### Package Installation ----------------------------------------------------
# NOTE: These packages need to be installed once. Uncomment and run if needed.
# See README.md for detailed installation instructions.

# # Check Seurat version
# packageVersion("Seurat")
# 
# # Install CRAN packages
# install.packages(c("Seurat", "umap", "magrittr", "dplyr", "writexl", 
#                    "devtools", "rJava", "xlsx", 
#                    "harmony", "metap", "cowplot", "reticulate", "tidyverse", 
#                    "igraph", "R.utils",
#                    "ggpmisc", "scCustomize", "fastTopics"))
# 
# # Install BiocManager and Bioconductor packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c("SingleCellExperiment", "S4Vectors", 
#                        "DESeq2", "multtest", "impute", "preprocessCore", 
#                        "GO.db", "UCell", "ComplexHeatmap",
#                        "BiocGenerics", "DelayedArray", "DelayedMatrixStats",
#                        "limma", "lme4", "SingleCellExperiment",
#                        "SummarizedExperiment", "batchelor", "Matrix.utils",
#                        "HDF5Array", "terra", "ggrastr"))
# 
# # Install GitHub packages
# library(devtools)
# remotes::install_github("mojaveazure/seurat-disk")
# remotes::install_github("satijalab/seurat-wrappers")
# devtools::install_github("cole-trapnell-lab/monocle3")
# remotes::install_github("czarnewski/niceRplots")
# devtools::install_github("kkdey/CountClust")
# devtools::install_github("KlugerLab/DAseq")
# 
# # Install archived package (specific version)
# install.packages('https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz', 
#                  repos = NULL, type = "source")

### Load Required Packages --------------------------------------------------
library(Seurat)
library(scCustomize)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(EnhancedVolcano)
library(DAseq)
library(readxl)
library(SingleCellExperiment)
library(flowCore)
library(CATALYST)
library(ggplot2)
library(cowplot)
library(writexl)