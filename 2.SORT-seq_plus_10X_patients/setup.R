## Install packages
install.packages(c("BiocManager", "devtools", "rmarkdown", "Rcpp", "RcppArmadillo", "xlsx"))
BiocManager::install(c("Seurat", "rDGIdb", "anndata", "Signac", "ggplot2", "patchwork", "dplyr", "reshape2", "pheatmap", "varhandle", "ggsci", "limma", "ggpubr", "AnnotationDbi", "org.Hs.eg.db", "msigdbr", "dittoSeq", "EGSEA", "biomaRt", "escape", "IRanges", "gage", "pcaMethods", "ggsignif","qusage", "GEOquery", "umap", "EnhancedVolcano"))
BiocManager::install(c("pathviewr", "drc", "AUCell", "mixtools", "decoupleR", "OmnipathR", "DESeq2","hgug4112a.db", "DropletUtils" ))
devtools::install_github("onofriAndreaPG/aomisc")

Sys.unsetenv("GITHUB_PAT")
devtools::install_github("mojaveazure/seurat-disk")
devtools::install_github("velocyto-team/velocyto.R")
remotes::install_github('satijalab/seurat-wrappers')

## Monocle 3:
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
 devtools::install_github('cole-trapnell-lab/leidenbase')
# devtools::install_github('cole-trapnell-lab/monocle3')
#
# renv::activate()
# renv::snapshot()  


## Load Packages
library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(dplyr)
library(plyr)
library(tibble)
library(tidyr)
library(reshape2)
library(pheatmap)
library(ggrepel)
library(ggsci)
library(SeuratDisk)
library(limma)
library(ggpubr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EGSEA)
library(GSEABase)
library(biomaRt)
library(velocyto.R)
library(monocle3)
library(SeuratWrappers)
library(ggsignif)
library(xlsx)
library(qusage)
library(escape)
library(msigdbr)
library(dittoSeq)
library(rstatix)
library(varhandle)
library(rDGIdb)
library(scales)
library(reticulate)
library(anndata)
library(decoupleR)
library(DESeq2)
library(BiocParallel)
library(anndata)
library(EnhancedVolcano)
library(AUCell)
library(mixtools)
library(RColorBrewer)
library(hgug4112a.db)
library(drc)
library(nlme)
library(aomisc)
library(pathviewr)
library(maptools)
register(MulticoreParam(6))

# Load custom functions
source("code/functions.R")

#Set up gene ontology objects for pathway analyses
ont_sets <- readList("cp_and_h.gmt.txt")

#Initialize bioMarts
marts <- listMarts()
mart <- marts[1,1]

human <- useMart(mart, dataset = "hsapiens_gene_ensembl")
mouse <- useMart(mart, dataset = "mmusculus_gene_ensembl")

## Order of execution
# setup.R
# load_data.R
# normalise_and_filter.R
# integrate.R
# auto_map_cell_types.R
# DE_Calling_and_cluster_annotation.R
# Myeloid_cells.R
# PBMC_monocytes.R
# plaque_macrophages.R
# full_set_43p_integratioin.R
# from_full_myeloid_43p_integration.R
# myeloid_43p_integration.R
# marker_selection.R
# final mac populations.R
# population plots.R
# population sizes.R
# metadata_correlations.R
# compare To Pan data.R
# Dib.R
# monocle.R
# velocyto.R
# decoupleR.R
# GMCSF macs.R
# Python setup.R
# cellrank.py


# myeloid_Cell_identity_back_mapping.R
