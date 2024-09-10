#=========================================================================================================================
## DE Calling and cluster annotation
##========================================================================================================================
dir.create("results", showWarnings = F)
dir.create("results/clusters", showWarnings = F)
source("code/functions.R")

## Scale all genes from RNA assay for diff gene calling and visualisation
# Set default to RNA
DefaultAssay(patients.seurat) <- "RNA"

# Normalize
patients.seurat <- NormalizeData(patients.seurat, normalization.method = "LogNormalize")

# Scale
patients.seurat <- ScaleData(patients.seurat, features = row.names(patients.seurat))

# Correlate Ident with predicted cell type
pred.map            <- data.frame(Ident = Idents(patients.seurat), pred = patients.seurat$predicted.celltype.l2)
pred.map            <- pred.map[order(pred.map$Ident),]
d                   <- duplicated(pred.map$Ident)
pred.map            <- pred.map[!d,]
row.names(pred.map) <- NULL


##========================================================================================================================
## Define marker genes per cluster
all.patients.seurat.markers <- FindAllMarkers(object = patients.seurat, only.pos = TRUE, min.pct = 0.25, thresh = 0.25)

# Save the object
saveRDS(all.patients.seurat.markers, "Seurat_Objects/main.patients.markers.RDS")

# Save the top9 markers per cluster
sep.markers <- all.patients.seurat.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
num.clus <- length(unique(Idents(patients.seurat)))
for(i in 0:(num.clus-1)){
  # Make dot plots
  customDot(object   = patients.seurat,
            features = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
            name     = paste("results/clusters/Cluster ",i, " - ", pred.map[as.numeric(i)+1,"pred"], " - markers", sep = ""))
}

# Save markers per cluster to disk
for(i in levels(Idents(patients.seurat))){
  x <- all.patients.seurat.markers[which(all.patients.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("results/clusters/Cluster ", i, " - ", pred.map[as.numeric(i)+1,"pred"], " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

##========================================================================================================================
# Look at some genes
bunchOfCustomPlots(object    = patients.seurat,
                   features  = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                "CD79A", "ACTA2", "CD34", "KIT", "OLR1" ), 
                   ncol      = 5,
                   name      = "results/Cell type defining genes",
                   Vln.width = 25)

# NK cell markers
bunchOfCustomPlots(object   = patients.seurat,
                   features = c("CD3E", "NCAM1", "KLRD1", "FCGR3A", "IL2RB", "IL7R"), 
                   name     = "results/NK genes")

# Adipocyte markers
bunchOfCustomPlots(object   = patients.seurat,
                   features = c("LEP", "HOXC8", "HOXC9", "UCP1", "CIDEA","PRDM16"), 
                   name     = "results/Adipocyte genes")

# SMC Monocyte genes
bunchOfCustomPlots(object   = patients.seurat,
                   features = c("ACTA2", "MYH11", "CD14", "CD68", "KLF4", "CD4", "TAGLN", "ADGRE1", "LGALS3"), 
                   name     = "results/SMC mono genes")

#Extra mono cell markers
bunchOfCustomPlots(object   = patients.seurat,
                   features = c("FCGR3A", "CCR2", "CCR5", "SELL", "ITGAX", "CX3CR1"), 
                   name     = "results/Extra mono marker genes")

# T Cell markers
bunchOfCustomPlots(object   = patients.seurat,
                   features = c("TBX21","FOXP3","GATA3","CD27","CD28", "RORA", "SELL","TIGIT","STAT3","FOS","ZEB2","PRF1","GZMA","GZMK","FOSB","IL7R"), 
                   name     = "results/CD4 T genes", 
                   ncol     = 4)

bunchOfCustomPlots(object   = patients.seurat,
                   features = c("CD3E", "CD8A", "GZMA", "NKG7", "PRF1", "GNLY", "CCL5", "CX3CR1", "FGFBP2", "TBX21", "ZNF683","GZMK", "CD69", "IL7R", "SELL", "GZMB"), 
                   name     = "results/CD8 T genes", 
                   ncol     = 4)

# Myeloid Markers
bunchOfCustomPlots(object   = patients.seurat,
                   features = c("CD200R1", "TGM2", "CD163"), 
                   name     = "results/M2 genes")

bunchOfCustomPlots(object   = patients.seurat,
                   features = c("CD2","ITGAM", "PECAM1", "NCAM1", "SELL", "CSF1R"), 
                   name     = "results/Monocyte genes")

bunchOfCustomPlots(object   = patients.seurat,
                   features = c("CD68","CD14","ITGAM","CD1C", "CLEC10A", "FCER1A", "CD2", "CD3E", "CD4"), 
                   name     = "results/Myeloid genes", 
                   ncol     = 3)

bunchOfCustomPlots(object   = patients.seurat,
                   features = c("CASP1","IL1B","TLR4","ABCA1","OLR1", "TREM2","CASP4","KLF4","TNF","IL18","ABCG1", "CD9"), 
                   name     = "results/Foam vs. pro inf genes", 
                   ncol     = 6)

bunchOfCustomPlots(object   = patients.seurat,
                   features = c("CD14","CD68","CD86","CD1C", "FCGR2B", "CLEC10A", "FCER1A", "HLA-DQA2", "HLA-DQB2"), 
                   name     = "results/DC genes", 
                   ncol     = 3)


#=========================================================================================================================
## Clean up
##========================================================================================================================
rm(list = ls()[!ls() %in% c("patients.seurat", "vars.to.regress", "all.patients.seurat.markers")])
source("code/functions.R")
save.image()