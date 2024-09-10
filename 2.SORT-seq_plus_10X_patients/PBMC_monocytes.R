#=========================================================================================================================
## Analyze the PBMC monocytes
##========================================================================================================================
dir.create("PBMC_monocytes_results",            showWarnings = F)
dir.create("PBMC_monocytes_results/clusters",   showWarnings = F)

##==========================================================================================================
## First let's subset out the pbmc only cells
# Subset the seurat object (remove the DC idents)
pbmc.mo.cells         <- WhichCells(mye.patients.seurat, expression = predicted.celltype.l1 == "Mono" & Tissue == "PBMC")
pbmc.monocytes.seurat <- subset(mye.patients.seurat, cells = pbmc.mo.cells)

# Remove trace contaminants of DCs and platelets
DefaultAssay(pbmc.monocytes.seurat) <- "RNA"
pbmc.monocytes.seurat <- subset(pbmc.monocytes.seurat, subset = CD1C == 0 & CLEC10A == 0 & GP9 == 0 & ITGB3 == 0 & GP1BB == 0 & CAVIN2 == 0)

## Recalibrate clustering
# Work with the integrated assay
DefaultAssay(pbmc.monocytes.seurat) <- "integrated"

# Run PCA and UMAP embedding
pbmc.monocytes.seurat <- RunPCA(pbmc.monocytes.seurat, verbose = T)
pbmc.monocytes.seurat <- RunUMAP(pbmc.monocytes.seurat, dims = 1:30, verbose = T)

# Show some stats
pbmc.monocytes.seurat

## Clustering
set.seed(1)
pbmc.monocytes.seurat <- FindNeighbors(pbmc.monocytes.seurat, dims = 1:30, verbose = T, force.recalc = T)
pbmc.monocytes.seurat <- FindClusters(pbmc.monocytes.seurat, verbose = T, resolution = 0.5, random.seed = 666)

# Visualise
DimPlot(pbmc.monocytes.seurat, label = TRUE, pt.size = 3, label.size = 10, shuffle = T) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("PBMC_monocytes_results/UMAP clusters.pdf", width = 15, height = 15)

DimPlot(pbmc.monocytes.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Patient", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("PBMC_monocytes_results/UMAP Patients.pdf", width = 15, height = 15)

p1 <- DimPlot(pbmc.monocytes.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Patient", shuffle = T, cells.highlight = list(P1 = WhichCells(pbmc.monocytes.seurat, expression = Patient == "P1")) ) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(pbmc.monocytes.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Patient", shuffle = T, cells.highlight = list(P2 = WhichCells(pbmc.monocytes.seurat, expression = Patient == "P2")) ) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p3 <- DimPlot(pbmc.monocytes.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Patient", shuffle = T, cells.highlight = list(P3 = WhichCells(pbmc.monocytes.seurat, expression = Patient == "P3")) ) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p1 + p2 + p3
ggsave("PBMC_monocytes_results/UMAP Patients split.pdf", width = 30, height = 10)

DimPlot(pbmc.monocytes.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Tissue", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("PBMC_monocytes_results/UMAP Tissue.pdf", width = 15, height = 15)

DimPlot(pbmc.monocytes.seurat, label = F, pt.size = 2, label.size = 10, group.by = "predicted.celltype.l1", shuffle = T, reduction = "umap") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("PBMC_monocytes_results/UMAP Predicted Cell Type L1.pdf", width = 15, height = 15)

DimPlot(pbmc.monocytes.seurat, label = F, pt.size = 2, label.size = 10, group.by = "predicted.celltype.l2", shuffle = T, reduction = "umap") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("PBMC_monocytes_results/UMAP Predicted Cell Type L2.pdf", width = 15, height = 15)


# Merge classical clusters and rename the non-classical ones while at it
pbmc.monocytes.seurat <- RenameIdents(pbmc.monocytes.seurat, "0" = "CD14+CD16-CD64+SLAN- Classical Monocytes",
                                       "1" = "CD14+CD16-CD64+SLAN- Classical Monocytes", 
                                       "2" = "CD14+CD16-CD64+SLAN- Classical Monocytes",
                                       "3" = "CD14+CD16+CD64-SLAN+ Non-Classical Monocytes",
                                       "4" = "CD14+CD16+CD64+SLAN- Intermediate Monocytes") 

num.clus        <- length(unique(Idents(pbmc.monocytes.seurat)))
DimPlot(pbmc.monocytes.seurat, label = TRUE, pt.size = 3, label.size = 10, shuffle = T, cols = c("cyan3","blue4", "blue1")) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("PBMC_monocytes_results/UMAP clusters names.pdf", width = 15, height = 15)


## Scale all genes from RNA assay for diff gene calling and visualisation
# Set default to RNA
DefaultAssay(pbmc.monocytes.seurat) <- "RNA"

# Normalize
pbmc.monocytes.seurat <- NormalizeData(pbmc.monocytes.seurat, normalization.method = "LogNormalize")

# Scale
pbmc.monocytes.seurat <- ScaleData(pbmc.monocytes.seurat, features = row.names(pbmc.monocytes.seurat))


## Save the myeloid pbmc seurat object
saveRDS(pbmc.monocytes.seurat, "Seurat_Objects/pbmc.monocytes.seurat.RDS")
#pbmc.monocytes.seurat <- readRDS(file = "Seurat_Objects/pbmc.monocytes.seurat.RDS")


##========================================================================================================================
## Define marker genes per cluster
pbmc.monocytes.seurat.markers <- FindAllMarkers(object = pbmc.monocytes.seurat, only.pos = TRUE, min.pct = 0.25, thresh = 0.25)

# Save the top markers per cluster
sep.markers <- pbmc.monocytes.seurat.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
num.clus <- length(unique(Idents(pbmc.monocytes.seurat)))
for(i in levels(Idents(pbmc.monocytes.seurat))){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = pbmc.monocytes.seurat,
                     name      = paste("PBMC_monocytes_results/clusters/Cluster ",i, " markers", sep = ""),
                     Vln.color = c("cyan3","blue4", "blue1"),
                     Vln.width = 15, Vln.height = 25)
}

# Save markers per cluster to disk
for(i in levels(Idents(pbmc.monocytes.seurat))){
  x <- pbmc.monocytes.seurat.markers[which(pbmc.monocytes.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("PBMC_monocytes_results/clusters/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

##========================================================================================================================
# Look at some genes
bunchOfCustomPlots(features   = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                  "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                  "CD79A", "ACTA2", "CD34", "KIT", "OLR1" ), 
                   object     = pbmc.monocytes.seurat,
                   ncol       = 5,
                   name       = "PBMC_monocytes_results/Cell type defining genes",
                   Vln.width  = 25, 
                   Vln.height = 25,
                   Vln.color  = c("cyan3","blue4", "blue1"))

# SMC Monocyte genes
bunchOfCustomPlots(features = c("ACTA2", "MYH11", "CD14", "CD68", "KLF4", "CD4", "TAGLN", "ADGRE1", "LGALS3"), 
                   object   = pbmc.monocytes.seurat,
                   name     = "PBMC_monocytes_results/SMC mono genes",
                   Vln.width = 15, Vln.height = 25,
                   Vln.color  = c("cyan3","blue4", "blue1"))

#Extra mono cell markers
bunchOfCustomPlots(features = c("FCGR3A", "CCR2", "CCR5", "SELL", "ITGAX", "CX3CR1"), 
                   object   = pbmc.monocytes.seurat,
                   name     = "PBMC_monocytes_results/Extra mono marker genes",
                   Vln.color  = c("cyan3","blue4", "blue1"))

# Myeloid Markers
bunchOfCustomPlots(features = c("CD200R1", "TGM2", "CD163"), 
                   object   = pbmc.monocytes.seurat,
                   name     = "PBMC_monocytes_results/M2 genes")

bunchOfCustomPlots(features = c("CD2","ITGAM", "PECAM1", "NCAM1", "SELL", "CSF1R"), 
                   object    = pbmc.monocytes.seurat,
                   name     = "PBMC_monocytes_results/Monocyte genes",
                   Vln.color  = c("cyan3","blue4", "blue1"))

bunchOfCustomPlots(features = c("CD68","CD14","ITGAM","CD1C", "CLEC10A", "FCER1A", "CD2", "CD3E", "CD4"),
                   object   = pbmc.monocytes.seurat,
                   name     = "PBMC_monocytes_results/Myeloid genes", 
                   ncol     = 3,
                   Vln.width = 15, Vln.height = 25,
                   Vln.color  = c("cyan3","blue4", "blue1"))

bunchOfCustomPlots(features = c("CASP1","IL1B","TLR4","ABCA1","OLR1", "TREM2","CASP4","KLF4","TNF","IL18","ABCG1", "CD9"), 
                   object   = pbmc.monocytes.seurat,
                   name     = "PBMC_monocytes_results/Foam vs. pro inf genes", 
                   ncol     = 6,
                   Vln.width = 30, Vln.height = 25,
                   Vln.color  = c("cyan3","blue4", "blue1"))

bunchOfCustomPlots(features = c("CD14","CD68","CD86","CD1C", "FCGR2B", "CLEC10A", "FCER1A", "HLA-DQA2", "HLA-DQB2"), 
                   object   = pbmc.monocytes.seurat,
                   name     = "PBMC_monocytes_results/DC genes", 
                   ncol     = 3,
                   Vln.width = 15, Vln.height = 25,
                   Vln.color  = c("cyan3","blue4", "blue1"))

bunchOfCustomPlots(features = c("NLRP3", "IL1B", "CASP1", "CASP4", "IL18", "PYCARD"), 
                   object   = pbmc.monocytes.seurat,
                   name     = "PBMC_monocytes_results/Inflammasome genes", 
                   ncol     = 3,
                   Vln.color  = c("cyan3","blue4", "blue1"))

bunchOfCustomPlots(features = c("CD14", "FCGR3A", "FCGR1A", "CCR2", "PADI4", "SECISBP2L", "CD86", "ITGAX", "LYZ", "CD33", "CD36", "ITGAM"), 
                   object   = pbmc.monocytes.seurat,
                   name     = "PBMC_monocytes_results/Mono subtype genes CD14 CD16 CD64 CCR2 PADI4 SLAN CD86 CD11c LYZ CD33 CD36 CD11b", 
                   ncol     = 3,
                   Vln.width = 15, Vln.height = 25,
                   Vln.color  = c("cyan3","blue4", "blue1"))

bunchOfCustomPlots(features = c("NR4A1"), 
                   object   = pbmc.monocytes.seurat,
                   name     = "PBMC_monocytes_results/NR4A1",
                   ncol     = 1,
                   Vln.color  = c("cyan3","blue4", "blue1"))

bunchOfCustomPlots(features  = c("CD14","FCGR3A", "CD68", 'CD86', "SELL", "CSF1R", "ITGAM", "ITGAX", "PECAM1"), 
                   object    = pbmc.monocytes.seurat,
                   ncol      = 3, Vln.stack = T,
                   name      = "PBMC_monocytes_results/macro and mono genes",
                   Vln.width = 10, Vln.height = 25,
                   Vln.color  = c("cyan3","blue4", "blue1"))

##==================================================================================
##==================================================================================
## Define unique and significant markers
uniq.pbmc.monocytes.seurat.markers <- pbmc.monocytes.seurat.markers[!(duplicated(pbmc.monocytes.seurat.markers$gene) | duplicated(pbmc.monocytes.seurat.markers$gene, fromLast = T)),]

# Save unique markers per cluster to disk
for(i in levels(Idents(pbmc.monocytes.seurat))){
  x <- uniq.pbmc.monocytes.seurat.markers[which(uniq.pbmc.monocytes.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("PBMC_monocytes_results/clusters/Cluster ", i, " unique var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

# Save significant unique markers per cluster to disk
for(i in levels(Idents(pbmc.monocytes.seurat))){
  x <- uniq.pbmc.monocytes.seurat.markers[which(uniq.pbmc.monocytes.seurat.markers$cluster == i & uniq.pbmc.monocytes.seurat.markers$p_val_adj < 0.1),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("PBMC_monocytes_results/clusters/Cluster ", i, " unique and significant var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

# Save significant markers per cluster to disk
for(i in levels(Idents(pbmc.monocytes.seurat))){
  x <- pbmc.monocytes.seurat.markers[which(pbmc.monocytes.seurat.markers$cluster == i & pbmc.monocytes.seurat.markers$p_val_adj < 0.1),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("PBMC_monocytes_results/clusters/Cluster ", i, " significant var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}


# plot the top9 unique markers per cluster
for(i in unique(Idents(pbmc.monocytes.seurat))){
  # Make feature, violin, an dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(uniq.pbmc.monocytes.seurat.markers[uniq.pbmc.monocytes.seurat.markers$cluster == i, "gene"]))[1:9], 
                     object    = pbmc.monocytes.seurat,
                     Vln.stack = T,
                     name      = paste("PBMC_monocytes_results/clusters/Cluster ",i, " unique markers", sep = ""),
                     Vln.color  = c("cyan3","blue4", "blue1"))
}


##==================================================================================
##==================================================================================
## Get ontologies
# All markers
for(i in levels(Idents(pbmc.monocytes.seurat))){
  x <- pbmc.monocytes.seurat.markers[which(pbmc.monocytes.seurat.markers$cluster == i),-6]
  get_ontology(x, name = i, outdir = "PBMC_monocytes_results/clusters", universe = pbmc.monocytes.seurat, full_GSEA = F, volcano.plot = T) 
}

# Unique markers
for(i in levels(Idents(pbmc.monocytes.seurat))){
  x <- uniq.pbmc.monocytes.seurat.markers[which(uniq.pbmc.monocytes.seurat.markers$cluster == i),-6]
  get_ontology(x, name = paste(i, ".unique_markers", sep = ""), outdir = "PBMC_monocytes_results/clusters", universe = pbmc.monocytes.seurat, full_GSEA = F, volcano.plot = T) 
}

