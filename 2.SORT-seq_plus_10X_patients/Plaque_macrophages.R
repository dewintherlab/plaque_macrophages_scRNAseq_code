#=========================================================================================================================
## Analyze the Plaque macrophages
##========================================================================================================================
dir.create("Plaque_macrophages_results",            showWarnings = F)
dir.create("Plaque_macrophages_results/clusters",   showWarnings = F)

##==========================================================================================================
## First let's subset out the plaque only cells
# Subset the seurat object (remove the DC idents)
plaque.ma.cells           <- WhichCells(mye.patients.seurat, expression = predicted.celltype.l1 == "Mono" & predicted.celltype.l2 != "cDC2"  & Tissue == "plaque")
plaque.macrophages.seurat <- subset(mye.patients.seurat, cells = plaque.ma.cells)

# Remove trace contaminants
DefaultAssay(plaque.macrophages.seurat) <- "RNA"
plaque.macrophages.seurat <- subset(plaque.macrophages.seurat, subset = CD1C == 0 & CLEC10A == 0)

## Recalibrate clustering
# Work with the integrated assay
DefaultAssay(plaque.macrophages.seurat) <- "integrated"

# Run PCA and UMAP embedding
plaque.macrophages.seurat <- RunPCA(plaque.macrophages.seurat, verbose = T)
plaque.macrophages.seurat <- RunUMAP(plaque.macrophages.seurat, dims = 1:30, verbose = T)

# Show some stats
plaque.macrophages.seurat

## Clustering
set.seed(1)
plaque.macrophages.seurat <- FindNeighbors(plaque.macrophages.seurat, dims = 1:30, verbose = T, force.recalc = T)
plaque.macrophages.seurat <- FindClusters(plaque.macrophages.seurat, verbose = T, resolution = 0.6, random.seed = 666)

# Visualise
DimPlot(plaque.macrophages.seurat, label = TRUE, pt.size = 3, label.size = 10, shuffle = T) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("Plaque_macrophages_results/UMAP clusters.pdf", width = 15, height = 15)

# Refine manually
plot         <- DimPlot(plaque.macrophages.seurat, reduction = "umap")
select.cells <- CellSelector(plot = plot)
Idents(plaque.macrophages.seurat, cells = select.cells) <- "6"

select.cells <- CellSelector(plot = plot)
Idents(plaque.macrophages.seurat, cells = select.cells) <- "1.1"

select.cells <- CellSelector(plot = plot)
Idents(plaque.macrophages.seurat, cells = select.cells) <- "1.1"

select.cells <- CellSelector(plot = plot)
Idents(plaque.macrophages.seurat, cells = select.cells) <- "2.1"

# Visualise
DimPlot(plaque.macrophages.seurat, label = TRUE, pt.size = 3, label.size = 10, shuffle = T) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("Plaque_macrophages_results/UMAP clusters refined.pdf", width = 15, height = 15)

DimPlot(plaque.macrophages.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Patient", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("Plaque_macrophages_results/UMAP Patients.pdf", width = 15, height = 15)

p1 <- DimPlot(plaque.macrophages.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Patient", shuffle = T, cells.highlight = list(P1 = WhichCells(plaque.macrophages.seurat, expression = Patient == "P1")) ) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(plaque.macrophages.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Patient", shuffle = T, cells.highlight = list(P2 = WhichCells(plaque.macrophages.seurat, expression = Patient == "P2")) ) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p3 <- DimPlot(plaque.macrophages.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Patient", shuffle = T, cells.highlight = list(P3 = WhichCells(plaque.macrophages.seurat, expression = Patient == "P3")) ) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p1+p2+p3
ggsave("Plaque_macrophages_results/UMAP Patients split.pdf", width = 30, height = 10)

DimPlot(plaque.macrophages.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Tissue", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("Plaque_macrophages_results/UMAP Tissue.pdf", width = 15, height = 15)

DimPlot(plaque.macrophages.seurat, label = F, pt.size = 2, label.size = 10, group.by = "predicted.celltype.l1", shuffle = T, reduction = "umap") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("Plaque_macrophages_results/UMAP Predicted Cell Type L1.pdf", width = 15, height = 15)

DimPlot(plaque.macrophages.seurat, label = F, pt.size = 2, label.size = 10, group.by = "predicted.celltype.l2", shuffle = T, reduction = "umap") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("Plaque_macrophages_results/UMAP Predicted Cell Type L2.pdf", width = 15, height = 15)


## Scale all genes from RNA assay for diff gene calling and visualisation
# Set default to RNA
DefaultAssay(plaque.macrophages.seurat) <- "RNA"

# Normalize
plaque.macrophages.seurat <- NormalizeData(plaque.macrophages.seurat, normalization.method = "LogNormalize")

# Scale
plaque.macrophages.seurat <- ScaleData(plaque.macrophages.seurat, features = row.names(plaque.macrophages.seurat))


## Save the myeloid pbmc seurat object
saveRDS(plaque.macrophages.seurat, "Seurat_Objects/plaque.macrophages.seurat.RDS")
#plaque.macrophages.seurat <- readRDS(file = "Seurat_Objects/plaque.macrophages.seurat.RDS")


##========================================================================================================================
## Define marker genes per cluster
plaque.macrophages.seurat.markers <- FindAllMarkers(object = plaque.macrophages.seurat, only.pos = TRUE, min.pct = 0.25, thresh = 0.25)

# Save the top9 markers per cluster
sep.markers <- plaque.macrophages.seurat.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
num.clus <- length(unique(Idents(plaque.macrophages.seurat)))
for(i in levels(Idents(plaque.macrophages.seurat))){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = plaque.macrophages.seurat,
                     name      = paste("Plaque_macrophages_results/clusters/Cluster ",i, " markers", sep = ""),
                     
                     Vln.width = 25, Vln.height = 25)
}

# Save markers per cluster to disk
for(i in levels(Idents(plaque.macrophages.seurat))){
  x <- plaque.macrophages.seurat.markers[which(plaque.macrophages.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("Plaque_macrophages_results/clusters/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

##========================================================================================================================
# Look at some genes
bunchOfCustomPlots(features   = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                  "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                  "CD79A", "ACTA2", "CD34", "KIT", "OLR1" ), 
                   object     = plaque.macrophages.seurat,
                   ncol       = 5,
                   name       = "Plaque_macrophages_results/Cell type defining genes",
                   Vln.width  = 25, 
                   Vln.height = 15
                   )

# SMC Monocyte genes
bunchOfCustomPlots(features = c("ACTA2", "MYH11", "CD14", "CD68", "KLF4", "CD4", "TAGLN", "ADGRE1", "LGALS3"), 
                   object   = plaque.macrophages.seurat,
                   name     = "Plaque_macrophages_results/SMC mono genes",
                   Vln.width = 25, Vln.height = 15
                   )

#Extra mono cell markers
bunchOfCustomPlots(features = c("FCGR3A", "CCR2", "CCR5", "SELL", "ITGAX", "CX3CR1"), 
                   object   = plaque.macrophages.seurat,
                   name     = "Plaque_macrophages_results/Extra mono marker genes",
                   )

# Myeloid Markers
bunchOfCustomPlots(features = c("CD200R1", "TGM2", "CD163"), 
                   object   = plaque.macrophages.seurat,
                   name     = "Plaque_macrophages_results/M2 genes")

bunchOfCustomPlots(features = c("CD2","ITGAM", "PECAM1", "NCAM1", "SELL", "CSF1R"), 
                   object    = plaque.macrophages.seurat,
                   name     = "Plaque_macrophages_results/Monocyte genes",
                   )

bunchOfCustomPlots(features = c("CD68","CD14","ITGAM","CD1C", "CLEC10A", "FCER1A", "CD2", "CD3E", "CD4"),
                   object   = plaque.macrophages.seurat,
                   name     = "Plaque_macrophages_results/Myeloid genes", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
                   )

bunchOfCustomPlots(features = c("CASP1","IL1B","TLR4","ABCA1","OLR1", "TREM2","CASP4","KLF4","TNF","IL18","ABCG1", "CD9"), 
                   object   = plaque.macrophages.seurat,
                   name     = "Plaque_macrophages_results/Foam vs. pro inf genes", 
                   ncol     = 3,
                   Vln.width = 30, Vln.height = 15
                   )

bunchOfCustomPlots(features = c("CD14","CD68","CD86","CD1C", "FCGR2B", "CLEC10A", "FCER1A", "HLA-DQA2", "HLA-DQB2"), 
                   object   = plaque.macrophages.seurat,
                   name     = "Plaque_macrophages_results/DC genes", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
                   )

bunchOfCustomPlots(features = c("NLRP3", "IL1B", "CASP1", "CASP4", "IL18", "PYCARD"), 
                   object   = plaque.macrophages.seurat,
                   name     = "Plaque_macrophages_results/Inflammasome genes", 
                   ncol     = 3,
                   )

bunchOfCustomPlots(features = c("CD14", "FCGR3A", "FCGR1A", "CCR2", "PADI4", "SECISBP2L", "CD86", "ITGAX", "LYZ", "CD33", "CD36", "ITGAM"), 
                   object   = plaque.macrophages.seurat,
                   name     = "Plaque_macrophages_results/Mono subtype genes CD14 CD16 CD64 CCR2 PADI4 SLAN CD86 CD11c LYZ CD33 CD36 CD11b", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
                   )

bunchOfCustomPlots(features = c("NR4A1"), 
                   object   = plaque.macrophages.seurat,
                   name     = "Plaque_macrophages_results/NR4A1",
                   ncol     = 1,
                   )

bunchOfCustomPlots(features  = c("CD14","FCGR3A", "CD68", 'CD86', "SELL", "CSF1R", "ITGAM", "ITGAX", "PECAM1"), 
                   object    = plaque.macrophages.seurat,
                   ncol      = 3, Vln.stack = T,
                   name      = "Plaque_macrophages_results/macro and mono genes",
                   Vln.width = 20, Vln.height = 20
                   )

bunchOfCustomPlots(features  = c("TREM2", "CD9", "GPNMB", "OLR1", "ABCA1", "ABCG1", "IL1B", "TNF", "CASP1", "CD14", "FCGR1A", "FCGR3A"), 
                   object    = plaque.macrophages.seurat,
                   ncol      = 3,
                   name      = "Plaque_macrophages_results/resident - foam - inflammatory - general markers",
                   Vln.width = 20, Vln.height = 20
)


bunchOfCustomPlots(features = c("TIMD4", "CCR2", "LYVE1", "FOLR2"), 
                   object   = mye.patients.seurat, assay = "RNA",
                   name     = "Plaque_macrophages_results/core_Gene_macs",
                   ncol     = 4, 
                   Vln.width = 25
)


customFeature(object = plaque.macrophages.seurat, features = c("MKI67"), name = "Plaque_macrophages_results/MKI67.pdf")

DefaultAssay(mye.patients.seurat) <- "integrated"
FeaturePlot(mye.patients.seurat, features = c("TIMD4", "LYVE1"))


##==================================================================================
##==================================================================================
## Define unique and significant markers
uniq.plaque.macrophages.seurat.markers <- plaque.macrophages.seurat.markers[!(duplicated(plaque.macrophages.seurat.markers$gene) | duplicated(plaque.macrophages.seurat.markers$gene, fromLast = T)),]

# Save unique markers per cluster to disk
for(i in levels(Idents(plaque.macrophages.seurat))){
  x <- uniq.plaque.macrophages.seurat.markers[which(uniq.plaque.macrophages.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("Plaque_macrophages_results/clusters/Cluster ", i, " unique var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

# Save significant unique markers per cluster to disk
for(i in levels(Idents(plaque.macrophages.seurat))){
  x <- uniq.plaque.macrophages.seurat.markers[which(uniq.plaque.macrophages.seurat.markers$cluster == i & uniq.plaque.macrophages.seurat.markers$p_val_adj < 0.1),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("Plaque_macrophages_results/clusters/Cluster ", i, " unique and significant var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

# Save significant markers per cluster to disk
for(i in levels(Idents(plaque.macrophages.seurat))){
  x <- plaque.macrophages.seurat.markers[which(plaque.macrophages.seurat.markers$cluster == i & plaque.macrophages.seurat.markers$p_val_adj < 0.1),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("Plaque_macrophages_results/clusters/Cluster ", i, " significant var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}


# plot the top9 unique markers per cluster
for(i in unique(Idents(plaque.macrophages.seurat))){
  # Make feature, violin, an dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(uniq.plaque.macrophages.seurat.markers[uniq.plaque.macrophages.seurat.markers$cluster == i, "gene"]))[1:9], 
                     object    = plaque.macrophages.seurat,
                     Vln.stack = T,
                     name      = paste("Plaque_macrophages_results/clusters/Cluster ",i, " unique markers", sep = ""),
                     Vln.color  = c("cyan3","blue4", "blue1"))
}


##==================================================================================
##==================================================================================
## Get ontologies
# All markers
for(i in levels(Idents(plaque.macrophages.seurat))){
  x <- plaque.macrophages.seurat.markers[which(plaque.macrophages.seurat.markers$cluster == i),-6]
  get_ontology(x, name = i, outdir = "Plaque_macrophages_results/clusters", universe = plaque.macrophages.seurat, full_GSEA = F, volcano.plot = T) 
}

# Top 8
for(i in levels(Idents(plaque.macrophages.seurat))){
  x <- plaque.macrophages.seurat.markers[which(plaque.macrophages.seurat.markers$cluster == i),-6]
  get_ontology(x, name = paste("Top 8", i, sep = " "), outdir = "Plaque_macrophages_results/clusters", universe = plaque.macrophages.seurat, full_GSEA = F, volcano.plot = F, plot.top.n = 8) 
}

# Unique markers
for(i in levels(Idents(plaque.macrophages.seurat))){
  x <- uniq.plaque.macrophages.seurat.markers[which(uniq.plaque.macrophages.seurat.markers$cluster == i),-6]
  get_ontology(x, name = paste(i, ".unique_markers", sep = ""), outdir = "Plaque_macrophages_results/clusters", universe = plaque.macrophages.seurat, full_GSEA = F, volcano.plot = T) 
}

##==================================================================================
##==================================================================================
# Rename the populations
plaque.macrophages.seurat <- StashIdent(plaque.macrophages.seurat, save.name = "numbered.idents")
plaque.macrophages.seurat <- RenameIdents(plaque.macrophages.seurat, 
                                         "0"   = "CD14+IL1B+SELL+S100A8+ Migrating Antigen-presenting Inflammatory Monocyte-derived Macrophages",
                                         "3"   = "CD14+IL1B++SELL+OLR1+ Migrating Highly Inflammatory Foamy Monocyte-derived Macrophages",
                                         "6"   = "CD14+IL1B+SELL+ABCA+  Migrating Antigen-presenting Inflammatory Monocyte-derived Macrophages",
                                         "1.1" = "CD14+IL1B+FN1+OLR1+ Inflammatory Foamy Macrophages",
                                         "4"   = "CD14+-IL1B+SELL+CD16++ Migrating Antigen-presenting Inflammatory Monocyte-derived Macrophages",
                                         "2"   = "CD14+IL1B-TREM2+ABCA1- Anti-inflammatory Resident-like Macrophages",
                                         "2.1" = "CD14+IL1B-TREM2+ABCA1+ Anti-inflammatory Foamy Resident-like Macrophages",
                                         "1"   = "CD14+IL1B+TREM2+ABCA1+ Inflammatory Foamy Resident-like Macrophages",
                                         "5"   = "CD14+IL1B+TREM2-HSPA6+ Lipid-stress Activated Foamy Macrophages"
                                         )

DimPlot(plaque.macrophages.seurat, label = TRUE, pt.size = 3, label.size = 8, shuffle = T, cols = M_refined.pop.colors[4:12]) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("Plaque_macrophages_results/UMAP clusters names.pdf", width = 15, height = 15)

# Plot with a legend at the bottom, use colors defined later in map_back
DimPlot(plaque.macrophages.seurat, label = F, pt.size = 3, shuffle = T, cols = M_refined.pop.colors[4:12]) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "bottom") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold"),
        legend.text      = element_text(size = 18, face = "bold")
  )
ggsave("Plaque_macrophages_results/UMAP clusters names legend.pdf", width = 15, height = 15)


## Save the seurat object
saveRDS(plaque.macrophages.seurat, "Seurat_Objects/plaque.macrophages.seurat.RDS")
#plaque.macrophages.seurat <- readRDS(file = "Seurat_Objects/plaque.macrophages.seurat.RDS")


## Re-plot some markers
customVln(features  = c("TREM2", "CD9", "GPNMB", "OLR1", "ABCA1", "ABCG1", "IL1B", "TNF", "CASP1", "CD14", "FCGR1A", "FCGR3A"), 
                   object    = plaque.macrophages.seurat,
                   ncol      = 3,
                   name      = "Plaque_macrophages_results/resident - foam - inflammatory - general markers - names.pdf", 
                   cols      = c("#B7701C" , "#A1451B","#8B1A1A", "#A37CC4","#CD9B1D", "#7FFF00", "#61C400", "#458B00","#483D8B"),
                   width = 20, height = 50
)


customVln(features = c("CASP1","IL1B","TLR4","ABCA1","OLR1", "TREM2","CASP4","KLF4","TNF","IL18","ABCG1", "CD9"), 
                   object = plaque.macrophages.seurat,
                   name   = "Plaque_macrophages_results/Foam vs. pro inf genes - names.pdf", 
                   ncol   = 3,
                   cols   = c("#B7701C" , "#A1451B","#8B1A1A", "#A37CC4","#CD9B1D", "#7FFF00", "#61C400", "#458B00","#483D8B"),
                   width  = 30, height = 50
)

customVln(features = c("IL1B", "FCGR3A"), 
          object = plaque.macrophages.seurat,
          name   = "Plaque_macrophages_results/zoom in IL1B FCGR3A.pdf", 
          ncol   = 2,
          width  = 20, height = 10
)

bunchOfCustomPlots(object    = plaque.macrophages.seurat, 
                   features  = c("TREM2", "CD9", "GPNMB", "OLR1", "ABCA1", "ABCG1", "IL1B", "TNF", "CASP1", "CD14", "SELL", "FCGR3A"), 
                   name      = "Plaque_macrophages_results/Marker roundup", 
                   Vln.stack = T, Vln.height = 20, Dot.width = 20, dot.scale = 10,
                   ncol      = 3)

customDot(object    = plaque.macrophages.seurat, 
          features  = c("TREM2", "CD9", "GPNMB", "OLR1", "ABCA1", "ABCG1", "IL1B", "TNF", "CASP1", "CD14", "SELL", "FCGR3A"), 
          name      = "Plaque_macrophages_results/Marker roundup - dot plot.pdf", dot.scale = 10, width = 20, cluster.idents = T)

customFeature(object    = plaque.macrophages.seurat, 
          features  = c("MTOR"), 
          name      = "Plaque_macrophages_results/MTOR - feature plot.pdf")


## IL1B production vs. foaminess
# Look at IL1B, LXR, and PPARG pathways
bunchOfCustomPlots(object    = plaque.macrophages.seurat, 
                   features  = c("IL1B", "NR1H3", "PPARG", "RELA", "TNF", "CCL2", "ABCA1", "FASN", "RXRA","APOE" ,"CREBBP", "EP300", "CEBPB", "SREBF1", "STAT3", "EGR2", "FRA2", "RARA"), 
                   name      = "Plaque_macrophages_results/IL1B LXR PPARG roundup", 
                   Vln.stack = T, Vln.height = 20, Dot.width = 20, dot.scale = 10,
                   ncol      = 3,
                   cols      = M_refined.pop.colors)

# Look at specific Faomy markers
bunchOfCustomPlots(object   = plaque.macrophages.seurat,
          features   = c("IL1B", "NR1H3", "PPARG", "ABCA1", "OLR1", "SREBF1"), 
          ncol       = 3, 
          name       = "Plaque_macrophages_results/IL1B vs. foam", 
          Vln.height = 30, 
          Vln.width  = 20, 
          Vln.color  = M_refined.pop.colors)

FeaturePlot(object = plaque.macrophages.seurat, features = c("IL1B", "OLR1"), pt.size = 3, blend = T)
ggsave(filename = "Plaque_macrophages_results/IL1B vs. OLR1.pdf", width = 30, height = 10)
FeaturePlot(object = plaque.macrophages.seurat, features = c("IL1B", "ABCA1"), pt.size = 3, blend = T)
ggsave(filename = "Plaque_macrophages_results/IL1B vs. ABCA1.pdf", width = 30, height = 10)
FeaturePlot(object = plaque.macrophages.seurat, features = c("IL1B", "SREBF1"), pt.size = 3, blend = T)
ggsave(filename = "Plaque_macrophages_results/IL1B vs. SREBF1.pdf", width = 30, height = 10)

##======================================================================================================
## Separate out IL1B positive cells
IL1B <- GetAssayData(plaque.macrophages.seurat, assay = "RNA")["IL1B",]
summary(IL1B)
IL1B <- names(IL1B)[IL1B > 0]

plaque.macrophages.seurat <- AddMetaData(plaque.macrophages.seurat, metadata = rep("FALSE", ncol(plaque.macrophages.seurat)), col.name = "IL1B_expr")
plaque.macrophages.seurat$IL1B_expr[names(plaque.macrophages.seurat$IL1B_expr) %in% IL1B] <- "TRUE"

# IL1B pos cells sonly
bunchOfCustomPlots(object     = subset(plaque.macrophages.seurat, cells = IL1B),
                   features   = c("IL1B", "NR1H3", "PPARG", "ABCA1", "OLR1", "SREBF1"), 
                   ncol       = 3, 
                   name       = "Plaque_macrophages_results/IL1B Pos cells vs.foam", 
                   Vln.height = 30, 
                   Vln.width  = 20, 
                   Vln.color  = M_refined.pop.colors)

# IL1B pos neg cells split violin
P1 <- VlnPlot(object     = plaque.macrophages.seurat,
              features   = c("NR1H3", "PPARG", "ABCA1", "ABCG1", "OLR1", "SREBF1"),
              split.plot = T, 
              split.by   = "IL1B", 
              ncol       = 3, 
              cols       = c("grey", "firebrick"))

p1[[1]] <- p1[[1]] + theme(axis.text.x = element_blank())
p1[[2]] <- p1[[2]] + theme(axis.text.x = element_blank())
p1[[3]] <- p1[[3]] + theme(axis.text.x = element_blank())
p1[[1]] <- p1[[1]] + theme(legend.position = "none")
p1[[2]] <- p1[[2]] + theme(legend.position = "none")
p1[[3]] <- p1[[3]] + theme(legend.position = "none")
p1[[4]] <- p1[[4]] + theme(legend.position = "none")
p1[[5]] <- p1[[5]] + theme(legend.position = "none")

wrap_elements(grid::textGrob(' ')) | p1
ggsave("Plaque_macrophages_results/IL1B vs. foam split violin.pdf", width = 15, height = 20)

IL1B.freq <- vector()
for (theIdent in unique(Idents(plaque.macrophages.seurat))){
  IL1B.freq[[theIdent]] <- ncol(subset(plaque.macrophages.seurat, ident = theIdent, IL1B == TRUE))
}

total.freq <- vector()
for (theIdent in unique(Idents(plaque.macrophages.seurat))){
  total.freq[[theIdent]] <- ncol(subset(plaque.macrophages.seurat, ident = theIdent))
}
freq.ratio <- IL1B.freq / total.freq

##================================================================================================================================
## Stratify by expression of some key genes
dir.create("Plaque_macrophages_results/Stratification", recursive = T, showWarnings = F)

plaque.macrophages.seurat <- stratifyByExpression(object        = plaque.macrophages.seurat, 
                                                  strat.by      = "IL1B", 
                                                  file.name     = "Plaque_macrophages_results/Stratification/IL1B stratified", 
                                                  verbose       = F)

plaque.macrophages.seurat <- stratifyByExpression(object        = plaque.macrophages.seurat, 
                                                  strat.by      = "OLR1", 
                                                  file.name     = "Plaque_macrophages_results/Stratification/OLR1 stratified", 
                                                  verbose       = F)

plaque.macrophages.seurat <- stratifyByExpression(object        = plaque.macrophages.seurat, 
                                                  strat.by      = "ABCA1", 
                                                  file.name     = "Plaque_macrophages_results/Stratification/ABCA1 stratified", 
                                                  verbose       = F)

plaque.macrophages.seurat <- stratifyByExpression(object        = plaque.macrophages.seurat, 
                                                  strat.by      = "TREM2", 
                                                  file.name     = "Plaque_macrophages_results/Stratification/TREM2 stratified", 
                                                  verbose       = F)

plaque.macrophages.seurat <- stratifyByExpression(object        = plaque.macrophages.seurat, 
                                                  strat.by      = "CD9", 
                                                  file.name     = "Plaque_macrophages_results/Stratification/CD9 stratified", 
                                                  verbose       = F)

plaque.macrophages.seurat <- stratifyByExpression(object        = plaque.macrophages.seurat, 
                                                  strat.by      = "FOLR2", 
                                                  file.name     = "Plaque_macrophages_results/Stratification/FOLR2 stratified", 
                                                  verbose       = F)

plaque.macrophages.seurat <- stratifyByExpression(object        = plaque.macrophages.seurat, 
                                                  strat.by      = "MRC1", 
                                                  file.name     = "Plaque_macrophages_results/Stratification/MRC1 stratified", 
                                                  verbose       = F)

plaque.macrophages.seurat <- stratifyByExpression(object        = plaque.macrophages.seurat, 
                                                  strat.by      = "LYVE1", 
                                                  file.name     = "Plaque_macrophages_results/Stratification/LYVE1 stratified", 
                                                  verbose       = F)

plaque.macrophages.seurat <- stratifyByExpression(object        = plaque.macrophages.seurat, 
                                                  strat.by      = "GPNMB", 
                                                  file.name     = "Plaque_macrophages_results/Stratification/GPNMB stratified", 
                                                  verbose       = F)

plaque.macrophages.seurat <- stratifyByExpression(object        = plaque.macrophages.seurat, 
                                                  strat.by      = "SPP1", 
                                                  file.name     = "Plaque_macrophages_results/Stratification/SPP1 stratified", 
                                                  verbose       = F)

plaque.macrophages.seurat <- stratifyByExpression(object        = plaque.macrophages.seurat, 
                                                  strat.by      = "TNFAIP3", 
                                                  file.name     = "Plaque_macrophages_results/Stratification/TNFAIP3 stratified", 
                                                  verbose       = F)

plaque.macrophages.seurat <- stratifyByExpression(object        = plaque.macrophages.seurat, 
                                                  strat.by      = "TNF", 
                                                  file.name     = "Plaque_macrophages_results/Stratification/TNF stratified", 
                                                  verbose       = F)


##=============================================================================================================
## Get markers for IL1B stratification
# High vs. Zero
IL1B.strat.markers <- FindMarkers(object = plaque.macrophages.seurat, min.pct = 0.25, logfc.threshold = 0.25, ident.1 = "High", ident.2 = "Zero", group.by = "IL1B_expr", only.pos = T)

bunchOfCustomPlots(object    = plaque.macrophages.seurat, 
                   features  = row.names(IL1B.strat.markers)[1:9], 
                   group.by  = "IL1B_expr", 
                   name      = "Plaque_macrophages_results/IL1B expression level stratification markers", 
                   pt.size   = 0, 
                   dot.scale = 10,
                   Vln.color = c("grey", "bisque", "coral", "firebrick"), 
                   ncol      = 3)

x <- subset(IL1B.strat.markers, p_val_adj < 0.05)
x$gene <- row.names(x)
nrow(IL1B.strat.markers)
nrow(x)
get_ontology(x, name = "IL1B expression level stratification padj 0.05", outdir = "Plaque_macrophages_results/", full_GSEA = F, universe = plaque.macrophages.seurat, plot.top.n = 10)

x <- IL1B.strat.markers
x$gene <- row.names(x)
nrow(IL1B.strat.markers)
nrow(x)
get_ontology(x, name = "IL1B expression level stratification full set", outdir = "Plaque_macrophages_results/", full_GSEA = F, universe = plaque.macrophages.seurat, plot.top.n = 10)

# Zero vs. High
IL1B.strat.markers <- FindMarkers(object = plaque.macrophages.seurat, min.pct = 0.25, logfc.threshold = 0.25, ident.2 = "High", ident.1 = "Zero", group.by = "IL1B_expr", only.pos = T)

bunchOfCustomPlots(object    = plaque.macrophages.seurat, 
                   features  = row.names(IL1B.strat.markers)[1:9], 
                   group.by  = "IL1B_expr", 
                   name      = "Plaque_macrophages_results/IL1B expression level stratification markers zero vs. high", 
                   pt.size   = 0, 
                   dot.scale = 10,
                   Vln.color = c("grey", "bisque", "coral", "firebrick"), 
                   ncol      = 3)

x <- subset(IL1B.strat.markers, p_val_adj < 0.05)
x$gene <- row.names(x)
nrow(IL1B.strat.markers)
nrow(x)
get_ontology(x, name = "IL1B expression level stratification zero vs. high padj 0.05", outdir = "Plaque_macrophages_results/", full_GSEA = F, universe = plaque.macrophages.seurat, plot.top.n = 10)

x <- IL1B.strat.markers
x$gene <- row.names(x)
nrow(IL1B.strat.markers)
nrow(x)
get_ontology(x, name = "IL1B expression level stratification zero vs. high full set", outdir = "Plaque_macrophages_results/", full_GSEA = F, universe = plaque.macrophages.seurat, plot.top.n = 10)

# High vs. Low
IL1B.strat.markers <- FindMarkers(object = plaque.macrophages.seurat, min.pct = 0.25, logfc.threshold = 0.25, ident.2 = "High", ident.1 = "Low", group.by = "IL1B_expr", only.pos = T)

bunchOfCustomPlots(object    = plaque.macrophages.seurat, 
                   features  = row.names(IL1B.strat.markers)[1:9], 
                   group.by  = "IL1B_expr", 
                   name      = "Plaque_macrophages_results/IL1B expression level stratification markers high vs. low", 
                   pt.size   = 0, 
                   dot.scale = 10,
                   Vln.color = c("grey", "bisque", "coral", "firebrick"), 
                   ncol      = 3)

x <- subset(IL1B.strat.markers, p_val_adj < 0.05)
x$gene <- row.names(x)
nrow(IL1B.strat.markers)
nrow(x)

x <- IL1B.strat.markers
x$gene <- row.names(x)
nrow(IL1B.strat.markers)
nrow(x)
get_ontology(x, name = "IL1B expression level stratification high vs low full set", outdir = "Plaque_macrophages_results/", full_GSEA = F, universe = plaque.macrophages.seurat, plot.top.n = 10)

# Low vs. High
IL1B.strat.markers <- FindMarkers(object = plaque.macrophages.seurat, min.pct = 0.25, logfc.threshold = 0.25, ident.1 = "High", ident.2 = "Low", group.by = "IL1B_expr", only.pos = T)

bunchOfCustomPlots(object    = plaque.macrophages.seurat, 
                   features  = row.names(IL1B.strat.markers)[1:9], 
                   group.by  = "IL1B_expr", 
                   name      = "Plaque_macrophages_results/IL1B expression level stratification markers low vs. high", 
                   pt.size   = 0, 
                   dot.scale = 10,
                   Vln.color = c("grey", "bisque", "coral", "firebrick"), 
                   ncol      = 3)

x <- subset(IL1B.strat.markers, p_val_adj < 0.05)
x$gene <- row.names(x)
nrow(IL1B.strat.markers)
nrow(x)
get_ontology(x, name = "IL1B expression level stratification low vs. high padj 0.05", outdir = "Plaque_macrophages_results/", full_GSEA = F, universe = plaque.macrophages.seurat, plot.top.n = 10)

x <- IL1B.strat.markers
x$gene <- row.names(x)
nrow(IL1B.strat.markers)
nrow(x)
get_ontology(x, name = "IL1B expression level stratification low vs. high full set", outdir = "Plaque_macrophages_results/", full_GSEA = F, universe = plaque.macrophages.seurat, plot.top.n = 10)



##==================================================================================
##==================================================================================
## Look at markers between the individual populations (comparing one vs. one)
dir.create("Plaque_macrophages_results/clusters/1_vs_1", showWarnings = F)
DefaultAssay(plaque.macrophages.seurat) <- "RNA"

# Populate a list with all possible combinations of populations and their markers
one_vs_one.markers <- list()
tmp.ident.list <- unique(Idents(plaque.macrophages.seurat))
for(theFirstIdent in unique(Idents(plaque.macrophages.seurat))){
  cat("Comparing", theFirstIdent, "with:\n")
  idents.not_self <- tmp.ident.list[!tmp.ident.list %in% theFirstIdent] # Don't compare an ident to itself
  for(theSecondIdent in idents.not_self){
    cat("\t", theSecondIdent,"...\n")
    one_vs_one.markers[[paste(theFirstIdent, "_", theSecondIdent, sep = "")]] <- FindMarkers(plaque.macrophages.seurat, 
                                                                                             ident.1         = theFirstIdent, 
                                                                                             ident.2         = theSecondIdent, 
                                                                                             min.pct         = 0.25, 
                                                                                             logfc.threshold = 0.25)
  }
  tmp.ident.list <- tmp.ident.list[-which(tmp.ident.list == theFirstIdent)] # Remove the current firstIdent from the list of second idents. This way we don't get the (useless) reciprocal comaprisons.
  cat("\n")
}
rm(tmp.ident.list)

# Write the markers to disk
lapply(seq_along(one_vs_one.markers), function(i) write.better.table(data = one_vs_one.markers[[i]], 
                                                                     file = paste("Plaque_macrophages_results/clusters/1_vs_1/", names(one_vs_one.markers)[[i]], " var genes.txt", sep = "")))

# Plot the markers
lapply(seq_along(one_vs_one.markers), function(i) customVln(object   = plaque.macrophages.seurat, 
                                                            name     = paste("Plaque_macrophages_results/clusters/1_vs_1/", names(one_vs_one.markers)[[i]], " top 9 markers violin plot.pdf", sep = ""), 
                                                            features = row.names(one_vs_one.markers[[i]])[1:9], 
                                                            ncol     = 3)
)


##==================================================================================
##==================================================================================
## Use new 10X - CEL-seq integrated populations
# Get new refined idents
int.idents <- Idents(integrated.mye.seurat)
plaque.macrophages.seurat <- AddMetaData(plaque.macrophages.seurat, metadata = int.idents[names(int.idents) %in% colnames(plaque.macrophages.seurat)], col.name = "method.int.idents")

# Plot
customUMAP(object     = plaque.macrophages.seurat, 
           group.by   = "method.int.idents", 
           cols       = M.int_refined.pop.colors, 
           title      = "43p refined populations", 
           legend.pos = "none", 
           plot.width = 15,
           file.name  = "Plaque_macrophages_results/UMAP 43p integrated populations.pdf" )

customUMAP(object     = plaque.macrophages.seurat, 
           group.by   = "method.int.idents", 
           cols       = M.int_refined.pop.colors, 
           title      = "43p refined populations", 
           legend.pos = "right", 
           plot.width = 15,
           file.name  = "Plaque_macrophages_results/UMAP 43p integrated populations legend.pdf" )



## Also check the pre-merged integreated clusters 2, 4, and 6
int.idents <- integrated.mye.seurat$method.int.idents.pre.merge
plaque.macrophages.seurat <- AddMetaData(plaque.macrophages.seurat, metadata = int.idents[names(int.idents) %in% colnames(plaque.macrophages.seurat)], col.name = "numbered.method.int.idents")

# Plot
customUMAP(object     = plaque.macrophages.seurat, 
           group.by   = "numbered.method.int.idents", 
           title      = "43p numbered populations", 
           legend.pos = "right", 
           plot.width = 15,
           file.name  = "Plaque_macrophages_results/UMAP 43p integrated numbered populations.pdf" )







