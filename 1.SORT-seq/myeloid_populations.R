## Analyses of myeloid populations
#==========================================================================================================
##=========================================================================================================
# Make directory for Myeloid_cell_subset_results
dir.create("Myeloid_cell_subset_results", showWarnings = FALSE)

# Extract myeloid clusters
pops <- unique(Idents(all.seur.combined))
CD68.pops <- factor(pops[grep("CD68", pops)])

# Make some plots
# FUT4 = CD15, CD11B = ITGAM
gen.mye.marks <- c("CD14", "CD68", "CD1C", "CD86", "KIT", "IL1B", "CASP1", "TLR4", "IL18", "TREM2", "OLR1", "ABCA1", "ITGAM", "CD33", "FUT4", "HLA-DRA", "ACTA2", "MYH11", "CD3E", "CD4") 

VlnPlot(all.seur.combined, idents = CD68.pops, features = gen.mye.marks, stack = T, same.y.lims = T, flip = T, fill.by = "ident") +
  theme(legend.position = "none")
ggsave("Myeloid_cell_subset_results/Marker genes all myeloid clusters violinplot.pdf", width = 10, height = 15)

RidgePlot(all.seur.combined, idents = CD68.pops, features = gen.mye.marks, stack = T, same.y.lims = T, fill.by = "ident") +
  theme(legend.position = "none")
ggsave("Myeloid_cell_subset_results/Marker genes all myeloid clusters Ridge ridgeplot.pdf", width = 20, height = 5)

DotPlot(all.seur.combined, idents = CD68.pops, features = gen.mye.marks, cluster.idents = T, dot.scale = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gsea()
ggsave("Myeloid_cell_subset_results/Marker genes all myeloid clusters dotplot.pdf", width = 15, height = 5)

#==========================================================================================================
##=========================================================================================================
## Continue with the macrophage populations
# Extract the mac clusters
M_clusters <- factor(CD68.pops[grep("macro|Foam", CD68.pops)])

# Make a subsetted seurat object
all.seur.combined.all_M_clusters <- subset(all.seur.combined, idents = M_clusters)

#==========================================================================================================
##=========================================================================================================
## Let's see if subclustering helps
## Renormalize and recluster
# run sctransform
all.seur.combined.all_M_clusters.reclustered <- SCTransform(all.seur.combined.all_M_clusters, vars.to.regress = vars.to.regress, verbose = FALSE)

# Run PCA and UMAP embedding
all.seur.combined.all_M_clusters.reclustered <- RunPCA(all.seur.combined.all_M_clusters.reclustered, verbose = T)
all.seur.combined.all_M_clusters.reclustered <- RunUMAP(all.seur.combined.all_M_clusters.reclustered, dims = 1:30, verbose = T)

# Clustering
set.seed(1)
all.seur.combined.all_M_clusters.reclustered <- FindNeighbors(all.seur.combined.all_M_clusters.reclustered, dims = 1:30, verbose = T, force.recalc = T)
all.seur.combined.all_M_clusters.reclustered <- FindClusters(all.seur.combined.all_M_clusters.reclustered, verbose = T, resolution = 0.75, random.seed = 666)
DimPlot(all.seur.combined.all_M_clusters.reclustered, label = TRUE, pt.size = 2, label.size = 10) + NoLegend() + theme(panel.background = element_blank(),
                                                                                                                   axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/UMAP reclustered.pdf", width = 10, height = 20)

VlnPlot(all.seur.combined.all_M_clusters.reclustered, features = gen.mye.marks, stack = T, same.y.lims = T, flip = T, fill.by = "ident") +
  theme(legend.position = "none")
ggsave("Myeloid_cell_subset_results/Marker genes reclustered macrophage clusters violinplot.pdf", width = 10, height = 15)


VlnPlot(all.seur.combined.all_M_clusters.reclustered, features = gen.mye.marks, stack = F, same.y.lims = T, ncol = 4, fill.by = "ident") +
  theme(legend.position = "none")
ggsave("Myeloid_cell_subset_results/Marker genes reclustered macrophage clusters violinplot no stack.pdf", width = 10, height = 15)

RidgePlot(all.seur.combined.all_M_clusters.reclustered, features = gen.mye.marks, stack = T, same.y.lims = T, fill.by = "ident") +
  theme(legend.position = "none")
ggsave("Myeloid_cell_subset_results/Marker genes reclustered macrophage clusters ridgeplot.pdf", width = 20, height = 5)

DotPlot(all.seur.combined.all_M_clusters.reclustered, features = gen.mye.marks, cluster.idents = T, dot.scale = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gsea()
ggsave("Myeloid_cell_subset_results/Marker genes reclustered macrophage clusters dotplot.pdf", width = 15, height = 5)

VlnPlot(all.seur.combined.all_M_clusters.reclustered, features = c("CD3E", "CD4", "CD8A"), stack = F, same.y.lims = T, ncol = 4, fill.by = "ident") +
  theme(legend.position = "none")
ggsave("Myeloid_cell_subset_results/Marker genes reclustered macrophage clusters T cell genes.pdf")

all.seur.combined.all_M_clusters <- all.seur.combined.all_M_clusters.reclustered
rm(all.seur.combined.all_M_clusters.reclustered)

##=========================================================================================================
##==========================================================================================================
## Inspect the mac clusters in isolation
## Scale all genes from RNA assay for diff gene calling and visualisation
# Set default to RNA
DefaultAssay(all.seur.combined.all_M_clusters) <- "RNA"

# Normalize
all.seur.combined.all_M_clusters <- NormalizeData(all.seur.combined.all_M_clusters, normalization.method = "LogNormalize")

# Scale
all.seur.combined.all_M_clusters <- ScaleData(all.seur.combined.all_M_clusters, vars.to.regress = vars.to.regress, features = row.names(all.seur.combined.all_M_clusters))

# Find marker genes
all.seur.combined.all_M_clusters.markers <- FindAllMarkers(all.seur.combined.all_M_clusters, only.pos = TRUE, min.pct = 0.25, thresh = 0.25)

#Save the top9 markers per cluster
sep.markers.M_clusters <- all.seur.combined.all_M_clusters.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

#plot the top9 markers per cluster
num.clus <- length(M_clusters)
for(i in unique(Idents(all.seur.combined.all_M_clusters))){
  VlnPlot(object = all.seur.combined.all_M_clusters,
          features = as.vector(unlist(sep.markers.M_clusters[sep.markers.M_clusters$cluster == i, "gene"]))
  ) + theme(panel.background = element_blank(),
            text             = element_text(size = 16, face = "bold"),
            axis.text        = element_text(size = 14, face = "plain")
  )
  ggsave(paste("Myeloid_cell_subset_results/mac cluster.",i, " markers violin.pdf", sep = ""), width = 10, height = 20)
}

#Plot top markers in a heatmap
top10 <- all.seur.combined.all_M_clusters.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(object = all.seur.combined.all_M_clusters, features = top10$gene, raster = F, label = F) + theme(legend.position  = "top",
                                                                                panel.background = element_blank(),
                                                                                text             = element_text(size = 16, face = "bold"),
                                                                                axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/mac clusters heatmap.pdf", width = 20, height =17)

#Write var genes per cluster to files
for (i in 0:(num.clus-1)){
  x <- all.seur.combined.all_M_clusters.markers[which(all.seur.combined.all_M_clusters.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("Myeloid_cell_subset_results/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}


VlnPlot(all.seur.combined.all_M_clusters, features = gen.mye.marks, stack = T, same.y.lims = T, flip = T, fill.by = "ident") +
  theme(legend.position = "none")
ggsave("Myeloid_cell_subset_results/Marker genes macrophage clusters violinplot.pdf", width = 10, height = 15)

VlnPlot(all.seur.combined.all_M_clusters, features = gen.mye.marks, stack = F, same.y.lims = T, ncol = 4, fill.by = "ident") +
  theme(legend.position = "none")
ggsave("Myeloid_cell_subset_results/Marker genes macrophage clusters violinplot no stack.pdf", width = 10, height = 15)

RidgePlot(all.seur.combined.all_M_clusters, features = gen.mye.marks, stack = T, same.y.lims = T, fill.by = "ident") +
  theme(legend.position = "none")
ggsave("Myeloid_cell_subset_results/Marker genes macrophage clusters Ridge ridgeplot.pdf", width = 20, height = 5)

DotPlot(all.seur.combined.all_M_clusters, features = gen.mye.marks, cluster.idents = T, dot.scale = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gsea()
ggsave("Myeloid_cell_subset_results/Marker genes macrophage clusters dotplot.pdf", width = 15, height = 5)

VlnPlot(all.seur.combined.all_M_clusters, features = c("CD3E", "CD4", "CD8A"), stack = F, same.y.lims = T, ncol = 4, fill.by = "ident") +
  theme(legend.position = "none")
ggsave("Myeloid_cell_subset_results/Marker genes  macrophage clusters T cell genes.pdf")


##==========================================================================================================
##==========================================================================================================
#Look at some genes
VlnPlot(object = all.seur.combined.all_M_clusters,
        features = c("CD4", "CD14", "CD68", "SPI1")
) + theme(panel.background = element_blank(),
          text             = element_text(size = 16, face = "bold"),
          axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/myeloid markers violin.pdf", width = 10, height = 20)

VlnPlot(object = all.seur.combined.all_M_clusters, 
        features = c("CASP1", "FOLR2", "TREM2", "IL1B", "TNF", "ABCG1"), ncol = 3
) + theme(panel.background = element_blank(),
          text             = element_text(size = 16, face = "bold"),
          axis.text        = element_text(size = 14, face = "plain", angle =45, hjust = 1),
          aspect.ratio = 1/2
)
ggsave("results/Mac markers.pdf", width = 10, height = 20)

VlnPlot(object = all.seur.combined.all_M_clusters,
        features = c("TNF", "CLEC10A", "IL1B", "CASP1","KLF4", "SELL", "S100A12", "OLR1","ABCA1","TREM2","MMP9","FABP4"), ncol = 4
) 
ggsave("Myeloid_cell_subset_results/Foam v proinf.pdf", width = 15, height = 15)


VlnPlot(object = all.seur.combined.all_M_clusters,
        features = c("CD200R1", "TGM2", "CD163")
) + theme(panel.background = element_blank(),
          text             = element_text(size = 16, face = "bold"),
          axis.text        = element_text(size = 14, face = "plain")
)

ggsave("Myeloid_cell_subset_results/m2 markers violin.pdf", width = 10, height = 20)

VlnPlot(object = all.seur.combined.all_M_clusters, 
        features = c("CD2","ITGAM", "PECAM1", "NCAM1", "SELL", "CSF1R")
) + theme(panel.background = element_blank(),
          text             = element_text(size = 16, face = "bold"),
          axis.text        = element_text(size = 14, face = "plain")
)

ggsave("Myeloid_cell_subset_results/monocyte markers violin.pdf", width = 10, height = 20)

VlnPlot(object = all.seur.combined.all_M_clusters,
        features = c("CD14", "CD68", "ITGAM", "CD33", "FUT4", "HLA-DRA"), ncol = 3
)
ggsave("Myeloid_cell_subset_results/MDSC genes.pdf", width = 25, height = 15)


#SMC Monocyte genes
VlnPlot(object = all.seur.combined.all_M_clusters, features = c("ACTA2", "MYH11", "CD14", "CD68", "KLF4", "CD4", "TAGLN", "ADGRE1", "LGALS3")) + 
  theme(panel.background = element_blank(),
        text             = element_text(size = 16, face = "bold"),
        axis.text        = element_text(size = 14, face = "plain",angle=45,hjust=1)
  )
ggsave("Myeloid_cell_subset_results/SMC mono genes violin.pdf", width = 10, height = 20)

ACTA2 <- GetAssayData(all.seur.combined.all_M_clusters[c("ACTA2"),])
all.seur.combined.all_M_clusters@meta.data$ACTA2_level <- as.vector(ACTA2 != 0)
VlnPlot(object = all.seur.combined.all_M_clusters, split.by = "ACTA2_level", split.plot = T,
        features = c("ACTA2", "CD68"))
ggsave("Myeloid_cell_subset_results/SMC mono genes violin ACTA2 Split.pdf", width = 10, height = 20)


VlnPlot(object = all.seur.combined.all_M_clusters, cols =c("#F2756D", "#B49F31", "#28B24B", "#00FF00"), features =  c("HLA-A","HLA-B","HLA-C","HLA-DMA",
                                                                                                                  "HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1",
                                                                                                                  "HLA-DPB1","HLA-DQA1","HLA-DRA","HLA-DRB1",
                                                                                                                  "HLA-DRB5","HLA-E","HLA-F","HLA-F-AS1"), ncol = 4) + 
  theme(panel.background = element_blank(),
        text             = element_text(size = 16, face = "bold"),
        axis.text        = element_text(size = 14, face = "plain",angle=45,hjust=1)
  )
ggsave("Myeloid_cell_subset_results/HLA genes myeloid clusters.pdf", width = 10, height = 20)

#Inflammasome genes
VlnPlot(object = all.seur.combined.all_M_clusters, features = c("NLRP3", "IL1R1", "IL18", "IL1B", "NEK7", "CASP1", "PYCARD", "CASP4", "NFKB1")) + 
  theme(panel.background = element_blank(),
        text             = element_text(size = 16, face = "bold"),
        axis.text        = element_text(size = 14, face = "plain",angle=45,hjust=1)
  )
ggsave("Myeloid_cell_subset_results/inflammasome violin.pdf", width = 10, height = 20)


#Rename clusters, replot UMAP, and rename marker genes and lists
new.ident <-  c("CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells", "CD68+CD3E+CD4+IL1B+ Inflammatory macrophages and T cells")
length(new.ident)
length(unique(new.ident))
names(new.ident)  <- 0:3
new.ident
all.seur.combined.all_M_clusters <- RenameIdents(all.seur.combined.all_M_clusters, new.ident)
DimPlot(all.seur.combined.all_M_clusters, label = T, pt.size = 2, label.size = 10) + theme(legend.position  = "none",
                                                                                       panel.background = element_blank(),
                                                                                       text             = element_text(size = 14, face = "bold"),
                                                                                       axis.text        = element_text(size = 14, face = "plain")
)
ggsave("Myeloid_cell_subset_results/UMAP clusters labels.pdf", width = 15, height = 15)

#Get cluster colors for subsequent plots
t                        <- DimPlot(all.seur.combined.all_M_clusters)
tbuild                   <- ggplot_build(t)
tdata                    <- tbuild$data[[1]] 
tdata                    <- tdata[order(tdata$group), ]
cluster_colours.M        <- unique(tdata$colour)
names(cluster_colours.M) <- new.ident

# Find marker genes
all.seur.combined.all_M_clusters.markers <- FindAllMarkers(all.seur.combined.all_M_clusters, only.pos = TRUE, min.pct = 0.25, 
                                                       thresh = 0.25)

#Write var genes per cluster to files
for(i in levels(Idents(all.seur.combined.all_M_clusters))){
  x <- all.seur.combined.all_M_clusters.markers[which(all.seur.combined.all_M_clusters.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("Myeloid_cell_subset_results/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

# Save the object
saveRDS(all.seur.combined.all_M_clusters, "Seurat Objects/myeloid.seurat.filtered.clustered.RDS")
