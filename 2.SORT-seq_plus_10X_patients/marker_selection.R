#=================================================================================================================================
## Define markers for the refined macrophage and monocyte populations
## Also define markers for the 'big 3' macrophage types (inflammatory, resident, foamy)
##================================================================================================================================
## per cluster
dir.create("marker_selection/clusters", recursive = T, showWarnings = F)

## Define marker genes per cluster
mye.patients.seurat.markers <- FindAllMarkers(object = mye.patients.seurat, only.pos = TRUE, min.pct = 0.25, thresh = 0.25)

# Save the top9 markers per cluster
sep.markers <- mye.patients.seurat.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
num.clus <- length(unique(Idents(mye.patients.seurat)))
for(i in levels(Idents(mye.patients.seurat))){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = mye.patients.seurat,
                     name      = paste("marker_selection/clusters/Cluster ",i, " markers", sep = ""),
                     Vln.width = 25, Vln.height = 50)
}

# Save markers per cluster to disk
for(i in levels(Idents(mye.patients.seurat))){
  x <- mye.patients.seurat.markers[which(mye.patients.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("marker_selection/clusters/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

##==================================================================================
##==================================================================================
## Get ontologies
# All markers
for(i in levels(Idents(mye.patients.seurat))){
  x <- mye.patients.seurat.markers[which(mye.patients.seurat.markers$cluster == i),-6]
  get_ontology(x, name = i, outdir = "marker_selection/clusters", universe = mye.patients.seurat, full_GSEA = F, volcano.plot = T) 
}

# Top 8
for(i in levels(Idents(mye.patients.seurat))){
  x <- mye.patients.seurat.markers[which(mye.patients.seurat.markers$cluster == i),-6]
  get_ontology(x, name = paste("Top 8", i, sep = " "), outdir = "marker_selection/clusters", universe = mye.patients.seurat, full_GSEA = F, volcano.plot = F, plot.top.n = 8) 
}


##==================================================================================
##==================================================================================
## Define unique markers
uniq.mye.patients.seurat.markers <- mye.patients.seurat.markers[!(duplicated(mye.patients.seurat.markers$gene) | duplicated(mye.patients.seurat.markers$gene, fromLast = T)),]

# Save the top9 markers per cluster
sep.markers <- uniq.mye.patients.seurat.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
num.clus <- length(unique(Idents(mye.patients.seurat)))
for(i in levels(Idents(mye.patients.seurat))){
  # Check if there are ANY uniqe markers to plot
  uniq.features <- as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"]))
  if(length(uniq.features != 0)){
    # Make feature, violin, and dot plots
    bunchOfCustomPlots(features  = uniq.features, 
                       object    = mye.patients.seurat,
                       name      = paste("marker_selection/clusters/Cluster ",i, " unique markers", sep = ""),
                       Vln.width = 25, Vln.height = 50)
    }
}

# Save markers per cluster to disk
for(i in levels(Idents(mye.patients.seurat))){
  x <- uniq.mye.patients.seurat.markers[which(uniq.mye.patients.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("marker_selection/clusters/Cluster ", i, " unique var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

# Unique markers
for(i in levels(Idents(mye.patients.seurat))){
  x <- uniq.mye.patients.seurat.markers[which(uniq.mye.patients.seurat.markers$cluster == i),-6]
  get_ontology(x, name = paste(i, ".unique_markers", sep = ""), outdir = "marker_selection/clusters", universe = mye.patients.seurat, full_GSEA = F, volcano.plot = T) 
}

# Top 8
for(i in levels(Idents(mye.patients.seurat))){
  x <- uniq.mye.patients.seurat.markers[which(uniq.mye.patients.seurat.markers$cluster == i),-6]
  get_ontology(x, name = paste("Top 8 ", i, ".unique_markers", sep = ""), outdir = "marker_selection/clusters", universe = mye.patients.seurat, full_GSEA = F, volcano.plot = T, plot.top.n = 8) 
}



##==================================================================================
##==================================================================================
## Per category
dir.create("marker_selection/meta_clusters", recursive = T, showWarnings = F)
ident.category <- as.vector(Idents(mye.patients.seurat))

# Rename idents
ident.category[grep("Monocytes",   ident.category)] <- "Monocytes"
ident.category[grep("Resident",    ident.category)] <- "Resident Macrophages"
ident.category[grep("Foamy Mac",   ident.category)] <- "Foamy Macrophages"
ident.category[grep("derived Mac", ident.category)] <- "Inflammatory Macrophages"

# Add to seurat object
mye.patients.seurat <- AddMetaData(mye.patients.seurat, ident.category, col.name = "meta_type")

# Set colors
mye.cat.colors        <- c("lightskyblue", "springgreen2", "firebrick3", "mediumslateblue")
names(mye.cat.colors) <- unique(mye.patients.seurat$meta_type)

# Plot the categories in a UMAP
DimPlot(mye.patients.seurat,
        label      = F, 
        pt.size    = 5,
        label.size = 10, 
        reduction  = "umap", 
        shuffle    = T, 
        group.by   = "meta_type",
        cols       = mye.cat.colors) + theme(legend.position  = "left",
                                                   panel.background = element_blank(),
                                                   legend.text      = element_text(size = 18, face = "bold"),
                                                   axis.text        = element_text(size = 14, face = "plain"),
                                                   plot.margin      = unit(c(1,1,1,4), units = "cm")
        ) 
ggsave("marker_selection/meta_clusters/Myeloid UMAP meta-clusters labels.pdf", width = 25, height = 15)

## Define marker genes per meta cluster
mye.patients.seurat.meta_markers <- data.frame()
for(theCat in unique(mye.patients.seurat$meta_type)){
  tmp.df         <- FindMarkers(object = mye.patients.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, group.by = "meta_type", ident.1 = theCat)
  tmp.df$cluster <- theCat
  mye.patients.seurat.meta_markers <- rbind(mye.patients.seurat.meta_markers, tmp.df)
}
mye.patients.seurat.meta_markers$gene <- row.names(mye.patients.seurat.meta_markers)

# Save the top9 meta_markers per meta cluster
sep.meta_markers <- mye.patients.seurat.meta_markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 meta_markers per meta cluster
num.clus <- length(unique(mye.patients.seurat$meta_type))
for(i in unique(mye.patients.seurat$meta_type)){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.meta_markers[sep.meta_markers$cluster == i, "gene"])), 
                     object    = mye.patients.seurat,
                     group.by  = "meta_type", 
                     Vln.color = mye.cat.colors, 
                     dot.scale = 15, Dot.width = 12,
                     name      = paste("marker_selection/meta_clusters/Meta Cluster ",i, " markers", sep = ""),
                     Vln.width = 20, Vln.height = 20)
}

# Save meta_markers per cluster to disk
for(i in unique(mye.patients.seurat$meta_type)){
  x <- mye.patients.seurat.meta_markers[which(mye.patients.seurat.meta_markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("marker_selection/meta_clusters/Meta Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

##==================================================================================
##==================================================================================
## Get ontologies
# All meta_markers
for(i in unique(mye.patients.seurat$meta_type)){
  x <- mye.patients.seurat.meta_markers[which(mye.patients.seurat.meta_markers$cluster == i),-6]
  get_ontology(x, name = i, outdir = "marker_selection/meta_clusters", universe = mye.patients.seurat, full_GSEA = F, volcano.plot = T) 
}

# Top 8
for(i in unique(mye.patients.seurat$meta_type)){
  x <- mye.patients.seurat.meta_markers[which(mye.patients.seurat.meta_markers$cluster == i),-6]
  get_ontology(x, name = paste("Top 8", i, sep = " "), outdir = "marker_selection/meta_clusters", universe = mye.patients.seurat, full_GSEA = F, volcano.plot = F, plot.top.n = 8) 
}


##==================================================================================
##==================================================================================
# Check some relevant genes in the meta clusters
## Re-plot some markers
customVln(features  = c("TREM2", "CD9", "GPNMB", "OLR1", "ABCA1", "ABCG1", "IL1B", "TNF", "CASP1", "CD14", "FCGR1A", "FCGR3A"), 
          object    = mye.patients.seurat,
          ncol      = 3,
          name      = "marker_selection/meta_clusters/resident - foam - inflammatory - general markers - names.pdf", 
          cols      = mye.cat.colors,
          group.by  = "meta_type",
          width = 20, height = 20
)


customVln(features = c("CASP1","IL1B","TLR4","ABCA1","OLR1", "TREM2","CASP4","KLF4","TNF","IL18","ABCG1", "CD9"), 
          object = mye.patients.seurat,
          name   = "marker_selection/meta_clusters/Foam vs. pro inf genes - names.pdf", 
          ncol   = 3,
          cols   = mye.cat.colors,
          group.by  = "meta_type",
          width  = 30, height = 20
)

customVln(features = c("IL1B", "FCGR3A"), 
          object = mye.patients.seurat,
          name   = "marker_selection/meta_clusters/zoom in IL1B FCGR3A.pdf", 
          ncol   = 2,
          cols   = mye.cat.colors,
          group.by  = "meta_type",
          width  = 10, height = 10
)

bunchOfCustomPlots(object    = mye.patients.seurat, 
                   features  = c("TREM2", "CD9", "GPNMB", "OLR1", "ABCA1", "ABCG1", "IL1B", "TNF", "CASP1", "CD14", "SELL", "FCGR3A"), 
                   name      = "marker_selection/meta_clusters/Marker roundup", 
                   Vln.stack = T, Vln.height = 20, Dot.width = 20, dot.scale = 15,
                   group.by  = "meta_type",
                   ncol      = 3)



