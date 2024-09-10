#=========================================================================================================================
## Integrate the myeloid cells with the 43 patient CEL-seq myeloid cells
##========================================================================================================================
dir.create("myeloid_43p_celseq_integration/First round of clustering/clusters", showWarnings = F, recursive = T)

## Load the data
mye.43p.seurat <- readRDS(file = "../SORTseq_Patients/Seurat Objects/myeloid.seurat.filtered.clustered.RDS")
mye.43p.seurat
unique(Idents(mye.43p.seurat))

## Clean up
# Remove the neutrophil cluster
mye.43p.seurat <- subset(mye.43p.seurat, idents = as.vector(unique(Idents(mye.43p.seurat))[grep("neutro", unique(Idents(mye.43p.seurat)), invert = T)]))

# Remove trace contaminants
DefaultAssay(mye.43p.seurat) <- "RNA"
mye.43p.seurat <- subset(mye.43p.seurat, subset = CD3E == 0 & CD8A == 0 & NCAM1 == 0 & CD79A == 0)

## Prep for integration
# Select the corect assay
DefaultAssay(plaque.macrophages.seurat) <- "RNA"

# Make a list of our objects to integrate
# Start from the macs only, and add in the monocytes later. Since the 43p data does not contain monos, this prevents overfitting and saves time while redefining the clusters
integrated.mye.seurat <- list(plaque.macrophages.seurat, mye.43p.seurat)

# Define var features
integrated.mye.seurat <- lapply(X = integrated.mye.seurat, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = integrated.mye.seurat)

## Integrate the data
# First select anchors for integration (by CCA reduction)
anchors <- FindIntegrationAnchors(object.list = integrated.mye.seurat, anchor.features = features)

# And do the actual integrating
integrated.mye.seurat <- IntegrateData(anchorset = anchors)

# Show stats
integrated.mye.seurat

#=========================================================================================================================
## Update metadata
head(integrated.mye.seurat@meta.data)
tail(integrated.mye.seurat@meta.data)

# Fix tissue column from 43p to 'plaque'
integrated.mye.seurat@meta.data[is.na(integrated.mye.seurat$Tissue),"Tissue"] <- "plaque"

# Fix sex column from 10X to 'male'
integrated.mye.seurat@meta.data[is.na(integrated.mye.seurat$Sex),"Sex"] <- "male"

# Add a new column 'Method' to distinguish the two datasets
method <- rep("CEL-seq", dim(integrated.mye.seurat)[2])
method[grep("^P", integrated.mye.seurat$Patient)] <- "10X"
integrated.mye.seurat <- AddMetaData(integrated.mye.seurat, method, col.name = "Method")

## Add a new column 'original.refined.ident' to capture the original idents of the unintegrated objects.
# Set up destination dataframe
ref.ident  <- data.frame(row.names = row.names(integrated.mye.seurat@meta.data))

# Fetch idents and stick 'em in a dataframe for merging
idents.43p <- data.frame(refined.ident = Idents(mye.43p.seurat))
idents.10X <- data.frame(refined.ident = Idents(plaque.macrophages.seurat))
idents.com <- rbind(idents.43p, idents.10X)

# Merge the frames together
ref.ident <- merge(ref.ident, idents.com, by = 0)
row.names(ref.ident) <- ref.ident$Row.names
ref.ident$Row.names <- NULL
tail(ref.ident)

# And add the metadata
integrated.mye.seurat <- AddMetaData(integrated.mye.seurat, metadata = ref.ident, col.name = "original.refined.ident")

## Expand the colorscheme to include the 43p clusters
int.mac.pop.colors <- M_refined.pop.colors

# Add 43p clusters
tmp.names                 <- names(int.mac.pop.colors)
int.mac.pop.colors        <- c(int.mac.pop.colors, c("aquamarine4","brown2","magenta3","magenta1", "brown3"))
names(int.mac.pop.colors) <- c(tmp.names, levels(unique(Idents(mye.43p.seurat))))

#=========================================================================================================================
## Run the embeddings, cluster, and have a look
##========================================================================================================================
# Run dimensionality reductions
integrated.mye.seurat <- ScaleData(object = integrated.mye.seurat, verbose = F)
integrated.mye.seurat <- RunPCA(   object = integrated.mye.seurat, verbose = F)
integrated.mye.seurat <- RunUMAP(  object = integrated.mye.seurat, reduction = "pca", dims = 1:30)

# Cluster
set.seed(1)
integrated.mye.seurat <- FindNeighbors(object = integrated.mye.seurat, reduction = "pca", dims = 1:30, force.recalc = T, verbose = T)
integrated.mye.seurat <- FindClusters( object = integrated.mye.seurat, resolution = 1.25, verbose = T, random.seed = 666)

## Visualize
# Plot Clusters
customUMAP(object     = integrated.mye.seurat, shuffle = T,
           title      = "Initial clustering", 
           legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/First round of clustering/Integration intitial clustering UMAP.pdf", 
           plot.width = 15)

# Refine manually
plot         <- DimPlot(integrated.mye.seurat, reduction = "umap")
select.cells <- CellSelector(plot = plot)
Idents(integrated.mye.seurat, cells = select.cells) <- "2.1"

plot         <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident")
select.cells <- CellSelector(plot = plot)
Idents(integrated.mye.seurat, cells = select.cells) <- "7.1"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", shuffle = T, title = "Integrated clusters", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/First round of clustering/Integration cluster UMAP.pdf", 
           plot.width = 15)

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "Patient", shuffle = T, title = "Patient distribution", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/First round of clustering/Integration patient UMAP.pdf", 
           plot.width = 15)

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "Method", shuffle = T, title = "Method distribution", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/First round of clustering/Integration method UMAP.pdf", 
           plot.width = 15)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "Patient", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "bottom") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "Method", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "bottom") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p3 <- DimPlot(integrated.mye.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p1 + p2 + p3
ggsave("myeloid_43p_celseq_integration/First round of clustering/Integration UMAP.pdf", width = 30, height = 15)

## Plot the original 10X and 43p idents
customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "Source clusters", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/First round of clustering/Integration original cluster UMAP.pdf", 
           plot.width = 15)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.mac.pop.colors[c(-1,-2,-3)], pt.size = 3) +
  ggtitle("Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T, pt.size = 3) +
  ggtitle("Integrated clusters") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("myeloid_43p_celseq_integration/First round of clustering/Integration original vs. integrated cluster UMAP.pdf", width = 20, height = 10)

# Highlight each source cluster
dir.create("myeloid_43p_celseq_integration/First round of clustering/Original source cluster highlights", showWarnings = F, recursive = T)
for(theIdent in unique(integrated.mye.seurat$original.refined.ident)){
  customUMAP(object          = integrated.mye.seurat, reduction = "umap", label = F, shuffle = T, title = theIdent, legend.pos = "none", 
             file.name       = paste("myeloid_43p_celseq_integration/First round of clustering/Original source cluster highlights/", theIdent, " UMAP.pdf", sep = ""),
             sizes.highlight = 4,
             cells.highlight = row.names(ref.ident)[which(ref.ident$refined.ident == theIdent)])
}

## Color only one source at a time
# 43p colored
int.43p.pop.colors <- int.mac.pop.colors
int.43p.pop.colors[grep("^CD68", names(int.43p.pop.colors), invert = T)] <- "grey"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "CEL-seq Source population identities", legend.pos = "right", 
           cols       = int.43p.pop.colors[c(-1,-2,-3)],
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/First round of clustering/Integration original 43p cluster UMAP.pdf", 
           plot.width = 15)

# 10X colored
int.10X.pop.colors <- int.mac.pop.colors
int.10X.pop.colors[grep("^CD68", names(int.10X.pop.colors), invert = F)] <- "grey"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "10X Source population identities", legend.pos = "right", 
           cols       = int.10X.pop.colors[c(-1,-2,-3)],
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/First round of clustering/Integration original 10X cluster UMAP.pdf", 
           plot.width = 15)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.43p.pop.colors[c(-1,-2,-3)], pt.size = 5) +
  ggtitle("CEL-seq Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.10X.pop.colors[c(-1,-2,-3)], pt.size = 5) +
  ggtitle("10X Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("myeloid_43p_celseq_integration/First round of clustering/Integration 43p vs. 10X clusters UMAP.pdf", width = 20, height = 10)

## Save the integrated seurat object
saveRDS(integrated.mye.seurat, "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")
#integrated.mye.seurat <- readRDS(file = "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")


#=========================================================================================================================
## Define the clusters: plot markers, ontologies, and gene of interest
##========================================================================================================================
## Scale all genes from RNA assay for diff gene calling and visualisation
# Set default to RNA
DefaultAssay(integrated.mye.seurat) <- "RNA"

# Normalize
integrated.mye.seurat <- NormalizeData(integrated.mye.seurat, normalization.method = "LogNormalize")

# Scale
integrated.mye.seurat <- ScaleData(integrated.mye.seurat, features = row.names(integrated.mye.seurat))

##========================================================================================================================
## Define marker genes per cluster
integrated.mye.seurat.markers <- FindAllMarkers(object = integrated.mye.seurat, only.pos = TRUE, min.pct = 0.25, thresh = 0.25)

# Save the top9 markers per cluster
sep.markers <- integrated.mye.seurat.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
num.clus <- length(unique(Idents(integrated.mye.seurat)))
for(i in levels(Idents(integrated.mye.seurat))){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = integrated.mye.seurat,
                     name      = paste("myeloid_43p_celseq_integration/First round of clustering/clusters/Cluster ",i, " markers", sep = ""),
                     
                     Vln.width = 25, Vln.height = 25)
}

# Save markers per cluster to disk
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("myeloid_43p_celseq_integration/First round of clustering/clusters/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

##========================================================================================================================
# Look at some genes
bunchOfCustomPlots(features   = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                  "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                  "CD79A", "ACTA2", "CD34", "KIT", "OLR1" ), 
                   object     = integrated.mye.seurat,
                   ncol       = 5,
                   name       = "myeloid_43p_celseq_integration/First round of clustering/Cell type defining genes",
                   Vln.width  = 25, 
                   Vln.height = 15
)

# SMC Monocyte genes
bunchOfCustomPlots(features = c("ACTA2", "MYH11", "CD14", "CD68", "KLF4", "CD4", "TAGLN", "ADGRE1", "LGALS3"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/First round of clustering/SMC mono genes",
                   Vln.width = 25, Vln.height = 15
)

#Extra mono cell markers
bunchOfCustomPlots(features = c("FCGR3A", "CCR2", "CCR5", "SELL", "ITGAX", "CX3CR1"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/First round of clustering/Extra mono marker genes",
)

# Myeloid Markers
bunchOfCustomPlots(features = c("CD200R1", "TGM2", "CD163"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/First round of clustering/M2 genes")

bunchOfCustomPlots(features = c("CD2","ITGAM", "PECAM1", "NCAM1", "SELL", "CSF1R"), 
                   object    = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/First round of clustering/Monocyte genes",
)

bunchOfCustomPlots(features = c("CD68","CD14","ITGAM","CD1C", "CLEC10A", "FCER1A", "CD2", "CD3E", "CD4"),
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/First round of clustering/Myeloid genes", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("CASP1","IL1B","TLR4","ABCA1","OLR1", "TREM2","CASP4","KLF4","TNF","IL18","ABCG1", "CD9"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/First round of clustering/Foam vs. pro inf genes", 
                   ncol     = 3,
                   Vln.width = 30, Vln.height = 15
)

bunchOfCustomPlots(features = c("CD14","CD68","CD86","CD1C", "FCGR2B", "CLEC10A", "FCER1A", "HLA-DQA2", "HLA-DQB2"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/First round of clustering/DC genes", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("NLRP3", "IL1B", "CASP1", "CASP4", "IL18", "PYCARD"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/First round of clustering/Inflammasome genes", 
                   ncol     = 3,
)

bunchOfCustomPlots(features = c("CD14", "FCGR3A", "FCGR1A", "CCR2", "PADI4", "SECISBP2L", "CD86", "ITGAX", "LYZ", "CD33", "CD36", "ITGAM"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/First round of clustering/Mono subtype genes CD14 CD16 CD64 CCR2 PADI4 SLAN CD86 CD11c LYZ CD33 CD36 CD11b", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("NR4A1"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/First round of clustering/NR4A1",
                   ncol     = 1,
)

bunchOfCustomPlots(features  = c("CD14","FCGR3A", "CD68", 'CD86', "SELL", "CSF1R", "ITGAM", "ITGAX", "PECAM1"), 
                   object    = integrated.mye.seurat,
                   ncol      = 3, Vln.stack = T,
                   name      = "myeloid_43p_celseq_integration/First round of clustering/macro and mono genes",
                   Vln.width = 20, Vln.height = 20
)

bunchOfCustomPlots(features  = c("TREM2", "CD9", "GPNMB", "OLR1", "ABCA1", "ABCG1", "IL1B", "TNF", "CASP1", "CD14", "FCGR1A", "FCGR3A"), 
                   object    = integrated.mye.seurat,
                   ncol      = 3,
                   name      = "myeloid_43p_celseq_integration/First round of clustering/resident - foam - inflammatory - general markers",
                   Vln.width = 20, Vln.height = 20
)


bunchOfCustomPlots(features = c("TIMD4", "CCR2", "LYVE1", "FOLR2"), 
                   object   = integrated.mye.seurat, assay = "RNA",
                   name     = "myeloid_43p_celseq_integration/First round of clustering/core_Gene_macs",
                   ncol     = 4, 
                   Vln.width = 25
)

bunchOfCustomPlots(features = c("CX3CR1", "LYVE1", "MRC1", "F13A1", "GAS6", "FOLR2", "TREM2","CD9", "IL1B"), 
                   object   = integrated.mye.seurat, assay = "RNA",
                   name     = "myeloid_43p_celseq_integration/First round of clustering/resident_macs",
                   ncol     = 3, 
                   Vln.width = 25
)

customFeature(object = integrated.mye.seurat, features = c("MKI67"), name = "myeloid_43p_celseq_integration/First round of clustering/MKI67.pdf")


##==================================================================================
## Get ontologies
# All markers
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  get_ontology(x, name = i, outdir = "myeloid_43p_celseq_integration/First round of clustering/clusters", universe = integrated.mye.seurat, full_GSEA = F, volcano.plot = T) 
}

# Top 8
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  get_ontology(x, name = paste("Top 8", i, sep = " "), outdir = "myeloid_43p_celseq_integration/First round of clustering/clusters", universe = integrated.mye.seurat, full_GSEA = F, volcano.plot = F, plot.top.n = 8) 
}


##==================================================================================
##==================================================================================
## Second round of clustering based on info gleaned from first round!
##==================================================================================
##==================================================================================
dir.create("myeloid_43p_celseq_integration/Second round of clustering/clusters", showWarnings = F, recursive = T)

## Visualize
# Plot Clusters
customUMAP(object     = integrated.mye.seurat, reduction = "umap", shuffle = T, title = "Initial clustering UMAP", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Second round of clustering/Integration original 43p cluster UMAP.pdf", 
           plot.width = 15)

# Refine manually
plot         <- DimPlot(integrated.mye.seurat, reduction = "umap")
select.cells <- CellSelector(plot = plot)
Idents(integrated.mye.seurat, cells = select.cells) <- "2"

plot         <- DimPlot(integrated.mye.seurat, reduction = "umap")
select.cells <- CellSelector(plot = plot)
Idents(integrated.mye.seurat, cells = select.cells) <- "0"

plot         <- DimPlot(integrated.mye.seurat, reduction = "umap")
select.cells <- CellSelector(plot = plot)
Idents(integrated.mye.seurat, cells = select.cells) <- "7.1"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", shuffle = T, title = "Initial clustering UMAP", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Second round of clustering/Integration cluster UMAP.pdf", 
           plot.width = 15)

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "Patient", shuffle = T, title = "Initial clustering UMAP", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Second round of clustering/Integration patient UMAP.pdf", 
           plot.width = 15)

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "Method", shuffle = T, title = "Initial clustering UMAP", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Second round of clustering/Integration method UMAP.pdf", 
           plot.width = 15)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "Patient", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "bottom") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "Method", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "bottom") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p3 <- DimPlot(integrated.mye.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p1 + p2 + p3
ggsave("myeloid_43p_celseq_integration/Second round of clustering/Integration UMAP.pdf", width = 30, height = 15)

## Plot the original 10X and 43p idents
customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "Source population identities", legend.pos = "right", 
           cols       = int.mac.pop.colors[c(-1,-2,-3)], 
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/Second round of clustering/Integration original cluster UMAP.pdf", 
           plot.width = 15)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.mac.pop.colors[c(-1,-2,-3)], pt.size = 3) +
  ggtitle("Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T, pt.size = 3) +
  ggtitle("Integrated clusters") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("myeloid_43p_celseq_integration/Second round of clustering/Integration original vs. integrated cluster UMAP.pdf", width = 20, height = 10)

# Highlight each source cluster
dir.create("myeloid_43p_celseq_integration/Second round of clustering/Original source cluster highlights", showWarnings = F, recursive = T)
for(theIdent in unique(integrated.mye.seurat$original.refined.ident)){
  customUMAP(object          = integrated.mye.seurat, reduction = "umap", label = F, shuffle = T, title = theIdent, legend.pos = "none", 
             file.name       = paste("myeloid_43p_celseq_integration/Second round of clustering/Original source cluster highlights/", theIdent, " UMAP.pdf", sep = ""),
             sizes.highlight = 4,
             cells.highlight = row.names(ref.ident)[which(ref.ident$refined.ident == theIdent)])
}

## Color only one source at a time
# 43p colored
int.43p.pop.colors <- int.mac.pop.colors
int.43p.pop.colors[grep("^CD68", names(int.43p.pop.colors), invert = T)] <- "grey"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "CEL-seq Source population identities", legend.pos = "right", 
           cols       = int.43p.pop.colors[c(-1,-2,-3)],
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/Second round of clustering/Integration original 43p cluster UMAP.pdf", 
           plot.width = 15)

# 10X colored
int.10X.pop.colors <- int.mac.pop.colors
int.10X.pop.colors[grep("^CD68", names(int.10X.pop.colors), invert = F)] <- "grey"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "10X Source population identities", legend.pos = "right", 
           cols       = int.10X.pop.colors[c(-1,-2,-3)],
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/Second round of clustering/Integration original 10X cluster UMAP.pdf", 
           plot.width = 15)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.43p.pop.colors[c(-1,-2,-3)], pt.size = 5) +
  ggtitle("CEL-seq Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.10X.pop.colors[c(-1,-2,-3)], pt.size = 5) +
  ggtitle("10X Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("myeloid_43p_celseq_integration/Second round of clustering/Integration 43p vs. 10X clusters UMAP.pdf", width = 20, height = 10)

## Save the integrated seurat object
saveRDS(integrated.mye.seurat, "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")
#integrated.mye.seurat <- readRDS(file = "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")


#=========================================================================================================================
## Define the clusters: plot markers, ontologies, and gene of interest
##========================================================================================================================
## Scale all genes from RNA assay for diff gene calling and visualisation
# Set default to RNA
DefaultAssay(integrated.mye.seurat) <- "RNA"

# Normalize
integrated.mye.seurat <- NormalizeData(integrated.mye.seurat, normalization.method = "LogNormalize")

# Scale
integrated.mye.seurat <- ScaleData(integrated.mye.seurat, features = row.names(integrated.mye.seurat))

##========================================================================================================================
## Define marker genes per cluster
integrated.mye.seurat.markers <- FindAllMarkers(object = integrated.mye.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save the top9 markers per cluster
sep.markers <- integrated.mye.seurat.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
num.clus <- length(unique(Idents(integrated.mye.seurat)))
for(i in levels(Idents(integrated.mye.seurat))){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = integrated.mye.seurat,
                     name      = paste("myeloid_43p_celseq_integration/Second round of clustering/clusters/Cluster ",i, " markers", sep = ""),
                     
                     Vln.width = 25, Vln.height = 25)
}

# Save markers per cluster to disk
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("myeloid_43p_celseq_integration/Second round of clustering/clusters/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

##========================================================================================================================
# Look at some genes
bunchOfCustomPlots(features   = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                  "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                  "CD79A", "ACTA2", "CD34", "KIT", "OLR1" ), 
                   object     = integrated.mye.seurat,
                   ncol       = 5,
                   name       = "myeloid_43p_celseq_integration/Second round of clustering/Cell type defining genes",
                   Vln.width  = 25, 
                   Vln.height = 15
)

# SMC Monocyte genes
bunchOfCustomPlots(features = c("ACTA2", "MYH11", "CD14", "CD68", "KLF4", "CD4", "TAGLN", "ADGRE1", "LGALS3"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Second round of clustering/SMC mono genes",
                   Vln.width = 25, Vln.height = 15
)

#Extra mono cell markers
bunchOfCustomPlots(features = c("FCGR3A", "CCR2", "CCR5", "SELL", "ITGAX", "CX3CR1"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Second round of clustering/Extra mono marker genes",
)

# Myeloid Markers
bunchOfCustomPlots(features = c("CD200R1", "TGM2", "CD163"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Second round of clustering/M2 genes")

bunchOfCustomPlots(features = c("CD2","ITGAM", "PECAM1", "NCAM1", "SELL", "CSF1R"), 
                   object    = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Second round of clustering/Monocyte genes",
)

bunchOfCustomPlots(features = c("CD68","CD14","ITGAM","CD1C", "CLEC10A", "FCER1A", "CD2", "CD3E", "CD4"),
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Second round of clustering/Myeloid genes", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("CASP1","IL1B","TLR4","ABCA1","OLR1", "TREM2","CASP4","KLF4","TNF","IL18","ABCG1", "CD9"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Second round of clustering/Foam vs. pro inf genes", 
                   ncol     = 3,
                   Vln.width = 30, Vln.height = 15
)

bunchOfCustomPlots(features = c("CD14","CD68","CD86","CD1C", "FCGR2B", "CLEC10A", "FCER1A", "HLA-DQA2", "HLA-DQB2"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Second round of clustering/DC genes", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("NLRP3", "IL1B", "CASP1", "CASP4", "IL18", "PYCARD"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Second round of clustering/Inflammasome genes", 
                   ncol     = 3,
)

bunchOfCustomPlots(features = c("CD14", "FCGR3A", "FCGR1A", "CCR2", "PADI4", "SECISBP2L", "CD86", "ITGAX", "LYZ", "CD33", "CD36", "ITGAM"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Second round of clustering/Mono subtype genes CD14 CD16 CD64 CCR2 PADI4 SLAN CD86 CD11c LYZ CD33 CD36 CD11b", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("NR4A1"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Second round of clustering/NR4A1",
                   ncol     = 1,
)

bunchOfCustomPlots(features  = c("CD14","FCGR3A", "CD68", 'CD86', "SELL", "CSF1R", "ITGAM", "ITGAX", "PECAM1"), 
                   object    = integrated.mye.seurat,
                   ncol      = 3, Vln.stack = T,
                   name      = "myeloid_43p_celseq_integration/Second round of clustering/macro and mono genes",
                   Vln.width = 20, Vln.height = 20
)

bunchOfCustomPlots(features  = c("TREM2", "CD9", "GPNMB", "LYVE1", "MRC1", "FOLR2" ,"OLR1", "ABCA1", "ABCG1", "IL1B", "TNF", "CASP1", "CD14", "FCGR1A", "FCGR3A"), 
                   object    = integrated.mye.seurat,
                   ncol      = 3,
                   name      = "myeloid_43p_celseq_integration/Second round of clustering/LAM - resident - foam - inflammatory - general markers",
                   Vln.width = 20, Vln.height = 20
)


bunchOfCustomPlots(features = c("TIMD4", "CCR2", "LYVE1", "FOLR2"), 
                   object   = integrated.mye.seurat, assay = "RNA",
                   name     = "myeloid_43p_celseq_integration/Second round of clustering/core_Gene_macs",
                   ncol     = 4, 
                   Vln.width = 25
)

bunchOfCustomPlots(features = c("CX3CR1", "LYVE1", "MRC1", "F13A1", "GAS6", "FOLR2", "TREM2","CD9", "IL1B"), 
                   object   = integrated.mye.seurat, assay = "RNA",
                   name     = "myeloid_43p_celseq_integration/Second round of clustering/resident_macs",
                   ncol     = 3, 
                   Vln.width = 25
)

customFeature(object = integrated.mye.seurat, features = c("MKI67"), name = "myeloid_43p_celseq_integration/Second round of clustering/MKI67.pdf")


##==================================================================================
## Get ontologies
# All markers
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  get_ontology(x, name = i, outdir = "myeloid_43p_celseq_integration/Second round of clustering/clusters", universe = integrated.mye.seurat, full_GSEA = F, volcano.plot = T) 
}

# Top 8
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  get_ontology(x, name = paste("Top 8", i, sep = " "), outdir = "myeloid_43p_celseq_integration/Second round of clustering/clusters", universe = integrated.mye.seurat, full_GSEA = F, volcano.plot = F, plot.top.n = 8) 
}


##==================================================================================
## Marker gene pathways heatmap
ont.M <- list()
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  ont.M[[i]] <- get_ontology(res = x, name = i, outdir = "myeloid_43p_celseq_integration/Second round of clustering/clusters", universe = integrated.mye.seurat, full_GSEA = F, return.data = T, plot.top.n = 10)
}

# Plot differential pathways from top 30 per population
# Collect all data into a neat 'frame
ont.M.df <- merge(x        = ont.M[[1]]$Pathways, 
                  y        = ont.M[[2]]$Pathways, 
                  by       = "name", 
                  all      = T, 
                  suffixes = paste(".", names(ont.M)[1:2], sep = ""))

for(i in names(ont.M[3:length(names(ont.M))])){
  ont.M.df <- merge(x        = ont.M.df, 
                    y        = ont.M[[i]]$Pathways, 
                    by       = "name", 
                    all      = T)
  
  colnames(ont.M.df)[(length(colnames(ont.M.df))-2):length(colnames(ont.M.df))] <- paste(colnames(ont.M.df)[(length(colnames(ont.M.df))-2):length(colnames(ont.M.df))], i, sep = ".")
}

# Clean up
ont.M.df[is.na(ont.M.df)] <- 1
row.names(ont.M.df)       <- ont.M.df$name
ont.M.df$name             <- NULL

# Keep only padj columns
ont.M.df <- ont.M.df[,seq(from = 1, to = length(colnames(ont.M.df)), by = 3)]
head(ont.M.df)            # Neat B)

# Log transform the data
log10.ont.M.df                        <- -log10(ont.M.df)
log10.ont.M.df[log10.ont.M.df == Inf] <- 0
head(log10.ont.M.df)

# Keep top 10 rows per population
top10.log10.ont.M.df <- log10.ont.M.df[order(log10.ont.M.df[,1], decreasing = T),][1:10,]
for(i in length(colnames(ont.M.df))){
  top10.log10.ont.M.df <- rbind(top10.log10.ont.M.df, log10.ont.M.df[order(log10.ont.M.df[,i], decreasing = T),][1:10,])
}
top10.log10.ont.M.df <- top10.log10.ont.M.df[grep("1$", row.names(top10.log10.ont.M.df), invert = T), ]

# Keep top 30 rows per population
top30.log10.ont.M.df <- log10.ont.M.df[order(log10.ont.M.df[,1], decreasing = T),][1:30,]
for(i in length(colnames(ont.M.df))){
  top30.log10.ont.M.df <- rbind(top30.log10.ont.M.df, log10.ont.M.df[order(log10.ont.M.df[,i], decreasing = T),][1:30,])
}
top30.log10.ont.M.df <- top30.log10.ont.M.df[grep("1$", row.names(top30.log10.ont.M.df), invert = T), ]

## And plot
# Absolute -log10 padj values
pheatmap(log10.ont.M.df,
         show_rownames = F,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 1, 
         border_color = NA, 
         color = colorSpacer(startcolor = "white", middlecolors = "goldenrod", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Second round of clustering/clusters/Differential pathways.pdf")

pheatmap(log10.ont.M.df,
         kmeans_k = 5,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA, 
         color = colorSpacer(startcolor = "white", middlecolors = "goldenrod", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Second round of clustering/clusters/Differential pathways k means.pdf")

pheatmap(top10.log10.ont.M.df,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA,
         color = colorSpacer(startcolor = "white", middlecolors = "goldenrod", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Second round of clustering/clusters/Differential pathways top 10.pdf")

pheatmap(top30.log10.ont.M.df,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA,
         color = colorSpacer(startcolor = "white", middlecolors = "goldenrod", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Second round of clustering/clusters/Differential pathways top 30.pdf")

# Row z-scaled-log10 padj values
pheatmap(log10.ont.M.df,
         show_rownames = F,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 1, 
         border_color = NA, 
         scale = "row",
         color = colorSpacer(startcolor = "dodgerblue", middlecolors = "white", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Second round of clustering/clusters/Differential pathways row_z.pdf")

pheatmap(log10.ont.M.df,
         kmeans_k = 5,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA, 
         scale = "row",
         color = colorSpacer(startcolor = "dodgerblue", middlecolors = "white", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Second round of clustering/clusters/Differential pathways k means row_z.pdf")

pheatmap(top10.log10.ont.M.df,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA,
         scale = "row",
         color = colorSpacer(startcolor = "dodgerblue", middlecolors = "white", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Second round of clustering/clusters/Differential pathways top 10 row_z.pdf")

pheatmap(top30.log10.ont.M.df,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA,   
         scale = "row",
         color = colorSpacer(startcolor = "dodgerblue", middlecolors = "white", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Second round of clustering/clusters/Differential pathways top 30 row_z.pdf")




##==================================================================================
##==================================================================================
## Third round of clustering based on info gleaned from first round!
##==================================================================================
##==================================================================================
dir.create("myeloid_43p_celseq_integration/Third round of clustering/clusters", showWarnings = F, recursive = T)

## Visualize
# Plot Clusters
customUMAP(object     = integrated.mye.seurat, reduction = "umap", shuffle = T, title = "Initial clustering UMAP", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Third round of clustering/Integration original 43p cluster UMAP.pdf", 
           plot.width = 15)

# Refine manually
plot         <- DimPlot(integrated.mye.seurat, reduction = "umap")
select.cells <- CellSelector(plot = plot)
Idents(integrated.mye.seurat, cells = select.cells) <- "6"

plot         <- DimPlot(integrated.mye.seurat, reduction = "umap")
select.cells <- CellSelector(plot = plot)
Idents(integrated.mye.seurat, cells = select.cells) <- "4"

plot         <- DimPlot(integrated.mye.seurat, reduction = "umap")
select.cells <- CellSelector(plot = plot)
Idents(integrated.mye.seurat, cells = select.cells) <- "7"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", shuffle = T, title = "Initial clustering UMAP", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Third round of clustering/Integration cluster UMAP.pdf", 
           plot.width = 11)

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "Patient", shuffle = T, title = "Initial clustering UMAP", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Third round of clustering/Integration patient UMAP.pdf", 
           plot.width = 11)

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "Method", shuffle = T, title = "Initial clustering UMAP", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Third round of clustering/Integration method UMAP.pdf", 
           plot.width = 11)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "Patient", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "bottom") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "Method", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "bottom") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p3 <- DimPlot(integrated.mye.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p1 + p2 + p3
ggsave("myeloid_43p_celseq_integration/Third round of clustering/Integration UMAP.pdf", width = 30, height = 15)

## Plot the original 10X and 43p idents
customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "Source population identities", legend.pos = "right", 
           cols       = int.mac.pop.colors[c(-1,-2,-3)], 
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/Third round of clustering/Integration original cluster UMAP.pdf", 
           plot.width = 20)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.mac.pop.colors[c(-1,-2,-3)], pt.size = 3) +
  ggtitle("Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T, pt.size = 3) +
  ggtitle("Integrated clusters") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("myeloid_43p_celseq_integration/Third round of clustering/Integration original vs. integrated cluster UMAP.pdf", width = 20, height = 10)

# Highlight each source cluster
dir.create("myeloid_43p_celseq_integration/Third round of clustering/Original source cluster highlights", showWarnings = F, recursive = T)
for(theIdent in unique(integrated.mye.seurat$original.refined.ident)){
  customUMAP(object          = integrated.mye.seurat, reduction = "umap", label = F, shuffle = T, title = theIdent, legend.pos = "none", 
             file.name       = paste("myeloid_43p_celseq_integration/Third round of clustering/Original source cluster highlights/", theIdent, " UMAP.pdf", sep = ""),
             sizes.highlight = 4,
             cells.highlight = row.names(ref.ident)[which(ref.ident$refined.ident == theIdent)])
}

## Color only one source at a time
# 43p colored
int.43p.pop.colors <- int.mac.pop.colors
int.43p.pop.colors[grep("^CD68", names(int.43p.pop.colors), invert = T)] <- "grey"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "CEL-seq Source population identities", legend.pos = "right", 
           cols       = int.43p.pop.colors[c(-1,-2,-3)],
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/Third round of clustering/Integration original 43p cluster UMAP.pdf", 
           plot.width = 20)

# 10X colored
int.10X.pop.colors <- int.mac.pop.colors
int.10X.pop.colors[grep("^CD68", names(int.10X.pop.colors), invert = F)] <- "grey"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "10X Source population identities", legend.pos = "right", 
           cols       = int.10X.pop.colors[c(-1,-2,-3)],
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/Third round of clustering/Integration original 10X cluster UMAP.pdf", 
           plot.width = 20)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.43p.pop.colors[c(-1,-2,-3)], pt.size = 5) +
  ggtitle("CEL-seq Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.10X.pop.colors[c(-1,-2,-3)], pt.size = 5) +
  ggtitle("10X Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("myeloid_43p_celseq_integration/Third round of clustering/Integration 43p vs. 10X clusters UMAP.pdf", width = 20, height = 10)

## Save the integrated seurat object
saveRDS(integrated.mye.seurat, "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")
#integrated.mye.seurat <- readRDS(file = "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")


#=========================================================================================================================
## Define the clusters: plot markers, ontologies, and gene of interest
##========================================================================================================================
## Scale all genes from RNA assay for diff gene calling and visualisation
# Set default to RNA
DefaultAssay(integrated.mye.seurat) <- "RNA"

# Normalize
integrated.mye.seurat <- NormalizeData(integrated.mye.seurat, normalization.method = "LogNormalize")

# Scale
integrated.mye.seurat <- ScaleData(integrated.mye.seurat, features = row.names(integrated.mye.seurat))

##========================================================================================================================
## Define marker genes per cluster
integrated.mye.seurat.markers <- FindAllMarkers(object = integrated.mye.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save the top9 markers per cluster
sep.markers <- integrated.mye.seurat.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
num.clus <- length(unique(Idents(integrated.mye.seurat)))
for(i in levels(Idents(integrated.mye.seurat))){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = integrated.mye.seurat,
                     name      = paste("myeloid_43p_celseq_integration/Third round of clustering/clusters/Cluster ",i, " markers", sep = ""),
                     
                     Vln.width = 25, Vln.height = 25)
}

# Save markers per cluster to disk
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("myeloid_43p_celseq_integration/Third round of clustering/clusters/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

##========================================================================================================================
# Look at some genes
bunchOfCustomPlots(features   = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                  "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                  "CD79A", "ACTA2", "CD34", "KIT", "OLR1" ), 
                   object     = integrated.mye.seurat,
                   ncol       = 5,
                   name       = "myeloid_43p_celseq_integration/Third round of clustering/Cell type defining genes",
                   Vln.width  = 25, 
                   Vln.height = 15
)

# SMC Monocyte genes
bunchOfCustomPlots(features = c("ACTA2", "MYH11", "CD14", "CD68", "KLF4", "CD4", "TAGLN", "ADGRE1", "LGALS3"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Third round of clustering/SMC mono genes",
                   Vln.width = 25, Vln.height = 15
)

#Extra mono cell markers
bunchOfCustomPlots(features = c("FCGR3A", "CCR2", "CCR5", "SELL", "ITGAX", "CX3CR1"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Third round of clustering/Extra mono marker genes",
)

# Myeloid Markers
bunchOfCustomPlots(features = c("CD200R1", "TGM2", "CD163"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Third round of clustering/M2 genes")

bunchOfCustomPlots(features = c("CD2","ITGAM", "PECAM1", "NCAM1", "SELL", "CSF1R"), 
                   object    = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Third round of clustering/Monocyte genes",
)

bunchOfCustomPlots(features = c("CD68","CD14","ITGAM","CD1C", "CLEC10A", "FCER1A", "CD2", "CD3E", "CD4"),
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Third round of clustering/Myeloid genes", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("CASP1","IL1B","TLR4","ABCA1","OLR1", "TREM2","CASP4","KLF4","TNF","IL18","ABCG1", "CD9"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Third round of clustering/Foam vs. pro inf genes", 
                   ncol     = 3,
                   Vln.width = 30, Vln.height = 15
)

bunchOfCustomPlots(features = c("CD14","CD68","CD86","CD1C", "FCGR2B", "CLEC10A", "FCER1A", "HLA-DQA2", "HLA-DQB2"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Third round of clustering/DC genes", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("NLRP3", "IL1B", "CASP1", "CASP4", "IL18", "PYCARD"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Third round of clustering/Inflammasome genes", 
                   ncol     = 3,
)

bunchOfCustomPlots(features = c("CD14", "FCGR3A", "FCGR1A", "CCR2", "PADI4", "SECISBP2L", "CD86", "ITGAX", "LYZ", "CD33", "CD36", "ITGAM"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Third round of clustering/Mono subtype genes CD14 CD16 CD64 CCR2 PADI4 SLAN CD86 CD11c LYZ CD33 CD36 CD11b", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("NR4A1"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Third round of clustering/NR4A1",
                   ncol     = 1,
)

bunchOfCustomPlots(features  = c("CD14","FCGR3A", "CD68", 'CD86', "SELL", "CSF1R", "ITGAM", "ITGAX", "PECAM1"), 
                   object    = integrated.mye.seurat,
                   ncol      = 3, Vln.stack = T,
                   name      = "myeloid_43p_celseq_integration/Third round of clustering/macro and mono genes",
                   Vln.width = 20, Vln.height = 20
)

bunchOfCustomPlots(features  = c("TREM2", "CD9", "GPNMB", "LYVE1", "MRC1", "FOLR2" ,"OLR1", "ABCA1", "ABCG1", "IL1B", "TNF", "CASP1", "CD14", "FCGR1A", "FCGR3A"), 
                   object    = integrated.mye.seurat,
                   ncol      = 3,
                   name      = "myeloid_43p_celseq_integration/Third round of clustering/LAM - resident - foam - inflammatory - general markers",
                   Vln.width = 20, Vln.height = 20
)


bunchOfCustomPlots(features = c("TIMD4", "CCR2", "LYVE1", "FOLR2"), 
                   object   = integrated.mye.seurat, assay = "RNA",
                   name     = "myeloid_43p_celseq_integration/Third round of clustering/core_Gene_macs",
                   ncol     = 4, 
                   Vln.width = 25
)

bunchOfCustomPlots(features = c("CX3CR1", "LYVE1", "MRC1", "F13A1", "GAS6", "FOLR2", "TREM2","CD9", "IL1B"), 
                   object   = integrated.mye.seurat, assay = "RNA",
                   name     = "myeloid_43p_celseq_integration/Third round of clustering/resident_macs",
                   ncol     = 3, 
                   Vln.width = 25
)

customFeature(object = integrated.mye.seurat, features = c("MKI67"), name = "myeloid_43p_celseq_integration/Third round of clustering/MKI67.pdf")


##==================================================================================
## Get ontologies
# All markers
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  get_ontology(x, name = i, outdir = "myeloid_43p_celseq_integration/Third round of clustering/clusters", universe = integrated.mye.seurat, full_GSEA = F, volcano.plot = T) 
}

# Top 8
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  get_ontology(x, name = paste("Top 8", i, sep = " "), outdir = "myeloid_43p_celseq_integration/Third round of clustering/clusters", universe = integrated.mye.seurat, full_GSEA = F, volcano.plot = F, plot.top.n = 8) 
}


##==================================================================================
## Marker gene pathways heatmap
ont.M <- list()
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  ont.M[[i]] <- get_ontology(res = x, name = i, outdir = "myeloid_43p_celseq_integration/Third round of clustering/clusters", universe = integrated.mye.seurat, full_GSEA = F, return.data = T, plot.top.n = 10)
}

# Plot differential pathways from top 30 per population
# Collect all data into a neat 'frame
ont.M.df <- merge(x        = ont.M[[1]]$Pathways, 
                  y        = ont.M[[2]]$Pathways, 
                  by       = "name", 
                  all      = T, 
                  suffixes = paste(".", names(ont.M)[1:2], sep = ""))

for(i in names(ont.M[3:length(names(ont.M))])){
  ont.M.df <- merge(x        = ont.M.df, 
                    y        = ont.M[[i]]$Pathways, 
                    by       = "name", 
                    all      = T)
  
  colnames(ont.M.df)[(length(colnames(ont.M.df))-2):length(colnames(ont.M.df))] <- paste(colnames(ont.M.df)[(length(colnames(ont.M.df))-2):length(colnames(ont.M.df))], i, sep = ".")
}

# Clean up
ont.M.df[is.na(ont.M.df)] <- 1
row.names(ont.M.df)       <- ont.M.df$name
ont.M.df$name             <- NULL

# Keep only padj columns
ont.M.df <- ont.M.df[,seq(from = 1, to = length(colnames(ont.M.df)), by = 3)]
head(ont.M.df)            # Neat B)

# Log transform the data
log10.ont.M.df                        <- -log10(ont.M.df)
log10.ont.M.df[log10.ont.M.df == Inf] <- 0
head(log10.ont.M.df)

# Keep top 10 rows per population
top10.log10.ont.M.df <- log10.ont.M.df[order(log10.ont.M.df[,1], decreasing = T),][1:10,]
for(i in length(colnames(ont.M.df))){
  top10.log10.ont.M.df <- rbind(top10.log10.ont.M.df, log10.ont.M.df[order(log10.ont.M.df[,i], decreasing = T),][1:10,])
}
top10.log10.ont.M.df <- top10.log10.ont.M.df[grep("1$", row.names(top10.log10.ont.M.df), invert = T), ]

# Keep top 30 rows per population
top30.log10.ont.M.df <- log10.ont.M.df[order(log10.ont.M.df[,1], decreasing = T),][1:30,]
for(i in length(colnames(ont.M.df))){
  top30.log10.ont.M.df <- rbind(top30.log10.ont.M.df, log10.ont.M.df[order(log10.ont.M.df[,i], decreasing = T),][1:30,])
}
top30.log10.ont.M.df <- top30.log10.ont.M.df[grep("1$", row.names(top30.log10.ont.M.df), invert = T), ]

## And plot
# Absolute -log10 padj values
pheatmap(log10.ont.M.df,
         show_rownames = F,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 1, 
         border_color = NA, 
         color = colorSpacer(startcolor = "white", middlecolors = "goldenrod", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Third round of clustering/clusters/Differential pathways.pdf")

pheatmap(log10.ont.M.df,
         kmeans_k = 5,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA, 
         color = colorSpacer(startcolor = "white", middlecolors = "goldenrod", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Third round of clustering/clusters/Differential pathways k means.pdf")

pheatmap(top10.log10.ont.M.df,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA,
         color = colorSpacer(startcolor = "white", middlecolors = "goldenrod", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Third round of clustering/clusters/Differential pathways top 10.pdf")

pheatmap(top30.log10.ont.M.df,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA,
         color = colorSpacer(startcolor = "white", middlecolors = "goldenrod", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Third round of clustering/clusters/Differential pathways top 30.pdf")

# Row z-scaled-log10 padj values
pheatmap(log10.ont.M.df,
         show_rownames = F,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 1, 
         border_color = NA, 
         scale = "row",
         color = colorSpacer(startcolor = "dodgerblue", middlecolors = "white", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Third round of clustering/clusters/Differential pathways row_z.pdf")

pheatmap(log10.ont.M.df,
         kmeans_k = 5,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA, 
         scale = "row",
         color = colorSpacer(startcolor = "dodgerblue", middlecolors = "white", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Third round of clustering/clusters/Differential pathways k means row_z.pdf")

pheatmap(top10.log10.ont.M.df,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA,
         scale = "row",
         color = colorSpacer(startcolor = "dodgerblue", middlecolors = "white", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Third round of clustering/clusters/Differential pathways top 10 row_z.pdf")

pheatmap(top30.log10.ont.M.df,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA,   
         scale = "row",
         color = colorSpacer(startcolor = "dodgerblue", middlecolors = "white", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Third round of clustering/clusters/Differential pathways top 30 row_z.pdf")



##==================================================================================
##==================================================================================
## Fourth round of clustering based on info gleaned from first round!
##==================================================================================
##==================================================================================
dir.create("myeloid_43p_celseq_integration/Fourth round of clustering/clusters", showWarnings = F, recursive = T)

## Visualize
# Plot Clusters
customUMAP(object     = integrated.mye.seurat, reduction = "umap", shuffle = T, title = "Initial clustering UMAP", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Fourth round of clustering/Integration original 43p cluster UMAP.pdf", 
           plot.width = 15)

# Save idents for later just in case
integrated.mye.seurat[["method.int.idents.pre.merge"]] <- Idents(integrated.mye.seurat)
  
# Merge 2, 4, and 6
select.cells <- WhichCells(integrated.mye.seurat, idents = c("4", "6"))
Idents(integrated.mye.seurat, cells = select.cells) <- "2"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", shuffle = T, title = "Machine idents", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Fourth round of clustering/Integration cluster UMAP.pdf", 
           plot.width = 11)

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "Patient", shuffle = T, title = "Patient", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Fourth round of clustering/Integration patient UMAP.pdf", 
           plot.width = 11)

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "Method", shuffle = T, title = "IMethod", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Fourth round of clustering/Integration method UMAP.pdf", 
           plot.width = 11)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "Patient", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "bottom") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "Method", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "bottom") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p3 <- DimPlot(integrated.mye.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p1 + p2 + p3
ggsave("myeloid_43p_celseq_integration/Fourth round of clustering/Integration UMAP.pdf", width = 30, height = 15)

## Plot the original 10X and 43p idents
customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "Source population identities", legend.pos = "right", 
           cols       = int.mac.pop.colors[c(-1,-2,-3)], 
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/Fourth round of clustering/Integration original cluster UMAP.pdf", 
           plot.width = 20)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.mac.pop.colors[c(-1,-2,-3)], pt.size = 3) +
  ggtitle("Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T, pt.size = 3) +
  ggtitle("Integrated clusters") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("myeloid_43p_celseq_integration/Fourth round of clustering/Integration original vs. integrated cluster UMAP.pdf", width = 20, height = 10)

# Highlight each source cluster
dir.create("myeloid_43p_celseq_integration/Fourth round of clustering/Original source cluster highlights", showWarnings = F, recursive = T)
for(theIdent in unique(integrated.mye.seurat$original.refined.ident)){
  customUMAP(object          = integrated.mye.seurat, reduction = "umap", label = F, shuffle = T, title = theIdent, legend.pos = "none", 
             file.name       = paste("myeloid_43p_celseq_integration/Fourth round of clustering/Original source cluster highlights/", theIdent, " UMAP.pdf", sep = ""),
             sizes.highlight = 4,
             cells.highlight = row.names(ref.ident)[which(ref.ident$refined.ident == theIdent)])
}

## Color only one source at a time
# 43p colored
int.43p.pop.colors <- int.mac.pop.colors
int.43p.pop.colors[grep("^CD68", names(int.43p.pop.colors), invert = T)] <- "grey"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "CEL-seq Source population identities", legend.pos = "right", 
           cols       = int.43p.pop.colors[c(-1,-2,-3)],
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/Fourth round of clustering/Integration original 43p cluster UMAP.pdf", 
           plot.width = 20)

# 10X colored
int.10X.pop.colors <- int.mac.pop.colors
int.10X.pop.colors[grep("^CD68", names(int.10X.pop.colors), invert = F)] <- "grey"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "10X Source population identities", legend.pos = "right", 
           cols       = int.10X.pop.colors[c(-1,-2,-3)],
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/Fourth round of clustering/Integration original 10X cluster UMAP.pdf", 
           plot.width = 20)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.43p.pop.colors[c(-1,-2,-3)], pt.size = 5) +
  ggtitle("CEL-seq Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.10X.pop.colors[c(-1,-2,-3)], pt.size = 5) +
  ggtitle("10X Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("myeloid_43p_celseq_integration/Fourth round of clustering/Integration 43p vs. 10X clusters UMAP.pdf", width = 20, height = 10)

## Save the integrated seurat object
saveRDS(integrated.mye.seurat, "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")
#integrated.mye.seurat <- readRDS(file = "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")


#=========================================================================================================================
## Define the clusters: plot markers, ontologies, and gene of interest
##========================================================================================================================
## Scale all genes from RNA assay for diff gene calling and visualisation
# Set default to RNA
DefaultAssay(integrated.mye.seurat) <- "RNA"

# Normalize
integrated.mye.seurat <- NormalizeData(integrated.mye.seurat, normalization.method = "LogNormalize")

# Scale
integrated.mye.seurat <- ScaleData(integrated.mye.seurat, features = row.names(integrated.mye.seurat))

##========================================================================================================================
## Define marker genes per cluster
integrated.mye.seurat.markers <- FindAllMarkers(object = integrated.mye.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save the top9 markers per cluster
sep.markers <- integrated.mye.seurat.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
num.clus <- length(unique(Idents(integrated.mye.seurat)))
for(i in levels(Idents(integrated.mye.seurat))){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = integrated.mye.seurat,
                     name      = paste("myeloid_43p_celseq_integration/Fourth round of clustering/clusters/Cluster ",i, " markers", sep = ""),
                     
                     Vln.width = 25, Vln.height = 25)
}

# Save markers per cluster to disk
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("myeloid_43p_celseq_integration/Fourth round of clustering/clusters/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

##========================================================================================================================
# Look at some genes
bunchOfCustomPlots(features   = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                  "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                  "CD79A", "ACTA2", "CD34", "KIT", "OLR1" ), 
                   object     = integrated.mye.seurat,
                   ncol       = 5,
                   name       = "myeloid_43p_celseq_integration/Fourth round of clustering/Cell type defining genes",
                   Vln.width  = 25, 
                   Vln.height = 15
)

# SMC Monocyte genes
bunchOfCustomPlots(features = c("ACTA2", "MYH11", "CD14", "CD68", "KLF4", "CD4", "TAGLN", "ADGRE1", "LGALS3"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Fourth round of clustering/SMC mono genes",
                   Vln.width = 25, Vln.height = 15
)

#Extra mono cell markers
bunchOfCustomPlots(features = c("FCGR3A", "CCR2", "CCR5", "SELL", "ITGAX", "CX3CR1"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Fourth round of clustering/Extra mono marker genes",
)

# Myeloid Markers
bunchOfCustomPlots(features = c("CD200R1", "TGM2", "CD163"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Fourth round of clustering/M2 genes")

bunchOfCustomPlots(features = c("CD2","ITGAM", "PECAM1", "NCAM1", "SELL", "CSF1R"), 
                   object    = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Fourth round of clustering/Monocyte genes",
)

bunchOfCustomPlots(features = c("CD68","CD14","ITGAM","CD1C", "CLEC10A", "FCER1A", "CD2", "CD3E", "CD4"),
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Fourth round of clustering/Myeloid genes", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("CASP1","IL1B","TLR4","ABCA1","OLR1", "TREM2","CASP4","KLF4","TNF","IL18","ABCG1", "CD9"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Fourth round of clustering/Foam vs. pro inf genes", 
                   ncol     = 3,
                   Vln.width = 30, Vln.height = 15
)

bunchOfCustomPlots(features = c("CD14","CD68","CD86","CD1C", "FCGR2B", "CLEC10A", "FCER1A", "HLA-DQA2", "HLA-DQB2"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Fourth round of clustering/DC genes", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("NLRP3", "IL1B", "CASP1", "CASP4", "IL18", "PYCARD"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Fourth round of clustering/Inflammasome genes", 
                   ncol     = 3,
)

bunchOfCustomPlots(features = c("CD14", "FCGR3A", "FCGR1A", "CCR2", "PADI4", "SECISBP2L", "CD86", "ITGAX", "LYZ", "CD33", "CD36", "ITGAM"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Fourth round of clustering/Mono subtype genes CD14 CD16 CD64 CCR2 PADI4 SLAN CD86 CD11c LYZ CD33 CD36 CD11b", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("NR4A1"), 
                   object   = integrated.mye.seurat,
                   name     = "myeloid_43p_celseq_integration/Fourth round of clustering/NR4A1",
                   ncol     = 1,
)

bunchOfCustomPlots(features  = c("CD14","FCGR3A", "CD68", 'CD86', "SELL", "CSF1R", "ITGAM", "ITGAX", "PECAM1"), 
                   object    = integrated.mye.seurat,
                   ncol      = 3, Vln.stack = T,
                   name      = "myeloid_43p_celseq_integration/Fourth round of clustering/macro and mono genes",
                   Vln.width = 20, Vln.height = 20
)

bunchOfCustomPlots(features  = c("TREM2", "CD9", "GPNMB", "LYVE1", "MRC1", "FOLR2" ,"OLR1", "ABCA1", "ABCG1", "IL1B", "TNF", "CASP1", "CD14", "FCGR1A", "FCGR3A"), 
                   object    = integrated.mye.seurat,
                   ncol      = 3,
                   name      = "myeloid_43p_celseq_integration/Fourth round of clustering/LAM - resident - foam - inflammatory - general markers",
                   Vln.width = 20, Vln.height = 20
)


bunchOfCustomPlots(features = c("TIMD4", "CCR2", "LYVE1", "FOLR2"), 
                   object   = integrated.mye.seurat, assay = "RNA",
                   name     = "myeloid_43p_celseq_integration/Fourth round of clustering/core_Gene_macs",
                   ncol     = 4, 
                   Vln.width = 25
)

bunchOfCustomPlots(features = c("CX3CR1", "LYVE1", "MRC1", "F13A1", "GAS6", "FOLR2", "TREM2","CD9", "IL1B"), 
                   object   = integrated.mye.seurat, assay = "RNA",
                   name     = "myeloid_43p_celseq_integration/Fourth round of clustering/resident_macs",
                   ncol     = 3, 
                   Vln.width = 25
)

customFeature(object = integrated.mye.seurat, features = c("MKI67"), name = "myeloid_43p_celseq_integration/Fourth round of clustering/MKI67.pdf")


##==================================================================================
## Get ontologies
# All markers
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  get_ontology(x, name = i, outdir = "myeloid_43p_celseq_integration/Fourth round of clustering/clusters", universe = integrated.mye.seurat, full_GSEA = F, volcano.plot = T) 
}

# Top 8
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  get_ontology(x, name = paste("Top 8", i, sep = " "), outdir = "myeloid_43p_celseq_integration/Fourth round of clustering/clusters", universe = integrated.mye.seurat, full_GSEA = F, volcano.plot = F, plot.top.n = 8) 
}


##==================================================================================
## Marker gene pathways heatmap
ont.M <- list()
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  ont.M[[i]] <- get_ontology(res = x, name = i, outdir = "myeloid_43p_celseq_integration/Fourth round of clustering/clusters", universe = integrated.mye.seurat, full_GSEA = F, return.data = T, plot.top.n = 10)
}

# Plot differential pathways from top 30 per population
# Collect all data into a neat 'frame
ont.M.df <- merge(x        = ont.M[[1]]$Pathways, 
                  y        = ont.M[[2]]$Pathways, 
                  by       = "name", 
                  all      = T, 
                  suffixes = paste(".", names(ont.M)[1:2], sep = ""))

for(i in names(ont.M[3:length(names(ont.M))])){
  ont.M.df <- merge(x        = ont.M.df, 
                    y        = ont.M[[i]]$Pathways, 
                    by       = "name", 
                    all      = T)
  
  colnames(ont.M.df)[(length(colnames(ont.M.df))-2):length(colnames(ont.M.df))] <- paste(colnames(ont.M.df)[(length(colnames(ont.M.df))-2):length(colnames(ont.M.df))], i, sep = ".")
}

# Clean up
ont.M.df[is.na(ont.M.df)] <- 1
row.names(ont.M.df)       <- ont.M.df$name
ont.M.df$name             <- NULL

# Keep only padj columns
ont.M.df <- ont.M.df[,seq(from = 1, to = length(colnames(ont.M.df)), by = 3)]
head(ont.M.df)            # Neat B)

# Log transform the data
log10.ont.M.df                        <- -log10(ont.M.df)
log10.ont.M.df[log10.ont.M.df == Inf] <- 0
head(log10.ont.M.df)

# Keep top 10 rows per population
top10.log10.ont.M.df <- log10.ont.M.df[order(log10.ont.M.df[,1], decreasing = T),][1:10,]
for(i in length(colnames(ont.M.df))){
  top10.log10.ont.M.df <- rbind(top10.log10.ont.M.df, log10.ont.M.df[order(log10.ont.M.df[,i], decreasing = T),][1:10,])
}
top10.log10.ont.M.df <- top10.log10.ont.M.df[grep("1$", row.names(top10.log10.ont.M.df), invert = T), ]

# Keep top 30 rows per population
top30.log10.ont.M.df <- log10.ont.M.df[order(log10.ont.M.df[,1], decreasing = T),][1:30,]
for(i in length(colnames(ont.M.df))){
  top30.log10.ont.M.df <- rbind(top30.log10.ont.M.df, log10.ont.M.df[order(log10.ont.M.df[,i], decreasing = T),][1:30,])
}
top30.log10.ont.M.df <- top30.log10.ont.M.df[grep("1$", row.names(top30.log10.ont.M.df), invert = T), ]

## And plot
# Absolute -log10 padj values
pheatmap(log10.ont.M.df,
         show_rownames = F,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 1, 
         border_color = NA, 
         color = colorSpacer(startcolor = "white", middlecolors = "goldenrod", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Fourth round of clustering/clusters/Differential pathways.pdf")

pheatmap(log10.ont.M.df,
         kmeans_k = 5,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA, 
         color = colorSpacer(startcolor = "white", middlecolors = "goldenrod", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Fourth round of clustering/clusters/Differential pathways k means.pdf")

pheatmap(top10.log10.ont.M.df,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA,
         color = colorSpacer(startcolor = "white", middlecolors = "goldenrod", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Fourth round of clustering/clusters/Differential pathways top 10.pdf")

pheatmap(top30.log10.ont.M.df,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA,
         color = colorSpacer(startcolor = "white", middlecolors = "goldenrod", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Fourth round of clustering/clusters/Differential pathways top 30.pdf")

# Row z-scaled-log10 padj values
pheatmap(log10.ont.M.df,
         show_rownames = F,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 1, 
         border_color = NA, 
         scale = "row",
         color = colorSpacer(startcolor = "dodgerblue", middlecolors = "white", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Fourth round of clustering/clusters/Differential pathways row_z.pdf")

pheatmap(log10.ont.M.df,
         kmeans_k = 5,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA, 
         scale = "row",
         color = colorSpacer(startcolor = "dodgerblue", middlecolors = "white", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Fourth round of clustering/clusters/Differential pathways k means row_z.pdf")

pheatmap(top10.log10.ont.M.df,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA,
         scale = "row",
         color = colorSpacer(startcolor = "dodgerblue", middlecolors = "white", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Fourth round of clustering/clusters/Differential pathways top 10 row_z.pdf")

pheatmap(top30.log10.ont.M.df,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = T,
         cellwidth = 30, 
         cellheight = 15, 
         border_color = NA,   
         scale = "row",
         color = colorSpacer(startcolor = "dodgerblue", middlecolors = "white", endcolor = "red3", steps = 255, return.colors = T),
         filename = "myeloid_43p_celseq_integration/Fourth round of clustering/clusters/Differential pathways top 30 row_z.pdf")


##==================================================================================
##==================================================================================
## Rename the final clusters
integrated.mye.seurat[["method_integrated.numbered.idents"]] <- Idents(object = integrated.mye.seurat)
integrated.mye.seurat <- RenameIdents(integrated.mye.seurat, 
                                          "0"   = "CD14+TREM2+FOLR2-ABCG+ Lipid Associated Macrophages",
                                          "1"   = "CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages",
                                          "2"   = "CD14+IL1B+SELL+S100A8+ Antigen-presenting Inflammatory Monocyte-derived Macrophages",
                                          "2.1" = "CD14+IL1B+SELL+MX1+ Interferon Activated Inflammatory Monocyte-derived Macrophages",
                                          "3"   = "CD14+TNF+TREM2+FOLR2+ Inflammatory Resident-like Lipid Associated Macrophages",
                                          "5"   = "CD14+-IL1B+SELL+CD16+ Migrating Inflammatory Monocyte-derived Macrophages",
                                          "7"   = "CD14+TREM2-OLR1+ABCA+ Foamy Macrophages",
                                          "7.1" = "CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages",
                                          "8"   = "CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages"
)


## Define useful pop colors
M.int_refined.pops <- levels(Idents(integrated.mye.seurat))

# Inflammatory Mo-Macs
M.int_refined.pop.colors        <- colorSpacer(startcolor = "firebrick4", endcolor = "goldenrod3", steps = 3, return.colors = T)
names(M.int_refined.pop.colors) <- M.int_refined.pops[grep(" Mono", M.int_refined.pops)]

# Foamy Cells
tmp.names                       <- names(M.int_refined.pop.colors)
M.int_refined.pop.colors        <- c(M.int_refined.pop.colors, colorSpacer(startcolor = "plum1", endcolor = "darkslateblue", steps = 3, return.colors = T))
names(M.int_refined.pop.colors) <- c(tmp.names, M.int_refined.pops[grep("Foamy", M.int_refined.pops)])

# Res Macs
tmp.names                       <- names(M.int_refined.pop.colors)
M.int_refined.pop.colors        <- c(M.int_refined.pop.colors, colorSpacer(startcolor = "chartreuse1", endcolor = "chartreuse4", steps = 3, return.colors = T))
names(M.int_refined.pop.colors) <- c(tmp.names, M.int_refined.pops[grep("Resident|ABCG", M.int_refined.pops)])

length(M.int_refined.pop.colors)

# Add mono colors
M.int_refined.pop.colors <- c(M_refined.pop.colors[1:3], M.int_refined.pop.colors)

# Plot a shiny new UMAP
customUMAP(object     = integrated.mye.seurat, label = TRUE, label.size = 8, shuffle = T, 
           cols       = M.int_refined.pop.colors, 
           legend.pos = "none",
           file.name  = "myeloid_43p_celseq_integration/Fourth round of clustering/UMAP clusters names labels.pdf")

customUMAP(object     = integrated.mye.seurat, label = FALSE, shuffle = T, plot.width = 18,
           cols       = M.int_refined.pop.colors, 
           legend.pos = "right",
           file.name  = "myeloid_43p_celseq_integration/Fourth round of clustering/UMAP clusters names legend.pdf")

## Save the integrated seurat object
saveRDS(integrated.mye.seurat, "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")
#integrated.mye.seurat <- readRDS(file = "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")




##==================================================================================
##==================================================================================
## Fifth Round of clustering
## Based on back-mapping to the 10X data, clus 6 actually clusters neatly there too, 
## and 4 makes up the bulk of the clas-mon infiltration. So:
## Back up for  a second and investigate cluster 2-4-6 diffs. 
##==================================================================================
##==================================================================================
dir.create("myeloid_43p_celseq_integration/Fifth round of clustering/clusters", recursive = T, showWarnings = F)

clus246.diff.markers <- list()
clus246.diff.markers[["2_v_4"]] <- FindMarkers(integrated.mye.seurat, group.by = "method.int.idents.pre.merge", ident.1 = 2, ident.2 = 4, logfc.threshold = 0.25, min.diff.pct = 0.25)
clus246.diff.markers[["2_v_6"]] <- FindMarkers(integrated.mye.seurat, group.by = "method.int.idents.pre.merge", ident.1 = 2, ident.2 = 6, logfc.threshold = 0.25, min.diff.pct = 0.25)
clus246.diff.markers[["4_v_6"]] <- FindMarkers(integrated.mye.seurat, group.by = "method.int.idents.pre.merge", ident.1 = 4, ident.2 = 6, logfc.threshold = 0.25, min.diff.pct = 0.25)

# Save the top9 markers per comparison
sep.markers <- list()
sep.markers[["2_v_4"]] <- subset(clus246.diff.markers[["2_v_4"]], avg_log2FC > 0)
sep.markers[["4_v_2"]] <- subset(clus246.diff.markers[["2_v_4"]], avg_log2FC < 0)

sep.markers[["2_v_6"]] <- subset(clus246.diff.markers[["2_v_6"]], avg_log2FC > 0)
sep.markers[["6_v_2"]] <- subset(clus246.diff.markers[["2_v_6"]], avg_log2FC < 0)

sep.markers[["4_v_6"]] <- subset(clus246.diff.markers[["4_v_6"]], avg_log2FC > 0)
sep.markers[["6_v_4"]] <- subset(clus246.diff.markers[["4_v_6"]], avg_log2FC < 0)


# plot the top9 markers per comparison
for(i in names(sep.markers)){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = head(row.names(sep.markers[[i]]), n = 9), 
                     object    = integrated.mye.seurat,
                     group.by  = "method.int.idents.pre.merge",
                     name      = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Cluster ",i, " markers", sep = ""),
                     Vln.width = 25, Vln.height = 25)
}

# Save markers per comparison to disk
for(i in names(sep.markers)){
  x      <- sep.markers[[i]]
  x$gene <- row.names(x)
  x      <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

# Get pathways per comparison
for(i in names(sep.markers)){
  x      <- sep.markers[[i]]
  x$gene <- row.names(x)
  get_ontology(res = x, name = i, outdir = "myeloid_43p_celseq_integration/Fifth round of clustering/clusters", universe = integrated.mye.seurat, full_GSEA = F, return.data = T, plot.top.n = 10)
}


##===============================================================================================================
# 4 and 6 don;t really differ meangingfully, let's try to merge those and keep 6 apart.
integrated.mye.seurat[["method_integrated.overmerged.idents"]] <- Idents(object = integrated.mye.seurat)
Idents(integrated.mye.seurat) <- integrated.mye.seurat$method.int.idents.pre.merge

integrated.mye.seurat <- RenameIdents(integrated.mye.seurat, 
                                      "0"   = "CD14+TREM2+FOLR2-ABCG+ Lipid Associated Macrophages",
                                      "1"   = "CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages",
                                      "2"   = "CD14+IL1B+SELL+S100A8+ Antigen-presenting Inflammatory Monocyte-derived Macrophages",
                                      "4"   = "CD14+IL1B+SELL+S100A8+ Antigen-presenting Inflammatory Monocyte-derived Macrophages",
                                      "6"   = "CD14+IL1B+SELL+CD16- Antigen-presenting Inflammatory Monocyte-derived Macrophages",
                                      "2.1" = "CD14+IL1B+SELL+MX1+ Interferon Activated Inflammatory Monocyte-derived Macrophages",
                                      "3"   = "CD14+TNF+TREM2+FOLR2+ Inflammatory Resident-like Lipid Associated Macrophages",
                                      "5"   = "CD14+-IL1B+SELL+CD16+ Migrating Inflammatory Monocyte-derived Macrophages",
                                      "7"   = "CD14+TREM2-OLR1+ABCA+ Foamy Macrophages",
                                      "7.1" = "CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages",
                                      "8"   = "CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages"
)


## Define useful pop colors
M.int_refined.pops <- levels(Idents(integrated.mye.seurat))

# Inflammatory Mo-Macs
M.int_refined.pop.colors        <- colorSpacer(startcolor = "goldenrod3", endcolor = "firebrick4", steps = 4, return.colors = T)
names(M.int_refined.pop.colors) <- M.int_refined.pops[grep(" Mono", M.int_refined.pops)]

# Foamy Cells
tmp.names                       <- names(M.int_refined.pop.colors)
M.int_refined.pop.colors        <- c(M.int_refined.pop.colors, colorSpacer(startcolor = "darkslateblue", endcolor = "plum1", steps = 3, return.colors = T))
names(M.int_refined.pop.colors) <- c(tmp.names, M.int_refined.pops[grep("Foamy", M.int_refined.pops)])

# Res Macs
tmp.names                       <- names(M.int_refined.pop.colors)
M.int_refined.pop.colors        <- c(M.int_refined.pop.colors, colorSpacer(startcolor = "chartreuse4", endcolor = "chartreuse1", steps = 3, return.colors = T))
names(M.int_refined.pop.colors) <- c(tmp.names, M.int_refined.pops[grep("Resident|ABCG", M.int_refined.pops)])

length(M.int_refined.pop.colors)

# Add mono colors
M.int_refined.pop.colors <- c(M_refined.pop.colors[1:3], M.int_refined.pop.colors)


## Add updated numbered clustered metadata for violins
tmp.idents <- integrated.mye.seurat$method.int.idents.pre.merge
tmp.idents[tmp.idents == 4] <- 2
integrated.mye.seurat <- AddMetaData(integrated.mye.seurat, metadata = tmp.idents, col.name = "method.int.idents.24.merge")


## Plot a shiny new UMAP
customUMAP(object     = integrated.mye.seurat, label = TRUE, label.size = 8, shuffle = T, 
           cols       = M.int_refined.pop.colors, 
           legend.pos = "none",
           file.name  = "myeloid_43p_celseq_integration/Fifth round of clustering/UMAP clusters names labels.pdf")

customUMAP(object     = integrated.mye.seurat, label = FALSE, shuffle = T, plot.width = 18,
           cols       = M.int_refined.pop.colors, 
           legend.pos = "right",
           file.name  = "myeloid_43p_celseq_integration/Fifth round of clustering/UMAP clusters names legend.pdf")

customUMAP(object     = integrated.mye.seurat, label = TRUE, label.size = 8, shuffle = T, 
           cols       = M.int_refined.pop.colors, 
           legend.pos = "none",
           file.name  = "myeloid_43p_celseq_integration/UMAP clusters names labels.pdf")

customUMAP(object     = integrated.mye.seurat, label = FALSE, shuffle = T, plot.width = 18,
           cols       = M.int_refined.pop.colors, 
           legend.pos = "right",
           file.name  = "myeloid_43p_celseq_integration/UMAP clusters names legend.pdf")

## Save the integrated seurat object
saveRDS(integrated.mye.seurat, "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")
#integrated.mye.seurat <- readRDS(file = "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")


## Plot some more UMAPs depedning on cell origin
customUMAP(object     = integrated.mye.seurat, reduction = "umap", shuffle = T, title = "Refined idents", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Fifth round of clustering/Integration cluster UMAP.pdf", 
           plot.width = 11)

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "Patient", shuffle = T, title = "Patient", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Fifth round of clustering/Integration patient UMAP.pdf", 
           plot.width = 11)

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "Method", shuffle = T, title = "IMethod", legend.pos = "right", 
           file.name  = "myeloid_43p_celseq_integration/Fifth round of clustering/Integration method UMAP.pdf", 
           plot.width = 11)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "Patient", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "bottom") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "Method", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "bottom") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p3 <- DimPlot(integrated.mye.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p1 + p2 + p3
ggsave("myeloid_43p_celseq_integration/Fifth round of clustering/Integration UMAP.pdf", width = 30, height = 15)

## Plot the original 10X and 43p idents
customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "Source population identities", legend.pos = "right", 
           cols       = int.mac.pop.colors[c(-1,-2,-3)], 
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/Fifth round of clustering/Integration original cluster UMAP.pdf", 
           plot.width = 20)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.mac.pop.colors[c(-1,-2,-3)], pt.size = 3) +
  ggtitle("Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T, pt.size = 3) +
  ggtitle("Integrated clusters") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("myeloid_43p_celseq_integration/Fifth round of clustering/Integration original vs. integrated cluster UMAP.pdf", width = 20, height = 10)

# Highlight each source cluster
dir.create("myeloid_43p_celseq_integration/Fifth round of clustering/Original source cluster highlights", showWarnings = F, recursive = T)
for(theIdent in unique(integrated.mye.seurat$original.refined.ident)){
  customUMAP(object          = integrated.mye.seurat, reduction = "umap", label = F, shuffle = T, title = theIdent, legend.pos = "none", 
             file.name       = paste("myeloid_43p_celseq_integration/Fifth round of clustering/Original source cluster highlights/", theIdent, " UMAP.pdf", sep = ""),
             sizes.highlight = 4,
             cells.highlight = row.names(ref.ident)[which(ref.ident$refined.ident == theIdent)])
}

## Color only one source at a time
# 43p colored
int.43p.pop.colors <- int.mac.pop.colors
int.43p.pop.colors[grep("^CD68", names(int.43p.pop.colors), invert = T)] <- "grey"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "CEL-seq Source population identities", legend.pos = "right", 
           cols       = int.43p.pop.colors[c(-1,-2,-3)],
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/Fifth round of clustering/Integration original 43p cluster UMAP.pdf", 
           plot.width = 20)

# 10X colored
int.10X.pop.colors <- int.mac.pop.colors
int.10X.pop.colors[grep("^CD68", names(int.10X.pop.colors), invert = F)] <- "grey"

customUMAP(object     = integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "10X Source population identities", legend.pos = "right", 
           cols       = int.10X.pop.colors[c(-1,-2,-3)],
           pt.size    = 5,
           file.name  = "myeloid_43p_celseq_integration/Fifth round of clustering/Integration original 10X cluster UMAP.pdf", 
           plot.width = 20)

p1 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.43p.pop.colors[c(-1,-2,-3)], pt.size = 5) +
  ggtitle("CEL-seq Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.mye.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, cols = int.10X.pop.colors[c(-1,-2,-3)], pt.size = 5) +
  ggtitle("10X Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("myeloid_43p_celseq_integration/Fifth round of clustering/Integration 43p vs. 10X clusters UMAP.pdf", width = 20, height = 10)


#=========================================================================================================================
## Define the clusters: plot markers and gene of interest
##========================================================================================================================
## Scale all genes from RNA assay for diff gene calling and visualisation
# Set default to RNA
DefaultAssay(integrated.mye.seurat) <- "RNA"

# Normalize
integrated.mye.seurat <- NormalizeData(integrated.mye.seurat, normalization.method = "LogNormalize")

# Scale
integrated.mye.seurat <- ScaleData(integrated.mye.seurat, features = row.names(integrated.mye.seurat))

##========================================================================================================================
## Define marker genes per cluster
integrated.mye.seurat.markers <- FindAllMarkers(object = integrated.mye.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save the top9 markers per cluster
sep.markers <- integrated.mye.seurat.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
num.clus <- length(unique(Idents(integrated.mye.seurat)))
for(i in levels(Idents(integrated.mye.seurat))){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = integrated.mye.seurat,
                     name      = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Cluster ",i, " markers", sep = ""),
                     group.by  = "method.int.idents.24.merge",
                     Vln.width = 25, Vln.height = 25)
}

# Save markers per cluster to disk
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  
  # In separate txt files
  write.table(x, file = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
  
  # In one spreadsheet file
  write.xlsx(x, file = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Marker genes.xlsx"), sheetName = i, col.names = T, row.names = F, append = T)
}


##==================================================================================
##==================================================================================
## Define unique and significant markers
uniq.integrated.mye.seurat.markers <- integrated.mye.seurat.markers[!(duplicated(integrated.mye.seurat.markers$gene) | duplicated(integrated.mye.seurat.markers$gene, fromLast = T)),]

# Save unique markers per cluster to disk
for(i in levels(Idents(integrated.mye.seurat))){
  x <- uniq.integrated.mye.seurat.markers[which(uniq.integrated.mye.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Cluster ", i, " unique var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
  
  # In one spreadsheet file
  write.xlsx(x, file = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Unique Marker genes.xlsx"), sheetName = i, col.names = T, row.names = F, append = T)

}

# Save significant unique markers per cluster to disk
for(i in levels(Idents(integrated.mye.seurat))){
  x <- uniq.integrated.mye.seurat.markers[which(uniq.integrated.mye.seurat.markers$cluster == i & uniq.integrated.mye.seurat.markers$p_val_adj < 0.1),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Cluster ", i, " unique and significant var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

# Save significant markers per cluster to disk
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i & integrated.mye.seurat.markers$p_val_adj < 0.1),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Cluster ", i, " significant var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

# Save sig markers per cluster to disk, sorted on log2FC
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i & integrated.mye.seurat.markers$p_val_adj < 0.1),-6]
  x <- x[,c(6,1,2,3,4,5)]
  
  x <- x[order(x$avg_log2FC, decreasing = T),]
  
  # In separate txt files
  write.table(x, file = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Cluster ", i, " sig var genes log2FC sorted.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
  
  # In one spreadsheet file
  write.xlsx(x, file = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Significant Marker genes log2FC sorted.xlsx"), sheetName = i, col.names = T, row.names = F, append = T)

  # Also save the top 10 per cluster only
  write.xlsx(x[1:10,], file = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Significant Marker genes log2FC sorted top 10.xlsx"), sheetName = i, col.names = T, row.names = F, append = T)

}

# Save top 10 sig markers to disk, sorted on log2FC, in one txt file
tmp.markers <- data.frame()
for(i in levels(Idents(integrated.mye.seurat))){
  x <- integrated.mye.seurat.markers[which(integrated.mye.seurat.markers$cluster == i & integrated.mye.seurat.markers$p_val_adj < 0.1),]
  x <- x[order(x$avg_log2FC, decreasing = T),][1:10,]
  x <- x[,c(6,7,1,2,3,4,5)]
  tmp.markers <- rbind(tmp.markers, as.data.frame(x))
}

tmp.markers <- tmp.markers[order(tmp.markers$cluster, tmp.markers$avg_log2FC, decreasing = T),]
write.table(tmp.markers, file = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Marker genes top10 per cluster log2FC sorted.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
rm(tmp.markers)

# plot the top9 unique markers per cluster
for(i in unique(Idents(integrated.mye.seurat))){
  # Make feature, violin, an dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(uniq.integrated.mye.seurat.markers[uniq.integrated.mye.seurat.markers$cluster == i, "gene"]))[1:9], 
                     object    = integrated.mye.seurat,
                     Vln.stack = T,
                     name      = paste("myeloid_43p_celseq_integration/Fifth round of clustering/clusters/Cluster ",i, " unique markers", sep = ""),
                     Vln.color  = c("cyan3","blue4", "blue1"))
}


##========================================================================================================================
# Look at some genes
bunchOfCustomPlots(features   = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                  "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                  "CD79A", "ACTA2", "CD34", "KIT", "OLR1" ), 
                   object     = integrated.mye.seurat,
                   group.by   = "method.int.idents.24.merge",
                   ncol       = 5,
                   name       = "myeloid_43p_celseq_integration/Fifth round of clustering/Cell type defining genes",
                   Vln.width  = 25, 
                   Vln.height = 25
)

# SMC Monocyte genes
bunchOfCustomPlots(features = c("ACTA2", "MYH11", "CD14", "CD68", "KLF4", "CD4", "TAGLN", "ADGRE1", "LGALS3"), 
                   object   = integrated.mye.seurat,
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/SMC mono genes",
                   Vln.width = 25, Vln.height = 15
)

#Extra mono cell markers
bunchOfCustomPlots(features = c("FCGR3A", "CCR2", "CCR5", "SELL", "ITGAX", "CX3CR1"), 
                   object   = integrated.mye.seurat,
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/Extra mono marker genes",
)

# Myeloid Markers
bunchOfCustomPlots(features = c("CD200R1", "TGM2", "CD163"), 
                   object   = integrated.mye.seurat,
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/M2 genes")

bunchOfCustomPlots(features = c("CD2","ITGAM", "PECAM1", "NCAM1", "SELL", "CSF1R"), 
                   object    = integrated.mye.seurat,
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/Monocyte genes",
)

bunchOfCustomPlots(features = c("CD68","CD14","ITGAM","CD1C", "CLEC10A", "FCER1A", "CD2", "CD3E", "CD4"),
                   object   = integrated.mye.seurat,
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/Myeloid genes", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("CASP1","IL1B","TLR4","ABCA1","OLR1", "TREM2","CASP4","KLF4","TNF","IL18","ABCG1", "CD9"), 
                   object   = integrated.mye.seurat,
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/Foam vs. pro inf genes", 
                   ncol     = 3,
                   Vln.width = 30, Vln.height = 15
)

bunchOfCustomPlots(features = c("CD14","CD68","CD86","CD1C", "FCGR2B", "CLEC10A", "FCER1A", "HLA-DQA2", "HLA-DQB2"), 
                   object   = integrated.mye.seurat,
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/DC genes", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("NLRP3", "IL1B", "CASP1", "CASP4", "IL18", "PYCARD"), 
                   object   = integrated.mye.seurat,
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/Inflammasome genes", 
                   ncol     = 3,
)

bunchOfCustomPlots(features = c("CD14", "FCGR3A", "FCGR1A", "CCR2", "PADI4", "SECISBP2L", "CD86", "ITGAX", "LYZ", "CD33", "CD36", "ITGAM"), 
                   object   = integrated.mye.seurat,
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/Mono subtype genes CD14 CD16 CD64 CCR2 PADI4 SLAN CD86 CD11c LYZ CD33 CD36 CD11b", 
                   ncol     = 3,
                   Vln.width = 25, Vln.height = 15
)

bunchOfCustomPlots(features = c("NR4A1"), 
                   object   = integrated.mye.seurat,
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/NR4A1",
                   ncol     = 1,
)

bunchOfCustomPlots(features  = c("CD14","FCGR3A", "CD68", 'CD86', "SELL", "CSF1R", "ITGAM", "ITGAX", "PECAM1"), 
                   object    = integrated.mye.seurat,
                   group.by = "method.int.idents.24.merge",
                   ncol      = 3, Vln.stack = T,
                   name      = "myeloid_43p_celseq_integration/Fifth round of clustering/macro and mono genes",
                   Vln.width = 20, Vln.height = 15
)

bunchOfCustomPlots(features  = c("TREM2", "CD9", "GPNMB", "LYVE1", "MRC1", "FOLR2" ,"OLR1", "ABCA1", "ABCG1", "IL1B", "TNF", "CASP1", "CD14", "FCGR1A", "FCGR3A"), 
                   object    = integrated.mye.seurat,
                   group.by = "method.int.idents.24.merge",
                   ncol      = 3,
                   name      = "myeloid_43p_celseq_integration/Fifth round of clustering/LAM - resident - foam - inflammatory - general markers",
                   Vln.width = 20, Vln.height = 15
)


bunchOfCustomPlots(features = c("TIMD4", "CCR2", "LYVE1", "FOLR2"), 
                   object   = integrated.mye.seurat, assay = "RNA",
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/core_Gene_macs",
                   ncol     = 4, 
                   Vln.width = 25
)

bunchOfCustomPlots(features = c("CX3CR1", "LYVE1", "MRC1", "F13A1", "GAS6", "FOLR2", "TREM2","CD9", "IL1B"), 
                   object   = integrated.mye.seurat, assay = "RNA",
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/resident_macs",
                   ncol     = 3, 
                   Vln.width = 25
)

bunchOfCustomPlots(features = c("CSF1R", "CDK1", "CCNB1", "MKI67", "RAD51", "BCCIP", "BIRC5"), 
                   object   = integrated.mye.seurat, assay = "RNA",
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/proliferation genes",
                   ncol     = 2
)

bunchOfCustomPlots(features = c("CCR2", "ITGAL", "SELL", "ITGA8", "SELPLG", "ITGAV"), 
                   object   = integrated.mye.seurat, assay = "RNA",
                   group.by = "method.int.idents.24.merge",
                   name     = "myeloid_43p_celseq_integration/Fifth round of clustering/motility genes",
                   ncol     = 2
)


customFeature(object = integrated.mye.seurat, features = c("MKI67"), name = "myeloid_43p_celseq_integration/Fifth round of clustering/MKI67.pdf")


##==================================================================================
##==================================================================================
##==================================================================================
## Stratify by expression
dir.create("myeloid_43p_celseq_integration/Stratification", recursive = T, showWarnings = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "IL1B", file.name = "myeloid_43p_celseq_integration/Stratification/IL1B stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "OLR1", file.name = "myeloid_43p_celseq_integration/Stratification/OLR1 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "TREM2", file.name = "myeloid_43p_celseq_integration/Stratification/TREM2 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "FOLR2", file.name = "myeloid_43p_celseq_integration/Stratification/FOLR2 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "FCGR3A", file.name = "myeloid_43p_celseq_integration/Stratification/FCGR3A stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "CD9", file.name = "myeloid_43p_celseq_integration/Stratification/CD9 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "ABCG1", file.name = "myeloid_43p_celseq_integration/Stratification/ABCG1 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "MRC1", file.name = "myeloid_43p_celseq_integration/Stratification/MRC1 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "LYVE1", file.name = "myeloid_43p_celseq_integration/Stratification/LYVE1 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "ABCA1", file.name = "myeloid_43p_celseq_integration/Stratification/ABCA1 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "TNF", file.name = "myeloid_43p_celseq_integration/Stratification/TNF stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "NLRP3", file.name = "myeloid_43p_celseq_integration/Stratification/NLRP3 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "STAT1", file.name = "myeloid_43p_celseq_integration/Stratification/STAT1 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "CASP1", file.name = "myeloid_43p_celseq_integration/Stratification/CASP1 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "FCGR1A", file.name = "myeloid_43p_celseq_integration/Stratification/CASP1 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "STAT3", file.name = "myeloid_43p_celseq_integration/Stratification/STAT3 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "MX1", file.name = "myeloid_43p_celseq_integration/Stratification/MX1 stratified", return.object = F)

stratifyByExpression(object   = integrated.mye.seurat,
                     strat.by = "IRF9", file.name = "myeloid_43p_celseq_integration/Stratification/IRF9 stratified", return.object = F)



### Plot for fig panel
# Change ident order for plotting
tmp.idents <- Idents(integrated.mye.seurat)
levels(tmp.idents)
tmp.idents <- factor(tmp.idents, levels = levels(tmp.idents)[c(7,5,4,1,6,3,2,10,9,8)])
Idents(integrated.mye.seurat) <- tmp.idents

# And plot
bunchOfCustomPlots(features   = c("TREM2", "LYVE1", "OLR1",  "IL1B",  "CD14", 
                                  "CD9",   "MRC1",  "ABCA1", "TNF",   "FCGR1A",
                                  "GPNMB", "FOLR2", "ABCG1", "CASP1", "FCGR3A"), 
                   object     = integrated.mye.seurat, 
                   Vln.color  = M.int_refined.pop.colors,
                   ncol       = 5,
                   name       = "myeloid_43p_celseq_integration/LAM - resident - foam - inflammatory - general markers pop colors",
                   Vln.draw.names = F,
                   Vln.width  = 30, Vln.height = 15
)


bunchOfCustomPlots(features = c("CSF1R", "CDK1", "CCNB1", "MKI67", "RAD51", "BCCIP", "BIRC5"), 
                   object   = integrated.mye.seurat,
                   Vln.color  = M.int_refined.pop.colors,
                   name     = "myeloid_43p_celseq_integration/proliferation genes",
                   Vln.draw.names = F,
                   ncol     = 3, Vln.width = 15
)

bunchOfCustomPlots(features = c("ITGAM", "ITGAL", "CD99", "ITGB1", "ITGB2", "ITGA4", "RAC1", "RAC2", "CDC42"), 
                   object   = integrated.mye.seurat,
                   Vln.color  = M.int_refined.pop.colors,
                   name     = "myeloid_43p_celseq_integration/motility genes",
                   Vln.draw.names = F, 
                   Vln.width = 20,
                   ncol     = 3
)

customDot(features = c("ITGAM", "ITGAL", "CD99", "ITGB1", "ITGB2", "ITGA4", "RAC1", "RAC2", "CDC42"), 
          object   = integrated.mye.seurat,
          cluster.idents = F,
          name     = "myeloid_43p_celseq_integration/motility genes - dotplot.pdf")

bunchOfCustomPlots(features = c("TMEM86A"), 
                   object   = integrated.mye.seurat,
                   Vln.color  = M.int_refined.pop.colors,
                   name     = "myeloid_43p_celseq_integration/TMEM86A",
                   Vln.draw.names = F)

bunchOfCustomPlots(features = c("TMEM86A"), 
                   object   = mye.43p.seurat,
                   Vln.color  = M.int_refined.pop.colors,
                   name     = "myeloid_43p_celseq_integration/TMEM86A 43p",
                   Vln.draw.names = F)

stratifyByExpression(object = integrated.mye.seurat, strat.by = "TMEM86A", file.name = "myeloid_43p_celseq_integration/TMEM86A quartiles.pdf", return.object = F)
stratifyByExpression(object = mye.43p.seurat, strat.by = "TMEM86A", file.name = "myeloid_43p_celseq_integration/TMEM86A 43p quartiles.pdf", return.object = F)
tmp.TMEM.strat <- stratifyByExpression(object = mye.43p.seurat, strat.by = "TMEM86A", file.name = "myeloid_43p_celseq_integration/TMEM86A 43p quartiles.pdf", return.object = T, do.plot = F)

stratPlots(object   = tmp.TMEM.strat, 
           group.by = "TMEM86A_expr", 
           features = c("ABCA1", "ABCG1", "OLR1","FOLR2", "MYLIP", "MARCHF6"), 
           cols     = c("grey", "bisque", "coral", "firebrick"),
           x.label  = paste("Stratified by", "TMEM86A" ,"expression", sep = " "), 
           name     = paste("myeloid_43p_celseq_integration/TMEM86A 43p", " - Foam and Resident genes 2", sep = ""),
           ncol     = 1)


customUMAP(object = mye.43p.seurat, file.name =  "myeloid_43p_celseq_integration/43p UMAP.pdf", legend.pos = "right", plot.width = 15, cols = int.43p.pop.colors[13:17])

customFeature(object = mye.43p.seurat, features = "TMEM86A", pt.size = 5, order = T, name = "myeloid_43p_celseq_integration/TMEM86A - feature plot.pdf")
customFeature(object = mye.43p.seurat, features = "MARCHF6", pt.size = 4, order = T, name = "myeloid_43p_celseq_integration/MARCHF6 - feature plot.pdf")
customFeature(object = mye.43p.seurat, features = "EEPD1", pt.size = 4, order = T, name = "myeloid_43p_celseq_integration/EEPD1 - feature plot.pdf")

## Check cell cycle
customUMAP(object = integrated.mye.seurat, group.by = "Phase", pt.size = 4, title = "Cell cycle", legend.pos = "right", file.name = "myeloid_43p_celseq_integration/cell cycle UMAP.pdf")

