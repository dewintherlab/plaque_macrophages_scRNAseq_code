#=========================================================================================================================
## Integrate  all cells with the 43 patient CEL-seq myeloid cells
##========================================================================================================================
dir.create("full_43p_celseq_integration/First round of clustering/clusters", showWarnings = F, recursive = T)

## Load the data
full.43p.seurat <- readRDS(file = "../SORTseq_Patients/Seurat Objects/main.seurat.filtered.clustered.RDS")
full.43p.seurat
unique(Idents(full.43p.seurat))


## Prep for integration
# Select the corect assay
DefaultAssay(full.43p.seurat) <- "RNA"
DefaultAssay(patients.seurat) <- "RNA"

# Make a list of our objects to integrate
integrated.full.seurat <- list(patients.seurat, full.43p.seurat)

# Define var features
integrated.full.seurat <- lapply(X = integrated.full.seurat, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = integrated.full.seurat)

# Scale data and run PCA using the selected features
integrated.full.seurat <- lapply(X = integrated.full.seurat, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

## Integrate the data
## Use RPCA  to speed things up.
# First select anchors for integration (by RPCA reduction)
anchors <- FindIntegrationAnchors(object.list = integrated.full.seurat, anchor.features = features, reduction = "rpca")

# And do the actual integrating
integrated.full.seurat <- IntegrateData(anchorset = anchors)

# Show stats
integrated.full.seurat

#=========================================================================================================================
## Update metadata
head(integrated.full.seurat@meta.data)
tail(integrated.full.seurat@meta.data)

# Fix tissue column from 43p to 'plaque'
integrated.full.seurat@meta.data[is.na(integrated.full.seurat$Tissue),"Tissue"] <- "plaque"

# Fix sex column from 10X to 'male'
integrated.full.seurat@meta.data[is.na(integrated.full.seurat$Sex),"Sex"] <- "male"

# Add a new column 'Method' to distinguish the two datasets
method <- rep("CEL-seq", dim(integrated.full.seurat)[2])
method[grep("^P", integrated.full.seurat$Patient)] <- "10X"
integrated.full.seurat <- AddMetaData(integrated.full.seurat, method, col.name = "Method")

## Add a new column 'original.refined.ident' to capture the original idents of the unintegrated objects.
# Set up destination dataframe
ref.ident  <- data.frame(row.names = row.names(integrated.full.seurat@meta.data))

# Fetch idents and stick 'em in a dataframe for merging
idents.43p <- data.frame(refined.ident = Idents(full.43p.seurat))
idents.10X <- data.frame(refined.ident = Idents(patients.seurat))
idents.com <- rbind(idents.43p, idents.10X)

# Merge the frames together
ref.ident <- merge(ref.ident, idents.com, by = 0)
row.names(ref.ident) <- ref.ident$Row.names
ref.ident$Row.names <- NULL
tail(ref.ident)

# And add the metadata
integrated.full.seurat <- AddMetaData(integrated.full.seurat, metadata = ref.ident, col.name = "original.refined.ident")
                                      

#=========================================================================================================================
## Run the embeddings, cluster, and have a look
##========================================================================================================================
# Run dimensionality reductions
integrated.full.seurat <- ScaleData(object = integrated.full.seurat, verbose = F)
integrated.full.seurat <- RunPCA(   object = integrated.full.seurat, verbose = F)
integrated.full.seurat <- RunUMAP(  object = integrated.full.seurat, reduction = "pca", dims = 1:30)

# Cluster
set.seed(1)
integrated.full.seurat <- FindNeighbors(object = integrated.full.seurat, reduction = "pca", dims = 1:30, force.recalc = T, verbose = T)
integrated.full.seurat <- FindClusters( object = integrated.full.seurat, resolution = 1.25, verbose = T, random.seed = 666)

## Visualize
# Plot Clusters
customUMAP(object     = integrated.full.seurat, shuffle = T,
           title      = "Initial clustering", 
           legend.pos = "right", 
           pt.size    = 2,
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration intitial clustering UMAP.pdf", 
           plot.width = 15)

# Plot Method
customUMAP(object     = integrated.full.seurat, shuffle = T,
           title      = "Initial clustering", 
           legend.pos = "right", group.by = "Method",
           pt.size    = 1,
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration Method UMAP.pdf", 
           plot.width = 15)

# Plot Patients
customUMAP(object     = integrated.full.seurat, shuffle = T,
           title      = "Initial clustering", 
           legend.pos = "right", group.by = "Patient",
           pt.size    = 1,
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration Patient UMAP.pdf", 
           plot.width = 15)

# Plot Tissue
customUMAP(object     = integrated.full.seurat, shuffle = T,
           title      = "Initial clustering", 
           legend.pos = "right", group.by = "Tissue",
           pt.size    = 1,
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration Tissue UMAP.pdf", 
           plot.width = 15)



## Plot the original 10X and 43p idents
customUMAP(object     = integrated.full.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "Source clusters", legend.pos = "none", 
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration original cluster UMAP.pdf", 
           plot.width = 15)

p1 <- DimPlot(integrated.full.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, pt.size = 3) +
  ggtitle("Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.full.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T, pt.size = 3) +
  ggtitle("Integrated clusters") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("full_43p_celseq_integration/First round of clustering/Integration original vs. integrated cluster UMAP.pdf", width = 20, height = 10)

# Highlight each source cluster
dir.create("full_43p_celseq_integration/First round of clustering/Original source cluster highlights", showWarnings = F, recursive = T)
for(theIdent in unique(integrated.full.seurat$original.refined.ident)){
  customUMAP(object          = integrated.full.seurat, reduction = "umap", label = F, shuffle = T, title = theIdent, legend.pos = "none", 
             file.name       = paste("full_43p_celseq_integration/First round of clustering/Original source cluster highlights/", theIdent, " UMAP.pdf", sep = ""),
             sizes.highlight = 4,
             cells.highlight = row.names(ref.ident)[which(ref.ident$refined.ident == theIdent)])
}

## Color only one source at a time
# 43p colored
cols.43p        <- colorSpacer(startcolor = "red", endcolor = "blue", steps = length(unique(Idents(subset(integrated.full.seurat, Method != "10X")))))
names(cols.43p) <- unique(Idents(subset(integrated.full.seurat, Method != "10X")))
customUMAP(object     = integrated.full.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "CEL-seq Source population identities", legend.pos = "none", 
           pt.size    = 5,
           cols       = cols.43p,
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration original 43p cluster UMAP.pdf", 
           plot.width = 15)

# 10X colored
cols.10X        <- colorSpacer(startcolor = "yellow", endcolor = "green", steps = length(unique(Idents(subset(integrated.full.seurat, Method == "10X")))))
names(cols.10X) <- unique(Idents(subset(integrated.full.seurat, Method == "10X")))
customUMAP(object     = integrated.full.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "10X Source population identities", legend.pos = "none", 
           pt.size    = 5,
           cols       = cols.10X,
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration original 10X cluster UMAP.pdf", 
           plot.width = 15)

# Double plot
p1 <- DimPlot(integrated.full.seurat, reduction = "umap", cols = cols.43p, group.by = "original.refined.ident", shuffle = T, pt.size = 5) +
  ggtitle("CEL-seq Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.full.seurat, reduction = "umap", cols = cols.10X, group.by = "original.refined.ident", shuffle = T, pt.size = 5) +
  ggtitle("10X Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("full_43p_celseq_integration/First round of clustering/Integration 43p vs. 10X clusters UMAP.pdf", width = 20, height = 10)

## Save the integrated seurat object
saveRDS(integrated.full.seurat, "Seurat_Objects/full.43p_10X.integrated.seurat.backup.RDS")
#integrated.full.seurat <- readRDS(file = "Seurat_Objects/full.43p_10X.integrated.seurat.backup.RDS")


#==================================================================================================================================
## Name all clusters so we can do a proper LRI analysis
# First let's make a bit less clusters (especially T cell cloud is too divided)
dir.create("full_43p_celseq_integration/Second round of clustering/", showWarnings = F, recursive = T)
set.seed(1)
integrated.full.seurat <- FindClusters( object = integrated.full.seurat, resolution = 1, verbose = T, random.seed = 666)

## Visualize
# Plot Clusters
customUMAP(object     = integrated.full.seurat, shuffle = T,
           title      = "Initial clustering", label = T, label.size = 4,
           legend.pos = "right", 
           pt.size    = 2,
           file.name  = "full_43p_celseq_integration/Second round of clustering/Integration intitial clustering UMAP.pdf", 
           plot.width = 15)

# Plot Method
customUMAP(object     = integrated.full.seurat, shuffle = T,
           title      = "Initial clustering", 
           legend.pos = "right", group.by = "Method",
           pt.size    = 1,
           file.name  = "full_43p_celseq_integration/Second round of clustering/Integration Method UMAP.pdf", 
           plot.width = 15)

# Plot Patients
customUMAP(object     = integrated.full.seurat, shuffle = T,
           title      = "Initial clustering", 
           legend.pos = "right", group.by = "Patient",
           pt.size    = 1,
           file.name  = "full_43p_celseq_integration/Second round of clustering/Integration Patient UMAP.pdf", 
           plot.width = 15)

# Plot Tissue
customUMAP(object     = integrated.full.seurat, shuffle = T,
           title      = "Initial clustering", 
           legend.pos = "right", group.by = "Tissue",
           pt.size    = 1,
           file.name  = "full_43p_celseq_integration/Second round of clustering/Integration Tissue UMAP.pdf", 
           plot.width = 15)


## Plot some marker genes
DefaultAssay(integrated.full.seurat) <- "RNA"
bunchOfCustomPlots(features   = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                  "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                  "CD79A", "ACTA2", "CD34", "KIT", "OLR1" ), 
                   object     = integrated.full.seurat,
                   ncol       = 3,
                   name       = "full_43p_celseq_integration/Second round of clustering/Cell type defining genes",
                   Vln.width  = 20, 
                   Vln.height = 25
)

## Plot some t cell and NK genes
bunchOfCustomPlots(features   = c("CD3E", "CD4", "CD8A", "FOXP3", "PTPRC", "CD27", "TRGC1", "TRGC2", "TRDC", "NCAM1", "KLRC1", "FCGR3A", "KLRD1", "CD79A", "CD14", "GNLY", "GZMB", "GZMA"), 
                   object     = integrated.full.seurat,
                   ncol       = 3,
                   name       = "full_43p_celseq_integration/Second round of clustering/T Cell type defining genes",
                   Vln.width  = 20, 
                   Vln.height = 25
)

## Plot some t cell and NK genes
bunchOfCustomPlots(features   = c("GYPA", "TFRC", "HBB"), 
                   object     = integrated.full.seurat,
                   ncol       = 3,
                   name       = "full_43p_celseq_integration/Second round of clustering/Ery type defining genes",
                   Vln.width  = 20, 
                   Vln.height = 25
)

## Plot some t cell and NK genes
bunchOfCustomPlots(features   = c("CD1C", "CLEC4A", "CLEC9A", "CLEC10A", "CD14","CD68","CD86", "FCGR2B", "FCER1A", "HLA-DQA2", "HLA-DQB2", "ITGAX"), 
                   object     = integrated.full.seurat,
                   ncol       = 3,
                   name       = "full_43p_celseq_integration/Second round of clustering/DC type defining genes",
                   Vln.width  = 20, 
                   Vln.height = 25
)


## Get markers per cluster
integrated.full.seurat.markers <- FindAllMarkers(integrated.full.seurat, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save the top9 markers per cluster
sep.markers <- integrated.full.seurat.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
dir.create("full_43p_celseq_integration/Second round of clustering/clusters/", showWarnings = F, recursive = T)
num.clus <- length(unique(Idents(integrated.full.seurat)))
for(i in levels(Idents(integrated.full.seurat))){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = integrated.full.seurat,
                     name      = paste("full_43p_celseq_integration/Second round of clustering/clusters/Cluster ",i, " markers", sep = ""),
                     
                     Vln.width = 25, Vln.height = 25)
}


## Use this info to rename all but the myeloid clusters
integrated.full.seurat <- RenameIdents(integrated.full.seurat, "0"  = "CD4+ T Cells",
                                                               "1"  = "CD4+ T Cells",
                                                               "2"  = "CD4+ T Cells",
                                                               "3"  = "CD8+ T Cells",
                                                               "4"  = "CD4+ T Cells",
                                                               "5"  = "CD8+ T Cells",
                                                               "7"  = "CD56-CD16+ NK Cells",
                                                               "9"  = "CD79+ Class-switched Memory B Cells",
                                                               "10" = "CD8+ T Cells",
                                                               "11" = "CD8+ T Cells",
                                                               "12" = "FOXP3+ T Cells",
                                                               "13" = "CD56+KLRC1+ NK Cells",
                                                               "15" = "CD1C+ cDC1",
                                                               "17" = "CD4+ T Cells",
                                                               "21" = "CD4+ T Cells",
                                                               "16" = "CD34+ Endothelial Cells",
                                                               "18" = "ACTA2+ Smooth Muscle Cells",
                                                               "19" = "CD3+MKI67+ Proliferating T Cells",
                                                               "20" = "GYPA+ Erythroid Cells",
                                                               "23" = "KIT+ Mast Cells",
                                                               "24" = "CD79+ Plasma B Cells",
                                                               "26" = "CD79+ B Cells",
                                                               "27" = "CLEC9A+ cDC2"
                                                               )

# Plot Clusters
customUMAP(object     = integrated.full.seurat, shuffle = T,
           title      = "Initial clustering", label = T, label.size = 4,
           legend.pos = "right", 
           pt.size    = 2,
           file.name  = "full_43p_celseq_integration/Second round of clustering/Integration named clustering no mye UMAP.pdf", 
           plot.width = 15)



## Incorporate our refined monocyte and macrophage identities
# Get myeloid idents
mye.idents        <- c(as.vector(Idents(integrated.mye.seurat)), 
                       as.vector(Idents(pbmc.monocytes.seurat)))
names(mye.idents) <- c(names(Idents(integrated.mye.seurat)),
                       names(Idents(pbmc.monocytes.seurat)))
mye.cells         <- names(mye.idents)


# Check where the mye idents fall on the full plot
p1 <- DimPlot(integrated.full.seurat, cells = mye.cells ) + theme(legend.position = "bottom")
p2 <- DimPlot(integrated.full.seurat) + theme(legend.position = "bottom")
p1 + p2
ggsave("full_43p_celseq_integration/Second round of clustering/Mye cell identities UMAP.pdf", width = 20, height = 10)

tmp.seurat <- AddMetaData(integrated.full.seurat, metadata = mye.idents, col.name = "mye.ident")
DimPlot(tmp.seurat, cells = mye.cells, group.by = "mye.ident", cols = M.int_refined.pop.colors) + theme(legend.position = "right")
ggsave("full_43p_celseq_integration/Second round of clustering/Original refined Mye cell identities UMAP.pdf", width = 20, height = 10)

# Set myeloid idents
Idents(integrated.full.seurat, cells = mye.cells) <- mye.idents

# Visualize
customUMAP(object = integrated.full.seurat, legend.pos = "right", plot.width = 20, pt.size = 1,
           file.name = "full_43p_celseq_integration/Second round of clustering/Refined mye pops UMAP.pdf")

## Refine populations by deleting cells left out of the 'pure' myeloid analyses
still.numbered        <- as.vector(Idents(integrated.full.seurat))
names(still.numbered) <- names(Idents(integrated.full.seurat))
still.numbered        <- still.numbered[still.numbered %in% c(6, 8, 14, 15, 22, 25, 27)]
customUMAP(object = integrated.full.seurat, legend.pos = "right", plot.width = 20, pt.size = 1, cells = names(still.numbered), file.name =  "full_43p_celseq_integration/Second round of clustering/Refined mye pops remaining numbers UMAP.pdf")

# Set idents
integrated.full.seurat <- RenameIdents(integrated.full.seurat, "6"  = "CD14+CD16-CD64+SLAN- Classical Monocytes")
customUMAP(object = integrated.full.seurat, legend.pos = "right", plot.width = 20, pt.size = 1, file.name =  "full_43p_celseq_integration/Second round of clustering/Refined mye pops added clas mons remaining numbers UMAP.pdf")

# Remove left over cells
integrated.full.seurat <- subset(integrated.full.seurat, cells = names(still.numbered), invert = T)

## Set up a new all inclusive color table
levels(Idents(integrated.full.seurat))[14:length(levels(Idents(integrated.full.seurat)))]
other.colors        <- c("CD4+ T Cells" = "cornsilk2", 
                         "CD8+ T Cells" = "burlywood3", 
                         "FOXP3+ T Cells" = "cornsilk4", 
                         "CD3+MKI67+ Proliferating T Cells" = "antiquewhite1", 
                         "CD79+ B Cells" = "deeppink3", 
                         "CD79+ Class-switched Memory B Cells" = "deeppink", 
                         "CD79+ Plasma B Cells" = "deeppink4", 
                         "CD56-CD16+ NK Cells" = "aquamarine3", 
                         "CD56+KLRC1+ NK Cells" = "aquamarine4", 
                         "CD1C+ cDC1" = "darkgoldenrod1", 
                         "CLEC9A+ cDC2" = "gold1", 
                         "GYPA+ Erythroid Cells" = "coral",
                         "KIT+ Mast Cells" = "coral3",
                         "CD34+ Endothelial Cells" = "chocolate2",
                         "ACTA2+ Smooth Muscle Cells" = "brown3")
full_set.colors     <- c(M.int_refined.pop.colors, other.colors)

## Visualize
customUMAP(object = integrated.full.seurat, legend.pos = "right", cols = full_set.colors, plot.width = 20, pt.size = 1, file.name =  "full_43p_celseq_integration/Second round of clustering/Refined mye pops cleaned UMAP.pdf")

## Draw a more 'top level' view
toplevelidents <- as.vector(Idents(integrated.full.seurat))
unique(toplevelidents)
toplevelidents[grep("B Cell",      toplevelidents)] <- "B Cells"
toplevelidents[grep("T Cells",     toplevelidents)] <- "T Cells"
toplevelidents[grep("NK Cells",    toplevelidents)] <- "NK Cells"
toplevelidents[grep("cDC",         toplevelidents)] <- "Dendritic Cells"
toplevelidents[grep("Macrophage",  toplevelidents)] <- "Macrophages"
toplevelidents[grep("Monocytes",   toplevelidents)] <- "Monocytes"
toplevelidents[grep("Endothelial", toplevelidents)] <- "Endothelial Cells"
toplevelidents[grep("Smooth",      toplevelidents)] <- "Smooth Muscle Cells"
toplevelidents[grep("Mast",        toplevelidents)] <- "Mast Cells"
toplevelidents[grep("Erythroid",   toplevelidents)] <- "Erythroid Cells"
unique(toplevelidents)

# Add to the seurat object
integrated.full.seurat <- AddMetaData(integrated.full.seurat, metadata = toplevelidents, col.name = "top.level.idents")

# Set colors
full_set.top.colors <- c("B Cells" = "yellow3", "Dendritic Cells" = "turquoise4", "Endothelial Cells" = "orange3", "Erythroid Cells" = "red2", "Macrophages" = "skyblue", "Mast Cells" = "purple4", "NK Cells" = "chartreuse4", "Smooth Muscle Cells" = "brown4", "T Cells" = "pink3", "Monocytes" = "blue")

# Plot total UMAP
customUMAP(object = integrated.full.seurat, legend.pos = "right", group.by = "top.level.idents", cols = full_set.top.colors, plot.width = 10, pt.size = 1, file.name =  "full_43p_celseq_integration/Second round of clustering/Top level Idents full UMAP.pdf")

# Plot plaque cell UMAP


## Save the integrated seurat object
saveRDS(integrated.full.seurat, "Seurat_Objects/full.43p_10X.integrated.cleaned.seurat.RDS")
#integrated.full.seurat <- readRDS(file = "Seurat_Objects/full.43p_10X.integrated.cleaned.seurat.RDS")
