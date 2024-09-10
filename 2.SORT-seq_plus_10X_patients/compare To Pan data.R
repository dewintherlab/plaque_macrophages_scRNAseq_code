#=================================================================================================================================
## Comparison of our data with Pan et al.
## Specifically, can we map our populations and do we see the same trajectories.
##================================================================================================================================
dir.create("Pan_validation_results", showWarnings = F)

## Load the pan data
# Load the object into the env
load("raw_data/Pan.seurat.RData", verbose = T)

# Clean up
pan.seurat <- seuset
rm(project.name)
rm(seuset)
pan.seurat <- UpdateSeuratObject(pan.seurat)
pan.seurat <- RunPCA(pan.seurat)
pan.seurat <- RunUMAP(pan.seurat, dims = 1:30, return.model = T)

# Check it out
pan.seurat
unique(Idents(pan.seurat))
row.names(GetAssayData(pan.seurat))

# Plot the given populations
DimPlot(pan.seurat)
customUMAP(object = pan.seurat, pt.size = 1.5, label = F, title = "Pan et al. data", shuffle = T, seed = 666, file.name = "Pan_validation_results/Pan.UMAP.pdf")


## Transfer our labels by reference mapping
# Find anchors
pan.anchors <- FindTransferAnchors(reference           = integrated.full.seurat,
                                   query               = pan.seurat,
                                   dims                = 1:30, 
                                   reference.reduction = "pca")

# Predict labels
predictions <- TransferData(anchorset = pan.anchors, 
                            refdata   = Idents(integrated.full.seurat),
                            dims      = 1:30)

# Add to the pan seurat object
pan.seurat <- AddMetaData(pan.seurat, metadata = predictions)

# Map to ref UMAP
integrated.full.seurat <- RunUMAP(integrated.full.seurat, dims = 1:50, return.model = T)
pan.seurat <- MapQuery(anchorset           = pan.anchors, 
                       reference           = integrated.full.seurat, 
                       query               = pan.seurat,
                       refdata             = Idents(integrated.full.seurat), 
                       reference.reduction = "pca", 
                       reduction.model     = "umap")

# Visualize
customUMAP(object = pan.seurat, group.by = "predicted.id", pt.size = 1.5, label = F, title = "Pan et al. ref mapped", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "Pan_validation_results/Pan.full_ref_map.UMAP.pdf", cols = full_set.colors, plot.width = 15)
customUMAP(object = pan.seurat, group.by = "predicted.id", reduction = "ref.umap", pt.size = 1.5, label = F, title = "Pan et al. ref UMAP mapped",
           shuffle = T, seed = 666, file.name = "Pan_validation_results/Pan.full_ref_umap.UMAP.pdf", legend.pos = "right", cols = full_set.colors, plot.width = 15)

## Now subset the macropahges and try again
# Subset the pan object
mye.pan.seurat <- subset(pan.seurat, idents = c("Macrophage 1", "Macrophage 2", "Macrophage 3"))

# Re-center
mye.pan.seurat <- RunPCA(mye.pan.seurat)
mye.pan.seurat <- RunUMAP(mye.pan.seurat, dims = 1:30, return.model = T)

# Check it out
mye.pan.seurat
unique(Idents(mye.pan.seurat))

# Plot the given populations
DimPlot(mye.pan.seurat)
customUMAP(object = mye.pan.seurat, pt.size = 1.5, label = F, title = "Pan et al. data", shuffle = T, seed = 666, file.name = "Pan_validation_results/Pan.mye.UMAP.pdf")


## Transfer our labels by reference mapping
# Find anchors
DefaultAssay(integrated.mye.seurat) <- "integrated"
mye.pan.anchors <- FindTransferAnchors(reference           = integrated.mye.seurat,
                                       query               = mye.pan.seurat,
                                       dims                = 1:50, 
                                       reference.reduction = "pca")

# Predict labels
predictions <- TransferData(anchorset = mye.pan.anchors, 
                            refdata   = Idents(integrated.mye.seurat),
                            dims      = 1:50)

# Add to the pan seurat object
mye.pan.seurat <- AddMetaData(mye.pan.seurat, metadata = predictions)

# Map to ref UMAP
integrated.mye.seurat <- RunUMAP(integrated.mye.seurat, dims = 1:50, return.model = T)
mye.pan.seurat <- MapQuery(anchorset       = mye.pan.anchors, 
                       reference           = integrated.mye.seurat, 
                       query               = mye.pan.seurat,
                       refdata             = Idents(integrated.mye.seurat), 
                       reference.reduction = "pca", 
                       reduction.model     = "umap")

# Visualize
customUMAP(object = mye.pan.seurat, group.by = "predicted.id", pt.size = 1.5, label = F, title = "Pan et al. ref mapped", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "Pan_validation_results/Pan.mye_ref_map.UMAP.pdf", cols = M.int_refined.pop.colors, plot.width = 15)
customUMAP(object = mye.pan.seurat, group.by = "predicted.id", reduction = "ref.umap", pt.size = 1.5, label = F, title = "Pan et al. ref UMAP mapped",
           shuffle = T, seed = 666, file.name = "Pan_validation_results/Pan.mye_ref_umap.UMAP.pdf", legend.pos = "right", cols = M.int_refined.pop.colors, plot.width = 15)
customUMAP(object = integrated.mye.seurat,  pt.size = 1.5, label = F, legend.pos = "right", 
           shuffle = T, seed = 666, file.name = "Pan_validation_results/Mye.UMAP.pdf", cols = M.int_refined.pop.colors, plot.width = 15)

bunchOfCustomPlots(features   = c("TREM2", "LYVE1", "OLR1",  "IL1B",  "CD14", 
                                  "CD9",   "MRC1",  "ABCA1", "TNF",   "FCGR1A",
                                  "GPNMB", "FOLR2", "ABCG1", "CASP1", "FCGR3A"), 
                   object     = mye.pan.seurat, 
                   group.by   = "predicted.id",
                   Vln.color  = M.int_refined.pop.colors,
                   ncol       = 5,
                   name       = "Pan_validation_results/LAM - resident - foam - inflammatory - general markers pop colors",
                   Vln.draw.names = F,
                   Vln.width  = 30, Vln.height = 15
)

stratifyByExpression(object = mye.pan.seurat, strat.by = "OLR1",  file.name = "Pan_validation_results/Mye.UMAP.OLR1.pdf",  onlyUMAP = T)
stratifyByExpression(object = mye.pan.seurat, strat.by = "ABCA1", file.name = "Pan_validation_results/Mye.UMAP.ABCA1.pdf", onlyUMAP = T)
stratifyByExpression(object = mye.pan.seurat, strat.by = "ABCG1", file.name = "Pan_validation_results/Mye.UMAP.ABCG1.pdf", onlyUMAP = T)
FeaturePlot(object = mye.pan.seurat, features = c("OLR1", "ABCA1", "ABCG1"), pt.size = 1, order = T)
ggsave("Pan_validation_results/Mye.UMAP.foamy.pdf")
