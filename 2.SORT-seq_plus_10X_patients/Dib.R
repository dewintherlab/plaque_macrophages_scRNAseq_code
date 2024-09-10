#======================================================================
## Compare with Dib et al. myeloid clusters
#======================================================================
## Create dir
dir.create("Dib_integration", showWarnings = F)

## Load the objects
# Seurat
workdir <- "raw_data/Dib_from_myel_clean_cluster.0.8.dir/"
dib.seurat <- readRDS(paste0(workdir,"begin.rds"))
dib.seurat

# UMAP coords
umap <- read.delim(paste0(workdir,"umap.0.5.tsv.gz"))
head(umap)

# Idents
cl <- read.table(paste0(workdir,"cluster_assignments.tsv.gz"), sep="\t", header=TRUE)
head(cl)
table(cl$cluster_id)

# UMAP coords embedding test
meta <- dib.seurat@meta.data
meta$barcode <- row.names(meta)
meta <- merge(meta, umap, by = 'barcode')
meta <- merge(meta, cl, by = 'barcode')
meta <- meta %>% mutate(cluster_id = factor(cluster_id))
head(meta)

ggplot(meta, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(aes(color=cluster_id), size=1.1) +
  theme_classic() +
  labs(x="UMAP 1", y="UMAP 2") +
  labs(color = "Cluster id")


# Add idents to Seurat object
cl.vec             <- cl$cluster_id
names(cl.vec)      <- cl$barcode
Idents(dib.seurat) <- cl.vec
head(Idents(dib.seurat))

## Prep for integration
dib.seurat <- UpdateSeuratObject(dib.seurat)
dib.seurat <- SCTransform(dib.seurat)
dib.seurat <- RunPCA(dib.seurat)
umap.mat <- as.matrix(meta[,c("UMAP_1", "UMAP_2")])
row.names(umap.mat) <- meta$barcode
dib.seurat[["umap"]] <- CreateDimReducObject(embeddings = umap.mat, assay = "RNA", key = "UMAP_")
 
# Plot the given populations
customUMAP(object = dib.seurat, pt.size = 1.2, label = F, title = "Dib et al. data", shuffle = T, seed = 666, plot.width = 10, file.name = "Dib_integration/dib.UMAP.pdf")

## Mac & Mono
## Transfer our labels by reference mapping
# Find anchors
DefaultAssay(dib.seurat) <- "RNA"
dib.anchors <- FindTransferAnchors(reference           = final.pop.call.from_full.integrated.mye.seurat,
                                   query               = dib.seurat,
                                   dims                = 1:30, 
                                   reference.reduction = "pca")

# Predict labels
predictions <- TransferData(anchorset = dib.anchors, 
                            refdata   = Idents(final.pop.call.from_full.integrated.mye.seurat),
                            dims      = 1:30)

# Add to the Dib seurat object
dib.seurat <- AddMetaData(dib.seurat, metadata = predictions)

# Map to ref UMAP
#final.pop.call.from_full.integrated.mye.seurat <- RunUMAP(final.pop.call.from_full.integrated.mye.seurat, dims = 1:30, reduction = "pca", return.model = TRUE)
dib.seurat <- MapQuery(anchorset           = dib.anchors, 
                       reference           = final.pop.call.from_full.integrated.mye.seurat, 
                       query               = dib.seurat,
                       refdata             = Idents(final.pop.call.from_full.integrated.mye.seurat), 
                       reference.reduction = "pca", 
                       reduction.model     = "umap")

head(Idents(dib.seurat))

# Visualize
customUMAP(object = final.pop.call.from_full.integrated.mye.seurat, pt.size = 1.2, label = F, title = "Mac & Mono UMAP", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "Dib_integration/full_ref.UMAP.pdf", cols = full_set.colors, plot.width = 15)
customUMAP(object = dib.seurat, group.by = "predicted.id", pt.size = 1.2, label = F, title = "Dib et al. ref mapped", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "Dib_integration/dib.full_ref_map.UMAP.pdf", cols = full_set.colors, plot.width = 25)
customUMAP(object = dib.seurat, pt.size = 1.2, label = F, title = "Dib et al. ref mapped", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "Dib_integration/dib.full_ref_map.dib_idents.UMAP.pdf", plot.width = 15)
customUMAP(object = dib.seurat, group.by = "predicted.id", reduction = "ref.umap", pt.size = 1.2, label = F, title = "Dib et al. ref UMAP mapped",
           shuffle = T, seed = 666, file.name = "Dib_integration/dib.full_ref_umap.UMAP.pdf", legend.pos = "right", cols = full_set.colors, plot.width = 15)
customUMAP(object = dib.seurat, reduction = "ref.umap", pt.size = 1.2, label = F, title = "Dib et al. ref UMAP mapped",
           shuffle = T, seed = 666, file.name = "Dib_integration/dib.full_ref_umap.dib_idents.UMAP.pdf", legend.pos = "right", plot.width = 10)

customVln(features   = c("TREM2", "PLTP", "PLIN2", "OLR1",  "IL1B",  "CD14", 
                         "CD9",   "MRC1",  "TREM1", "ABCA1", "TNF",   "FCGR1A",
                         "GPNMB", "FOLR2", "HSPA6", "ABCG1", "CASP1", "FCGR3A"), 
          object     = dib.seurat, 
          group.by   = "predicted.id",
          ncol       = 6,
          name       = "Dib_integration/dib.full_Ref LAM - resident - foam - foam - inflammatory - general markers Violin.pdf",
          draw.names = F,
          cols       = M.int_refined.pop.colors,
          width      = 30,
          height     = 15
)

customFeature(features   = c("TREM2", "PLTP", "PLIN2", "OLR1",  "IL1B",  "CD14", 
                             "CD9",   "MRC1",  "TREM1", "ABCA1", "TNF",   "FCGR1A",
                             "GPNMB", "FOLR2", "HSPA6", "ABCG1", "CASP1", "FCGR3A"), 
              object     = dib.seurat, 
              ncol       = 6,
              pt.size    = 1.1,
              name       = "Dib_integration/dib.full_Ref LAM - resident - foam - foam - inflammatory - general markers Feature.pdf",
              width      = 30,
              height     = 15)

customVln(features   = c("TREM2", "PLTP", "PLIN2", "OLR1",  "IL1B",  "CD14", 
                         "CD9",   "MRC1",  "TREM1", "ABCA1", "TNF",   "FCGR1A",
                         "GPNMB", "FOLR2", "HSPA6", "ABCG1", "CASP1", "FCGR3A"), 
          object     = dib.seurat, 
          ncol       = 6,
          name       = "Dib_integration/dib_Ref LAM - resident - foam - foam - inflammatory - general markers Violin.pdf",
          draw.names = T,
          width      = 30,
          height     = 15
)

customFeature(features   = c("TREM2", "PLTP", "PLIN2", "OLR1",  "IL1B",  "CD14", 
                             "CD9",   "MRC1",  "TREM1", "ABCA1", "TNF",   "FCGR1A",
                             "GPNMB", "FOLR2", "HSPA6", "ABCG1", "CASP1", "FCGR3A"), 
          object     = dib.seurat, 
          ncol       = 6,
          pt.size    = 1.1,
          name       = "Dib_integration/dib_Ref LAM - resident - foam - foam - inflammatory - general markers Feature.pdf",
          width      = 30,
          height     = 15
)

customVln(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "HMOX1", "FOLR2", "C1QA", "HLA-DRB5", "CD1C", "LAMP3", "CLEC9A", "SIGLEC6"), 
          object     = dib.seurat, 
          ncol       = 5,
          name       = "Dib_integration/dib_Ref Dib markers Violin.pdf",
          draw.names = T,
          width      = 30,
          height     = 15
)

customFeature(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "HMOX1", "FOLR2", "C1QA", "HLA-DRB5", "CD1C", "LAMP3", "CLEC9A", "SIGLEC6"), 
          object     = dib.seurat, 
          ncol       = 5,
          pt.size    = 1.1,
          name       = "Dib_integration/dib_Ref Dib markers Feature.pdf",
          width      = 30,
          height     = 15
)

customVln(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "HMOX1", "FOLR2", "C1QA", "HLA-DRB5", "CD1C", "LAMP3", "CLEC9A", "SIGLEC6"), 
          object     = dib.seurat, 
          group.by   = "predicted.id",
          ncol       = 5,
          name       = "Dib_integration/dib.full_Ref Dib markers Violin.pdf",
          draw.names = F,
          cols       = M.int_refined.pop.colors,
          width      = 30,
          height     = 15
)

customFeature(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "HMOX1", "FOLR2", "C1QA", "HLA-DRB5", "CD1C", "LAMP3", "CLEC9A", "SIGLEC6"), 
          object     = dib.seurat, 
          ncol       = 5,
          pt.size    = 1.1,
          name       = "Dib_integration/dib.full_Ref Dib markers Feature.pdf",
          width      = 30,
          height     = 15
)

## Mac only
## Transfer our labels by reference mapping
# Find anchors
DefaultAssay(dib.seurat) <- "RNA"
DefaultAssay(final.pop.call.integrated.mye.velocyto.seurat) <- "integrated"

# Fix counts assay for integrated assay
cm <- final.pop.call.integrated.mye.velocyto.seurat[["RNA"]]$counts
ft <- row.names(final.pop.call.integrated.mye.velocyto.seurat)
cm <- cm[row.names(cm) %in% ft,]
final.pop.call.integrated.mye.velocyto.seurat[["integrated"]]$counts <- cm

# Fix scale data
final.pop.call.integrated.mye.velocyto.seurat <- ScaleData(final.pop.call.integrated.mye.velocyto.seurat)

dib.anchors <- FindTransferAnchors(reference           = final.pop.call.integrated.mye.velocyto.seurat,
                                   query               = dib.seurat,
                                   dims                = 1:30, 
                                   reference.reduction = "pca")

# Predict labels
predictions <- TransferData(anchorset = dib.anchors, 
                            refdata   = Idents(final.pop.call.integrated.mye.velocyto.seurat),
                            dims      = 1:30)

# Add to the Dib seurat object
mac.dib.seurat <- AddMetaData(dib.seurat, metadata = predictions)
# Map to ref UMAP
#final.pop.call.integrated.mye.velocyto.seurat <- RunUMAP(final.pop.call.integrated.mye.velocyto.seurat, dims = 1:30, reduction = "pca", return.model = TRUE)
#final.pop.call.integrated.mye.velocyto.seurat <- UpdateSeuratObject(final.pop.call.integrated.mye.velocyto.seurat)
mac.dib.seurat <- MapQuery(anchorset           = dib.anchors, 
                           reference           = final.pop.call.integrated.mye.velocyto.seurat, 
                           query               = mac.dib.seurat,
                           refdata             = Idents(final.pop.call.integrated.mye.velocyto.seurat),
                           reference.reduction = "pca", 
                           reduction.model     = "umap")

DefaultAssay(mac.dib.seurat) <- "RNA"

DimPlot(mac.dib.seurat, group.by = "predicted.id", reduction = "ref.umap", cols = full_set.colors) + NoLegend()
DimPlot(final.pop.call.integrated.mye.velocyto.seurat) + NoLegend()

# Visualize
customUMAP(object = final.pop.call.integrated.mye.velocyto.seurat, pt.size = 1.2, label = F, title = "Mac UMAP", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "Dib_integration/mac_ref.UMAP.pdf", cols = full_set.colors, plot.width = 15)
customUMAP(object = dib.seurat, pt.size = 1.2, label = F, title = "Dib et al. ref mapped", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "Dib_integration/dib.mac_ref_map.dib_idents.UMAP.pdf",  plot.width = 15)
customUMAP(object = dib.seurat, group.by = "predicted.id", pt.size = 1.2, label = F, title = "Dib et al. ref mapped", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "Dib_integration/dib.mac_ref_map.UMAP.pdf", cols = full_set.colors, plot.width = 25)
customUMAP(object = dib.seurat, group.by = "predicted.id", reduction = "ref.umap", pt.size = 1.2, label = F, title = "Dib et al. ref UMAP mapped",
           shuffle = T, seed = 666, file.name = "Dib_integration/dib.mac_ref_umap.UMAP.pdf", legend.pos = "right", cols = full_set.colors, plot.width = 15)
customUMAP(object = dib.seurat, reduction = "ref.umap", pt.size = 1.2, label = F, title = "Dib et al. ref UMAP mapped",
           shuffle = T, seed = 666, file.name = "Dib_integration/dib.mac_ref_umap.dib_idents.UMAP.pdf", legend.pos = "right", plot.width = 10)

customVln(features   = c("TREM2", "PLTP", "PLIN2", "OLR1",  "IL1B",  "CD14", 
                         "CD9",   "MRC1",  "TREM1", "ABCA1", "TNF",   "FCGR1A",
                         "GPNMB", "FOLR2", "HSPA6", "ABCG1", "CASP1", "FCGR3A"), 
          object     = dib.seurat, 
          group.by   = "predicted.id",
          ncol       = 6,
          name       = "Dib_integration/dib.mac_Ref LAM - resident - foam - foam - inflammatory - general markers Violin.pdf",
          draw.names = F,
          cols       = M.int_refined.pop.colors,
          width      = 30,
          height     = 15
)

customFeature(features   = c("TREM2", "PLTP", "PLIN2", "OLR1",  "IL1B",  "CD14", 
                             "CD9",   "MRC1",  "TREM1", "ABCA1", "TNF",   "FCGR1A",
                             "GPNMB", "FOLR2", "HSPA6", "ABCG1", "CASP1", "FCGR3A"), 
          object     = dib.seurat, 
          ncol       = 6,
          pt.size    = 1.1,
          name       = "Dib_integration/dib.mac_Ref LAM - resident - foam - foam - inflammatory - general markers Feature.pdf",
          width      = 30,
          height     = 15
)

customVln(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "HMOX1", "FOLR2", "C1QA", "HLA-DRB5", "CD1C", "LAMP3", "CLEC9A", "SIGLEC6"), 
          object     = dib.seurat, 
          group.by   = "predicted.id",
          ncol       = 5,
          name       = "Dib_integration/dib.mac_Ref Dib markers Violin.pdf",
          draw.names = F,
          cols       = M.int_refined.pop.colors,
          width      = 30,
          height     = 15
)

customFeature(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "HMOX1", "FOLR2", "C1QA", "HLA-DRB5", "CD1C", "LAMP3", "CLEC9A", "SIGLEC6"), 
          object     = dib.seurat, 
          ncol       = 5,
          pt.size    = 1.1,
          name       = "Dib_integration/dib.mac_Ref Dib markers Feature.pdf",
          width      = 30,
          height     = 15
)

## Remove DC clusters for cleaner visuals
mac.dib.seurat <- subset(dib.seurat, idents = c(0,1,3,4,5,6,8,9))
customUMAP(object = mac.dib.seurat, pt.size = 1.2, label = F, title = "Dib et al. data - no DCs", shuffle = T, seed = 666, plot.width = 10, file.name = "Dib_integration/mac_only.dib.UMAP.pdf")

## Transfer our labels by reference mapping
# Find anchors
DefaultAssay(mac.dib.seurat) <- "RNA"
DefaultAssay(final.pop.call.integrated.mye.velocyto.seurat) <- "integrated"
dib.anchors <- FindTransferAnchors(reference           = final.pop.call.integrated.mye.velocyto.seurat,
                                   query               = mac.dib.seurat,
                                   dims                = 1:30, 
                                   reference.reduction = "pca")

# Predict labels
predictions <- TransferData(anchorset = dib.anchors, 
                            refdata   = Idents(final.pop.call.integrated.mye.velocyto.seurat),
                            dims      = 1:30)

# Add to the Dib seurat object
mac.dib.seurat <- AddMetaData(mac.dib.seurat, metadata = predictions)

# Map to ref UMAP
final.pop.call.integrated.mye.velocyto.seurat  <- RunUMAP(final.pop.call.integrated.mye.velocyto.seurat, dims = 1:30, reduction = "pca", return.model = TRUE)
        mac.dib.seurat <- MapQuery(anchorset       = dib.anchors, 
                               reference           = final.pop.call.integrated.mye.velocyto.seurat, 
                               query               = mac.dib.seurat,
                               refdata             = Idents(final.pop.call.integrated.mye.velocyto.seurat), 
                               reference.reduction = "pca", 
                               reduction.model     = "umap")

# Visualize
customUMAP(object = final.pop.call.integrated.mye.velocyto.seurat, pt.size = 1.2, label = F, title = "Mac UMAP", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "Dib_integration/mac_only_dib.mac_ref.UMAP.pdf", cols = full_set.colors, plot.width = 15)
customUMAP(object = mac.dib.seurat, pt.size = 1.2, label = F, title = "Dib et al. ref mapped", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "Dib_integration/mac_only_dib.dib.mac_ref_map.dib_idents.UMAP.pdf",  plot.width = 15)
customUMAP(object = mac.dib.seurat, group.by = "predicted.id", pt.size = 1.5, label = F, title = "Dib et al. ref mapped", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "Dib_integration/mac_only_dib.dib.mac_ref_map.UMAP.pdf", cols = full_set.colors, plot.width = 20)
customUMAP(object = mac.dib.seurat, group.by = "predicted.id", reduction = "ref.umap", pt.size = 1.2, label = F, title = "Dib et al. ref UMAP mapped",
           shuffle = T, seed = 666, file.name = "Dib_integration/mac_only_dib.dib.mac_ref_umap.UMAP.pdf", legend.pos = "right", cols = full_set.colors, plot.width = 15)
customUMAP(object = mac.dib.seurat, reduction = "ref.umap", pt.size = 1.2, label = F, title = "Dib et al. ref UMAP mapped",
           shuffle = T, seed = 666, file.name = "Dib_integration/mac_only_dib.dib.mac_ref_umap.dib_idents.UMAP.pdf", legend.pos = "right", plot.width = 10)

customVln(features   = c("TREM2", "PLTP", "PLIN2", "OLR1",  "IL1B",  "CD14", 
                         "CD9",   "MRC1",  "TREM1", "ABCA1", "TNF",   "FCGR1A",
                         "GPNMB", "FOLR2", "HSPA6", "ABCG1", "CASP1", "FCGR3A"), 
          object     = mac.dib.seurat, 
          group.by   = "predicted.id",
          ncol       = 6,
          name       = "Dib_integration/mac_only_dib.dib.mac_Ref LAM - resident - foam - foam - inflammatory - general markers Violin.pdf",
          draw.names = F,
          cols       = M.int_refined.pop.colors,
          width      = 30,
          height     = 15
)

customFeature(features   = c("TREM2", "PLTP", "PLIN2", "OLR1",  "IL1B",  "CD14", 
                             "CD9",   "MRC1",  "TREM1", "ABCA1", "TNF",   "FCGR1A",
                             "GPNMB", "FOLR2", "HSPA6", "ABCG1", "CASP1", "FCGR3A"), 
              object     = mac.dib.seurat, 
              ncol       = 6,
              pt.size    = 1.1,
              name       = "Dib_integration/mac_only_dib.dib.mac_Ref LAM - resident - foam - foam - inflammatory - general markers Feature.pdf",
              width      = 30,
              height     = 15
)

customVln(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "HMOX1", "FOLR2", "C1QA", "HLA-DRB5", "CD1C", "LAMP3", "CLEC9A", "SIGLEC6"), 
          object     = mac.dib.seurat, 
          group.by   = "predicted.id",
          ncol       = 5,
          name       = "Dib_integration/mac_only_dib.dib.mac_Ref Dib markers Violin.pdf",
          draw.names = F,
          cols       = M.int_refined.pop.colors,
          width      = 30,
          height     = 15
)

customVln(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "HMOX1", "FOLR2", "C1QA", "HLA-DRB5", "MRC1"), 
          object     = mac.dib.seurat, 
          group.by   = "predicted.id",
          ncol       = 3,
          name       = "Dib_integration/mac_only_dib.dib.mac_Ref Dib markers curated Violin.pdf",
          draw.names = F,
          cols       = M.int_refined.pop.colors,
          width      = 15,
          height     = 15
)

customFeature(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "HMOX1", "FOLR2", "C1QA", "HLA-DRB5", "CD1C", "LAMP3", "CLEC9A", "SIGLEC6"), 
              object     = mac.dib.seurat, 
              ncol       = 5,
              pt.size    = 1.1,
              name       = "Dib_integration/mac_only_dib.dib.mac_Ref Dib markers Feature.pdf",
              width      = 30,
              height     = 15
)

customVln(features   = c("TREM1", "PLIN2", "IL1B", "TNF", "HSPA6", "CASP1"), 
          object     = mac.dib.seurat, 
          group.by   = "predicted.id",
          cols       = M.int_refined.pop.colors,
          ncol       = 2,
          name       = "Dib_integration/ddata plabels foam dif markers Violin.pdf",
          draw.names = F,
          width      = 15,
          height     = 20
)

customVln(features   = c("PLTP", "C1QA", "FOLR2", "GPNMB", "TREM2", "CD9", "IL1B", "TNF"), 
          object     = mac.dib.seurat, 
          group.by   = "predicted.id",
          cols       = M.int_refined.pop.colors,
          ncol       = 2,
          name       = "Dib_integration/ddata plabels res dif markers Violin.pdf",
          draw.names = F,
          width      = 15,
          height     = 20
)

customVln(features   = c("S100A8", "SELL", "FCGR3A", "MX1"), 
          object     = mac.dib.seurat, 
          group.by   = "predicted.id",
          cols       = M.int_refined.pop.colors,
          ncol       = 2,
          name       = "Dib_integration/ddata plabels inf dif markers Violin.pdf",
          draw.names = F,
          width      = 15,
          height     = 10
)


customVln(features   = c("TREM1", "PLIN2", "IL1B", "TNF", "HSPA6", "CASP1"), 
          object     = mac.dib.seurat, 
          
          ncol       = 2,
          name       = "Dib_integration/ddata dlabels foam dif markers Violin.pdf",
          draw.names = T,
          width      = 15,
          height     = 20
)

customVln(features   = c("PLTP", "C1QA", "FOLR2", "GPNMB", "TREM2", "CD9", "IL1B", "TNF"), 
          object     = mac.dib.seurat, 
          ncol       = 2,
          name       = "Dib_integration/ddata dlabels res dif markers Violin.pdf",
          draw.names = T,
          width      = 15,
          height     = 20
)

customVln(features   = c("S100A8", "SELL", "FCGR3A", "MX1"), 
          object     = mac.dib.seurat, 
          ncol       = 2,
          name       = "Dib_integration/ddata dlabels inf dif markers Violin.pdf",
          draw.names = T,
          width      = 15,
          height     = 10
)


## Vice Versa -- map our stuff to the Dib types and UMAP
## Use the mac only object!
# Save the Dib idents by appending a 'D'
mac.dib.seurat         <- AddMetaData(mac.dib.seurat, paste0("D", Idents(mac.dib.seurat)), col.name = "dib.idents")
Idents(mac.dib.seurat) <- mac.dib.seurat$dib.idents

## Transfer our labels by reference mapping
# Find anchors
DefaultAssay(mac.dib.seurat) <- "RNA"
mac.dib.seurat <- FindVariableFeatures(mac.dib.seurat)
DefaultAssay(final.pop.call.integrated.mye.velocyto.seurat) <- "integrated"
mac.anchors <- FindTransferAnchors(reference           = mac.dib.seurat,
                                   query               = final.pop.call.integrated.mye.velocyto.seurat,
                                   dims                = 1:30, 
                                   reference.reduction = "pca")

# Predict labels
predictions <- TransferData(anchorset = mac.anchors, 
                            refdata   = Idents(mac.dib.seurat),
                            dims      = 1:30)

# Add to the seurat object
final.pop.call.integrated.mye.velocyto.seurat <- AddMetaData(final.pop.call.integrated.mye.velocyto.seurat, metadata = predictions)

# Map to ref UMAP
mac.dib.seurat <- RunUMAP(mac.dib.seurat, dims = 1:30, reduction = "pca", return.model = TRUE)
final.pop.call.integrated.mye.velocyto.seurat <- MapQuery(anchorset           = mac.anchors, 
                       reference           = mac.dib.seurat, 
                       query               = final.pop.call.integrated.mye.velocyto.seurat,
                       refdata             = Idents(mac.dib.seurat), 
                       reference.reduction = "pca", 
                       reduction.model     = "umap")

# Visualize
customUMAP(object = final.pop.call.integrated.mye.velocyto.seurat, group.by = "predicted.id", pt.size = 1.5, label = F, title = "Dib idents", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "Dib_integration/mac.dib_ref_map.UMAP.dib_idents.pdf", plot.width = 10)

customVln(features   = c("TREM2", "PLTP", "PLIN2", "OLR1",  "IL1B",  "CD14", 
                         "CD9",   "MRC1",  "TREM1", "ABCA1", "TNF",   "FCGR1A",
                         "GPNMB", "FOLR2", "HSPA6", "ABCG1", "CASP1", "FCGR3A"), 
                   object     = final.pop.call.integrated.mye.velocyto.seurat, 
                   group.by   = "predicted.id",
                   ncol       = 6,
                   name       = "Dib_integration/mac.dib_Ref LAM - resident - foam - foam - inflammatory - general markers Violin.pdf",
                   draw.names = T,
                   width      = 30,
                   height     = 15
)

customFeature(features   = c("TREM2", "PLTP", "PLIN2", "OLR1",  "IL1B",  "CD14", 
                             "CD9",   "MRC1",  "TREM1", "ABCA1", "TNF",   "FCGR1A",
                             "GPNMB", "FOLR2", "HSPA6", "ABCG1", "CASP1", "FCGR3A"), 
          object     = final.pop.call.integrated.mye.velocyto.seurat, 
          ncol       = 6,
          pt.size    = 1.1,
          name       = "Dib_integration/mac.dib_Ref LAM - resident - foam - foam - inflammatory - general markers Feature.pdf",
          width      = 30,
          height     = 15
)

customVln(features   = c("TREM2", "PLTP", "PLIN2", "OLR1",  "IL1B",  "CD14", 
                         "CD9",   "MRC1",  "TREM1", "ABCA1", "TNF",   "FCGR1A",
                         "GPNMB", "FOLR2", "HSPA6", "ABCG1", "CASP1", "FCGR3A"), 
          object     = final.pop.call.integrated.mye.velocyto.seurat, 
          ncol       = 6,
          cols       = M.int_refined.pop.colors,
          name       = "Dib_integration/mac_Ref LAM - resident - foam - foam - inflammatory - general markers Violin.pdf",
          draw.names = F,
          width      = 30,
          height     = 15
)

customFeature(features   = c("TREM2", "PLTP", "PLIN2", "OLR1",  "IL1B",  "CD14", 
                             "CD9",   "MRC1",  "TREM1", "ABCA1", "TNF",   "FCGR1A",
                             "GPNMB", "FOLR2", "HSPA6", "ABCG1", "CASP1", "FCGR3A"), 
          object     = final.pop.call.integrated.mye.velocyto.seurat, 
          ncol       = 6,
          pt.size    = 1.1,
          name       = "Dib_integration/mac_Ref LAM - resident - foam - foam - inflammatory - general markers Feature.pdf",
          width      = 30,
          height     = 15
)

customVln(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "HMOX1", "FOLR2", "C1QA", "HLA-DRB5", "CD1C", "LAMP3", "CLEC9A", "SIGLEC6"), 
          object     = final.pop.call.integrated.mye.velocyto.seurat, 
          group.by   = "predicted.id",
          ncol       = 5,
          name       = "Dib_integration/mac.dib_Ref Dib markers Violin.pdf",
          draw.names = T,
          width      = 30,
          height     = 15
)

customFeature(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "HMOX1", "FOLR2", "C1QA", "HLA-DRB5", "CD1C", "LAMP3", "CLEC9A", "SIGLEC6"), 
          object     = final.pop.call.integrated.mye.velocyto.seurat, 
          ncol       = 5, 
          pt.size    = 1.1,
          name       = "Dib_integration/mac.dib_Ref Dib markers Feature.pdf",
          width      = 30,
          height     = 15
)

customVln(features   = c("TREM1", "PLIN2", "IL1B", "TNF", "HSPA6", "CASP1"), 
          object     = final.pop.call.integrated.mye.velocyto.seurat, 
          group.by   = "predicted.id",
          ncol       = 2,
          name       = "Dib_integration/pdata dlabels foam dif markers Violin.pdf",
          draw.names = T,
          width      = 15,
          height     = 20
)

customVln(features   = c("PLTP", "C1QA", "FOLR2", "GPNMB", "TREM2", "CD9", "IL1B", "TNF"), 
          object     = final.pop.call.integrated.mye.velocyto.seurat, 
          group.by   = "predicted.id",
          ncol       = 2,
          name       = "Dib_integration/pdata dlabels res dif markers Violin.pdf",
          draw.names = T,
          width      = 15,
          height     = 20
)

customVln(features   = c("S100A8", "SELL", "FCGR3A", "MX1"), 
          object     = final.pop.call.integrated.mye.velocyto.seurat, 
          group.by   = "predicted.id",
          ncol       = 2,
          name       = "Dib_integration/pdata dlabels inf dif markers Violin.pdf",
          draw.names = T,
          width      = 15,
          height     = 10
)


customVln(features   = c("TREM1", "PLIN2", "IL1B", "TNF", "HSPA6", "CASP1"), 
          object     = final.pop.call.integrated.mye.velocyto.seurat, 
          cols       = M.int_refined.pop.colors,
          ncol       = 2,
          name       = "Dib_integration/pdata plabels foam dif markers Violin.pdf",
          draw.names = F,
          width      = 15,
          height     = 20
)

customVln(features   = c("PLTP", "C1QA", "FOLR2", "GPNMB", "TREM2", "CD9", "IL1B", "TNF"), 
          object     = final.pop.call.integrated.mye.velocyto.seurat, 
          cols       = M.int_refined.pop.colors,
          ncol       = 2,
          name       = "Dib_integration/pdata plabels res dif markers Violin.pdf",
          draw.names = F,
          width      = 15,
          height     = 20
)

customVln(features   = c("S100A8", "SELL", "FCGR3A", "MX1"), 
          object     = final.pop.call.integrated.mye.velocyto.seurat, 
          cols       = M.int_refined.pop.colors,
          ncol       = 2,
          name       = "Dib_integration/pdata plabels inf dif markers Violin.pdf",
          draw.names = F,
          width      = 15,
          height     = 10
)


## Recluster the Dib TREM1 macs
dir.create("Dib_integration/TREM1_cluster", showWarnings = F, recursive = T)
# Subset the TREM1 cluster
customUMAP(object = dib.seurat, file.name = "Dib_integration/TREM1_cluster/TREM1 cluster.pdf", cells.highlight = WhichCells(dib.seurat, idents = 4))
trem1.dib.seurat <- subset(dib.seurat, idents = 4)
trem1.dib.seurat

## Recluster
DefaultAssay(trem1.dib.seurat) <- "RNA"

# Run PCA and UMAP embedding
trem1.dib.seurat <- RunPCA(trem1.dib.seurat, verbose = T)
trem1.dib.seurat <- RunUMAP(trem1.dib.seurat, dims = 1:30, verbose = T)

# Transform
trem1.dib.seurat <- SCTransform(trem1.dib.seurat)

# Clustering
set.seed(1)
trem1.dib.seurat <- FindNeighbors(trem1.dib.seurat, dims = 1:30, verbose = T, force.recalc = T)
trem1.dib.seurat <- FindClusters(trem1.dib.seurat, verbose = T, resolution = 0.3, random.seed = 666)

# Visualise
DimPlot(trem1.dib.seurat, label = TRUE, pt.size = 3, label.size = 10, shuffle = T) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

customUMAP(object = trem1.dib.seurat, file.name = "Dib_integration/TREM1_cluster/TREM1 subclustered UMAP.pdf")

customVln(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "FOLR2", "C1QA"), 
          object     = trem1.dib.seurat, 
          ncol       = 3,
          name       = "Dib_integration/TREM1_cluster/markers Violin.pdf",
          draw.names = T,
          width      = 15,
          height     = 15
)

customFeature(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "FOLR2", "C1QA"), 
              object     = trem1.dib.seurat, 
              ncol       = 3, 
              pt.size    = 1.1,
              name       = "Dib_integration/TREM1_cluster/markers Feature.pdf",
              width      = 15,
              height     = 15
)


customVln(features   = c("TNF", "HSPA6", "ABCA1"), 
          object     = trem1.dib.seurat, 
          ncol       = 3,
          name       = "Dib_integration/TREM1_cluster/markers 2 Violin.pdf",
          draw.names = T,
          width      = 15,
          height     = 5
)

customFeature(features   = c("TNF", "HSPA6", "ABCA1"), 
              object     = trem1.dib.seurat, 
              ncol       = 3, 
              pt.size    = 1.1,
              name       = "Dib_integration/TREM1_cluster/markers 2 Feature.pdf",
              width      = 10,
              height     = 10
)


## Add our labels
pid              <- subset(mac.dib.seurat, idents = "D4")$predicted.id
pid              <- pid[row.names(trem1.dib.seurat@meta.data)]
unique(pid)
trem1.dib.seurat <- AddMetaData(trem1.dib.seurat, pid, col.name = "predicted.id")

customVln(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "FOLR2", "C1QA"), 
          object     = trem1.dib.seurat, 
          group.by   = "predicted.id",
          cols       = M.int_refined.pop.colors,
          ncol       = 3,
          name       = "Dib_integration/TREM1_cluster/Prange_labels markers Violin.pdf",
          draw.names = F,
          width      = 20,
          height     = 10
)

customUMAP(object     = trem1.dib.seurat, 
           group.by   = "predicted.id", 
           file.name  = "Dib_integration/TREM1_cluster/TREM1 subclustered Prange labels UMAP.pdf", 
           legend.pos = "right", 
           cols       = M.int_refined.pop.colors, 
           plot.width = 15
)

## Check DE markers
DefaultAssay(trem1.dib.seurat) <- "RNA"

# Normalize and scale
trem1.dib.seurat <- NormalizeData(trem1.dib.seurat, normalization.method = "LogNormalize")
trem1.dib.seurat <- ScaleData(trem1.dib.seurat, features = row.names(trem1.dib.seurat))

# Find markers
trem1.dib.markers <- FindAllMarkers(trem1.dib.seurat, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

# Save the top9 markers per cluster
sep.markers <- trem1.dib.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
num.clus <- length(unique(Idents(trem1.dib.seurat)))
for(i in levels(Idents(trem1.dib.seurat))){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = trem1.dib.seurat,
                     name      = paste("Dib_integration/TREM1_cluster/Cluster ",i, " markers", sep = ""),
                     
                     Vln.width = 25, Vln.height = 25)
}

# Save markers per cluster to disk
for(i in levels(Idents(trem1.dib.seurat))){
  x <- trem1.dib.markers[which(trem1.dib.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("Dib_integration/TREM1_cluster/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

## Check how our foam cells look if we subset them out
unique(Idents(final.pop.call.integrated.mye.velocyto.seurat))
foam.idents               <- unique(Idents(final.pop.call.integrated.mye.velocyto.seurat))[c(3,5,10)]
foam.integrated.my.seurat <- subset(final.pop.call.integrated.mye.velocyto.seurat, idents = foam.idents)
foam.integrated.my.seurat

## Recluster
DefaultAssay(trem1.dib.seurat) <- "RNA"

# Store idents
foam.integrated.my.seurat <- AddMetaData(foam.integrated.my.seurat, Idents(foam.integrated.my.seurat), col.name = "foam.idents")

# Run PCA and UMAP embedding
foam.integrated.my.seurat <- RunPCA(foam.integrated.my.seurat, verbose = T)
foam.integrated.my.seurat <- RunUMAP(foam.integrated.my.seurat, dims = 1:30, verbose = T)

# Transform
foam.integrated.my.seurat <- SCTransform(foam.integrated.my.seurat)

# Clustering
set.seed(1)
foam.integrated.my.seurat <- FindNeighbors(foam.integrated.my.seurat, dims = 1:30, verbose = T, force.recalc = T)
names(foam.integrated.my.seurat@graphs)
foam.integrated.my.seurat <- FindClusters(foam.integrated.my.seurat, verbose = T, resolution = 0.5, random.seed = 666, graph.name = "integrated_nn")

# Visualise
DimPlot(foam.integrated.my.seurat, label = TRUE, pt.size = 3, label.size = 10, shuffle = T) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

customUMAP(object = foam.integrated.my.seurat, file.name = "Dib_integration/TREM1_cluster/Prange Foam subclustered UMAP.pdf", legend.pos = "right")
customUMAP(object = foam.integrated.my.seurat, group.by = "foam.idents", file.name = "Dib_integration/TREM1_cluster/Prange Foam subclustered original label UMAP.pdf", legend.pos = "right", plot.width = 15)

customVln(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "FOLR2", "C1QA"), 
          object     = foam.integrated.my.seurat, 
          ncol       = 3,
          name       = "Dib_integration/TREM1_cluster/Prange Foam subclustered markers Violin.pdf",
          draw.names = T,
          width      = 15,
          height     = 15
)

customFeature(features   = c("S100A8", "IL10", "IL1B", "TREM1", "PLIN2", "TREM2", "CD9", "FOLR2", "C1QA"), 
              object     = foam.integrated.my.seurat, 
              ncol       = 3, 
              pt.size    = 1.1,
              name       = "Dib_integration/TREM1_cluster/Prange Foam subclustered markers Feature.pdf",
              width      = 15,
              height     = 15
)


customVln(features   = c("TNF", "HSPA6", "ABCA1"), 
          object     = foam.integrated.my.seurat, 
          ncol       = 3,
          name       = "Dib_integration/TREM1_cluster/Prange Foam subclustered markers 2 Violin.pdf",
          draw.names = T,
          width      = 15,
          height     = 5
)

customFeature(features   = c("TNF", "HSPA6", "ABCA1"), 
              object     = foam.integrated.my.seurat, 
              ncol       = 3, 
              pt.size    = 1.1,
              name       = "Dib_integration/TREM1_cluster/Prange Foam subclustered markers 2 Feature.pdf",
              width      = 10,
              height     = 10
)

## Check DE markers
DefaultAssay(foam.integrated.my.seurat) <- "RNA"

# Normalize and scale
foam.integrated.my.seurat <- NormalizeData(foam.integrated.my.seurat, normalization.method = "LogNormalize")
foam.integrated.my.seurat <- ScaleData(foam.integrated.my.seurat, features = row.names(foam.integrated.my.seurat))

# Find markers
foam.integrated.my.markers <- FindAllMarkers(foam.integrated.my.seurat, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

# Save the top9 markers per cluster
sep.markers <- foam.integrated.my.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
num.clus <- length(unique(Idents(foam.integrated.my.seurat)))
for(i in levels(Idents(foam.integrated.my.seurat))){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = foam.integrated.my.seurat,
                     name      = paste("Dib_integration/TREM1_cluster/Prange Foam subclustered Cluster ",i, " markers", sep = ""),
                     
                     Vln.width = 25, Vln.height = 25)
}

# Save markers per cluster to disk
for(i in levels(Idents(foam.integrated.my.seurat))){
  x <- foam.integrated.my.markers[which(foam.integrated.my.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  write.table(x, file = paste("Dib_integration/TREM1_cluster/Prange Foam subclustered Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

