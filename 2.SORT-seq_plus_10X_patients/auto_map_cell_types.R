#=========================================================================================================================
## Map to RNA + CITE reference model
## Ref model is PBMC so while this does give a rough overview, it will miss any plaque cell types that are not present in the PBMC reference
##========================================================================================================================
dir.create("reference_mapping_results", showWarnings = F)

# Load reference data (162K PBMCs with 228 antibodies, from Sajita lab)
reference <- LoadH5Seurat("Seurat_Objects/pbmc_multimodal.h5seurat")
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("reference_mapping_results/reference UMAP.pdf", width = 10, height = 10)

# Find anchors
anchors <- FindTransferAnchors(reference            = reference,
                               query                = patients.seurat,
                               normalization.method = "SCT",
                               reference.reduction  = "pca",
                                               dims = 1:50
)

# Map to reference
patients.seurat <- MapQuery(anchorset               = anchors,
                               query                = patients.seurat,
                               reference            = reference,
                               refdata              = list(celltype.l1   = "celltype.l1",
                                                           celltype.l2   = "celltype.l2",
                                                           predicted_ADT = "ADT"),
                                reference.reduction = "pca", 
                                reduction.model     = "wnn.umap"
)

#=========================================================================================================================
## Visualise projections
##========================================================================================================================
## Visualise projections onto the pre-exisitng reference UMAP (reference on-existent cells (e.g. macrophages) will be mapped towards the closest match)
DefaultAssay(patients.seurat) <- "integrated"
p1 = DimPlot(patients.seurat, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE, shuffle = T) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p2 = DimPlot(patients.seurat, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE, shuffle = T) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p1 + p2
ggsave("reference_mapping_results/Reference mapped reference UMAP.pdf", width = 20, height = 10)

# Also plot on th equery (excisintg UMAP)
p1 = DimPlot(patients.seurat, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE, shuffle = T) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p2 = DimPlot(patients.seurat, reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE, shuffle = T) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p1 + p2
ggsave("reference_mapping_results/Reference mapped query UMAP.pdf", width = 20, height = 10)

## Visualize ADT and predicted ADT from the reference
## First what we actually measured
DefaultAssay(patients.seurat) <- 'ADT'
rownames(patients.seurat)

# Ref mapped
FeaturePlot(patients.seurat, features = rownames(patients.seurat), reduction = "ref.umap", cols = c("lightgrey", "darkgreen"), ncol = 2, order = T) +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("reference_mapping_results/Reference mapped reference UMAP Antibody binding CD3 CD4 CD8 CD14.pdf")

# Query mapped
FeaturePlot(patients.seurat, features = rownames(patients.seurat), reduction = "umap", cols = c("lightgrey", "darkgreen"), ncol = 2, order = T) +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("reference_mapping_results/Reference mapped qurey UMAP Antibody binding CD3 CD4 CD8 CD14.pdf")

## Then the prediction by ref mapping
DefaultAssay(patients.seurat) <- 'predicted_ADT'
rownames(patients.seurat)

# Ref mapped Ab panel
FeaturePlot(patients.seurat, features = c("CD3-1", "CD4-1", "CD8", "CD14"), reduction = "ref.umap", cols = c("lightgrey", "darkgreen"), ncol = 2, order = T) +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("reference_mapping_results/Reference mapped reference UMAP predicted Antibody binding CD3 CD4 CD8 CD14.pdf")

# Query mapped Ab panel
FeaturePlot(patients.seurat, features = c("CD3-1", "CD4-1", "CD8", "CD14"), reduction = "umap", cols = c("lightgrey", "darkgreen"), ncol = 2, order = T) +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("reference_mapping_results/Reference mapped query UMAP predicted Antibody binding CD3 CD4 CD8 CD14.pdf")

# Ref mapped expanded panel
FeaturePlot(patients.seurat, features = c("CD3-1", "CD4-1", "CD8", "CD14", "CD34", "CD45RA", "CD1c", "CD16", "CD79a", "CD68", "CD86", "CD40"), reduction = "ref.umap", cols = c("lightgrey", "darkgreen"), ncol = 3, order = T) +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("reference_mapping_results/Reference mapped reference UMAP predicted Antibody binding bigger panel.pdf", width = 20, height = 20)

# Query mapped expanded panel
FeaturePlot(patients.seurat, features = c("CD3-1", "CD4-1", "CD8", "CD14", "CD34", "CD45RA", "CD1c", "CD16", "CD79a", "CD68", "CD86", "CD40"), reduction = "umap", cols = c("lightgrey", "darkgreen"), ncol = 3, order = T) +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("reference_mapping_results/Reference mapped query UMAP predicted Antibody binding bigger panel.pdf", width = 20, height = 20)

# Reset default assay
DefaultAssay(patients.seurat) <- "integrated"

## Create a de novo visualisation instead, to better capture all cell types not available in the reference set
# Merge reference and query
reference$id          <- 'reference'
patients.seurat$id     <- 'query'
refquery              <- merge(reference, patients.seurat)
refquery[["pca"]]     <- merge(reference[["pca"]], patients.seurat[["ref.pca"]])
refquery              <- RunUMAP(refquery, reduction = 'pca', dims = 1:50)

# Visualise
p1 <- DimPlot(refquery, group.by = 'id', shuffle = TRUE) +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p2 <- DimPlot(refquery, group.by = 'Tissue', shuffle = TRUE) +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p3 <- DimPlot(refquery, group.by = 'Patient', shuffle = TRUE) +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p1 + p2 + p3
ggsave("reference_mapping_results/Reference mapped de novo UMAP.pdf", width = 20, height = 10)
DimPlot(refquery, group.by = c('celltype.l2', 'predicted.celltype.l2'), shuffle = TRUE) +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("reference_mapping_results/Reference mapped de novo UMAP cell types.pdf", width = 20, height = 10)

#=========================================================================================================================
## Clean up
##========================================================================================================================
rm(list = ls()[!ls() %in% c("patients.seurat", "vars.to.regress")])
source("code/functions.R")
patients.seurat@meta.data$id <- NULL
save.image()

## Save the integrated seurat object
saveRDS(patients.seurat, "Seurat_Objects/main.patients.seurat.filtered.mapped.RDS")
#patients.seurat <- readRDS(file = "Seurat_Objects/main.patients.seurat.filtered.mapped.RDS")

