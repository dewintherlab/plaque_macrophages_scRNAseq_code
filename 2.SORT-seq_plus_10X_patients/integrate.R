#=========================================================================================================================
## Integrate SCTransformed, filtered samples
##========================================================================================================================
dir.create("integration_QC", showWarnings = F)

# Prep for integration
features        <- SelectIntegrationFeatures(object.list         = patients.seurat, nfeatures = 3000)
patients.seurat <- PrepSCTIntegration(       object.list         = patients.seurat, anchor.features = features)
anchors         <- FindIntegrationAnchors(   object.list         = patients.seurat, 
                                            normalization.method = "SCT",
                                            anchor.features      = features, 
                                            dims                 = 1:30, 
                                            reduction            = "rpca", 
                                            k.anchor             = 20)

# Integrate the samples
patients.seurat <- IntegrateData(anchorset           = anchors, 
                                normalization.method = "SCT", 
                                dims                 = 1:30)
# Show stats
patients.seurat

# Clean up to save memory
rm(list = ls()[!ls() %in% c("patients.seurat", "vars.to.regress")])
source("code/functions.R")

#=========================================================================================================================
## Run the embeddings, cluster, and have a look
##=========================================================================================================================
# Run dimensionality reductions
patients.seurat <- RunPCA( object = patients.seurat, verbose = FALSE)
patients.seurat <- RunUMAP(object = patients.seurat, reduction = "pca", dims = 1:30)

# Cluster
set.seed(1)
patients.seurat <- FindNeighbors(object = patients.seurat, reduction = "pca", dims = 1:30, force.recalc = T, verbose = T)
patients.seurat <- FindClusters( object = patients.seurat, resolution = 1.5, verbose = T, random.seed = 666)

# Visualize
DimPlot(patients.seurat, reduction = "umap", group.by = "Patient", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("integration_QC/Integration patient UMAP.pdf", width = 10, height = 10)

DimPlot(patients.seurat, reduction = "umap", group.by = "Tissue", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("integration_QC/Integration tissue UMAP.pdf", width = 10, height = 10)

DimPlot(patients.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("integration_QC/Integration cluster UMAP.pdf", width = 10, height = 10)

p1 <- DimPlot(patients.seurat, reduction = "umap", cells.highlight = list(P1 = WhichCells(patients.seurat, expression = Patient == "P1")), shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p2 <- DimPlot(patients.seurat, reduction = "umap", cells.highlight = list(P2 = WhichCells(patients.seurat, expression = Patient == "P2")), shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p3 <- DimPlot(patients.seurat, reduction = "umap", cells.highlight = list(P3 = WhichCells(patients.seurat, expression = Patient == "P3")), shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p1 + p2 + p3
ggsave("integration_QC/Integration patients separate UMAP.pdf", width = 30, height = 10)

p1 <- DimPlot(patients.seurat, reduction = "umap", group.by = "Patient", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "bottom") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p2 <- DimPlot(patients.seurat, reduction = "umap", group.by = "Tissue", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "bottom") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p3 <- DimPlot(patients.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p1 + p2 + p3
ggsave("integration_QC/Integration UMAP.pdf", width = 30, height = 10)


#=========================================================================================================================
## Have a look at the antibody capture data
##========================================================================================================================
DefaultAssay(patients.seurat) <- "ADT"

# Normalise
patients.seurat <- NormalizeData(object = patients.seurat, normalization.method = "CLR", margin = 2, verbose = T)

## Visualise
# Plot all
FeaturePlot(patients.seurat, row.names(patients.seurat@assays$ADT), pt.size = 1, cols = c("lightgrey", "darkgreen"), order = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("integration_QC/Antibody capture feature plots.pdf", width = 15, height = 15)

## Plot antibody overlaps
# CD3 vs CD4
FeaturePlot(patients.seurat, row.names(patients.seurat@assays$ADT)[c(1,2)], pt.size = 1, blend = T, order = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("integration_QC/Antibody capture CD3 vs CD4.pdf", width = 20, height = 5)

# CD3 vs CD8
FeaturePlot(patients.seurat, row.names(patients.seurat@assays$ADT)[c(1,3)], pt.size = 1, blend = T, order = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("integration_QC/Antibody capture CD3 vs CD8.pdf", width = 20, height = 5)

# CD4 vs CD8
FeaturePlot(patients.seurat, row.names(patients.seurat@assays$ADT)[c(2,3)], pt.size = 1, blend = T, order = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("integration_QC/Antibody capture CD4 vs CD8.pdf", width = 20, height = 5)

# CD3 vs CD14
FeaturePlot(patients.seurat, row.names(patients.seurat@assays$ADT)[c(1,4)], pt.size = 1, blend = T, order = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("integration_QC/Antibody capture CD3 vs CD14.pdf", width = 20, height = 5)

# CD4 vs CD14
FeaturePlot(patients.seurat, row.names(patients.seurat@assays$ADT)[c(2,4)], pt.size = 1, blend = T, order = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("integration_QC/Antibody capture CD4 vs CD14.pdf", width = 20, height = 5)

# CD8 vs CD14
FeaturePlot(patients.seurat, row.names(patients.seurat@assays$ADT)[c(3,4)], pt.size = 1, blend = T, order = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("integration_QC/Antibody capture CD8 vs CD14.pdf", width = 20, height = 5)

# Reset the default assay
DefaultAssay(patients.seurat) <- "integrated"

#=========================================================================================================================
## Clean up
##========================================================================================================================
rm(list = ls()[!ls() %in% c("patients.seurat", "vars.to.regress")])
source("code/functions.R")
save.image()

## Save the integrated seurat objects
saveRDS(patients.seurat, "Seurat_Objects/main.patients.seurat.integrated.RDS")
#patients.seurat <- readRDS(file = "Seurat_Objects/main.patients.seurat.integrated.RDS")


