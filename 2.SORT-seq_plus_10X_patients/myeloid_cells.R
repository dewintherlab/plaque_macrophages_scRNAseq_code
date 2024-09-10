#=========================================================================================================================
## Analyze the myeloid cells
##========================================================================================================================
dir.create("myeloid_results",            showWarnings = F)
dir.create("myeloid_results/clusters",   showWarnings = F)

##==========================================================================================================
## First let's have a closer look at some cell type defining genes
customVln(object = patients.seurat, features = c("CD1C", "CD3E", "CD4", "CD8A", "CD14", "FCGR3A", "CD68", "CD79B", "CD86", "NCAM1"), stack = T, name = "myeloid_results/All clusters - cell type defining genes.pdf", ncol = 3)

## Subset and analyze the myeloid cell populations
# Check cell types
unique(patients.seurat$predicted.celltype.l1)

# Extract meyloid cells
mye.cells <- row.names(patients.seurat@meta.data[grep("Mono|DC", patients.seurat$predicted.celltype.l2),])

# Subset the seurat object
mye.patients.seurat <- subset(patients.seurat, cells = mye.cells)

# Remove trace contaminants
DefaultAssay(mye.patients.seurat) <- "RNA"
mye.patients.seurat <- subset(mye.patients.seurat, subset = CD3E == 0 & CD8A == 0 & NCAM1 == 0 & CD79A == 0)

## Recalibrate clustering
# Work with the integrated assay
DefaultAssay(mye.patients.seurat) <- "integrated"

# Run PCA and UMAP embedding
mye.patients.seurat <- RunPCA(mye.patients.seurat, verbose = T)
mye.patients.seurat <- RunUMAP(mye.patients.seurat, dims = 1:30, verbose = T)

## Clustering
set.seed(1)
mye.patients.seurat <- FindNeighbors(mye.patients.seurat, dims = 1:30, verbose = T, force.recalc = T)
mye.patients.seurat <- FindClusters(mye.patients.seurat, verbose = T, resolution = 1, random.seed = 666)
num.clus          <- length(unique(Idents(mye.patients.seurat)))
DimPlot(mye.patients.seurat, label = TRUE, pt.size = 2, label.size = 10, shuffle = T) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("myeloid_results/UMAP clusters.pdf", width = 15, height = 15)

DimPlot(mye.patients.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Patient", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("myeloid_results/UMAP Patient.pdf", width = 15, height = 15)

DimPlot(mye.patients.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Patient", shuffle = T, cells.highlight = list(P1 = WhichCells(mye.patients.seurat, expression = Patient == "P1")) ) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("myeloid_results/UMAP Patient 1.pdf", width = 15, height = 15)

DimPlot(mye.patients.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Patient", shuffle = T, cells.highlight = list(P2 = WhichCells(mye.patients.seurat, expression = Patient == "P2")) ) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("myeloid_results/UMAP Patient 2.pdf", width = 15, height = 15)

DimPlot(mye.patients.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Patient", shuffle = T, cells.highlight = list(P3 = WhichCells(mye.patients.seurat, expression = Patient == "P3")) ) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("myeloid_results/UMAP Patient 3.pdf", width = 15, height = 15)

DimPlot(mye.patients.seurat, label = F, pt.size = 2, label.size = 10, group.by = "Tissue", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("myeloid_results/UMAP Tissue.pdf", width = 15, height = 15)

DimPlot(mye.patients.seurat, label = F, pt.size = 2, label.size = 10, group.by = "predicted.celltype.l1", reduction = "umap", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("myeloid_results/UMAP Predicted Celltype L1.pdf", width = 15, height = 15)

DimPlot(mye.patients.seurat, label = F, pt.size = 2, label.size = 10, group.by = "predicted.celltype.l2", reduction = "umap", shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("myeloid_results/UMAP Predicted Celltype L2.pdf", width = 15, height = 15)

## Save the myeloid seurat object
saveRDS(mye.patients.seurat, "Seurat_Objects/myeloid.cells.seurat.RDS")
#mye.patients.seurat <- readRDS(file = "Seurat_Objects/myeloid.cells.seurat.RDS")


customVln(mye.patients.seurat, features = c("SDC1", "SDC2","SDC3", "SDC4") , name = "myeloid_results/syndecan.pdf", width = 40, height = 20, ncol = 4)

