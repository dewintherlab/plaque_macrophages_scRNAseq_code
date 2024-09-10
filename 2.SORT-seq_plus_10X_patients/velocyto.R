#=================================================================================================================================
## Use velocyto to geretate cell trajectories based on spliced / unspliced read ratios
## velocyto loom object prepped from .bam files by running velocyto.py
##================================================================================================================================
dir.create("velocyto_results", showWarnings = F)

## Prep objects
# Read velocyto loom file
all.samples.loom <- ReadVelocity(file = "raw_data/merged.loom")

# Convert to Seurat object
all.samples.velocyto.seurat <- as.Seurat(x = all.samples.loom)

##================================================================================================================================
## Keep only the cells we still have in our filtered, final main object
# Check Cell names
colnames(refined.patients.seurat)[1:5]
tail(colnames(refined.patients.seurat))
colnames(all.samples.velocyto.seurat)[1:5]
tail(all.samples.velocyto.seurat)

# Convert velocyto object cell names to the format in our 'normal' objects
velo.cell.names <- colnames(all.samples.velocyto.seurat)
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g002_ranger_out:([A-Z]+)x", replacement = "PBMC_\\1-1_1")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g003_ranger_out:([A-Z]+)x", replacement = "plaque_\\1-1_1")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g004_ranger_out:([A-Z]+)x", replacement = "PBMC_\\1-1_2")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g005_ranger_out:([A-Z]+)x", replacement = "plaque_\\1-1_2")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g006_ranger_out:([A-Z]+)x", replacement = "PBMC_\\1-1_3")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g007_ranger_out:([A-Z]+)x", replacement = "plaque_\\1-1_3")

# Rename the velo cells
all.samples.velocyto.seurat <- RenameCells(all.samples.velocyto.seurat, new.names = velo.cell.names)

# And subset the object to contain only the filtered cells
all.samples.velocyto.seurat <- subset(all.samples.velocyto.seurat, cells = colnames(refined.patients.seurat))

# Also subset vice versa in case we are missing some cells in the velo object
refined.samples.velocyto.seurat <- subset(refined.patients.seurat, cells = colnames(all.samples.velocyto.seurat))

##================================================================================================================================
## Match on features as well
# Check feature names
DefaultAssay(refined.samples.velocyto.seurat) <- "RNA"
sort(row.names(refined.samples.velocyto.seurat))[1:5]
sort(row.names(all.samples.velocyto.seurat))[1:5]

# And subset the object
all.samples.velocyto.seurat <- subset(all.samples.velocyto.seurat, features = row.names(refined.samples.velocyto.seurat))

##================================================================================================================================
## Transfer velocyto data to our main object
# Assays
refined.samples.velocyto.seurat[["spliced"]]   <- CreateAssayObject(data = GetAssayData(all.samples.velocyto.seurat, assay = "spliced"))
refined.samples.velocyto.seurat[["unspliced"]] <- CreateAssayObject(data = GetAssayData(all.samples.velocyto.seurat, assay = "unspliced"))
refined.samples.velocyto.seurat[["ambiguous"]] <- CreateAssayObject(data = GetAssayData(all.samples.velocyto.seurat, assay = "ambiguous"))

# Metadata
refined.samples.velocyto.seurat <- AddMetaData(refined.samples.velocyto.seurat, metadata = all.samples.velocyto.seurat@meta.data)
head(refined.samples.velocyto.seurat@meta.data)
tail(refined.samples.velocyto.seurat@meta.data)

# Save the object
saveRDS(refined.samples.velocyto.seurat, "Seurat_Objects/main.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")
#refined.samples.velocyto.seurat <- readRDS(file = "Seurat_Objects/main.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")

##================================================================================================================================
## Transfer velocyto data to our myeloid object
# Subset the velo object to contain only the myeloid cells
mye.refined.samples.velocyto.seurat <- mye.patients.seurat
mye.samples.velocyto.seurat <- subset(all.samples.velocyto.seurat, cells = colnames(mye.patients.seurat))

# Also subset vice versa because we might be missing some cells in the velo object
mye.refined.samples.velocyto.seurat <- subset(mye.refined.samples.velocyto.seurat, cells = colnames(mye.samples.velocyto.seurat))

##================================================================================================================================
## Match on features as well
DefaultAssay(mye.refined.samples.velocyto.seurat) <- "RNA"

# Subset the object
mye.samples.velocyto.seurat <- subset(mye.samples.velocyto.seurat, features = row.names(mye.refined.samples.velocyto.seurat))

##================================================================================================================================
## Transfer velocyto data to our main object
# Assays
mye.refined.samples.velocyto.seurat[["spliced"]]   <- CreateAssayObject(data = GetAssayData(mye.samples.velocyto.seurat, assay = "spliced"))
mye.refined.samples.velocyto.seurat[["unspliced"]] <- CreateAssayObject(data = GetAssayData(mye.samples.velocyto.seurat, assay = "unspliced"))
mye.refined.samples.velocyto.seurat[["ambiguous"]] <- CreateAssayObject(data = GetAssayData(mye.samples.velocyto.seurat, assay = "ambiguous"))

# Metadata
mye.refined.samples.velocyto.seurat <- AddMetaData(mye.refined.samples.velocyto.seurat, metadata = mye.samples.velocyto.seurat@meta.data)
head(mye.refined.samples.velocyto.seurat@meta.data)
tail(mye.refined.samples.velocyto.seurat@meta.data)

# Save the object
saveRDS(mye.refined.samples.velocyto.seurat, "Seurat_Objects/myeloid.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")
#mye.refined.samples.velocyto.seurat <- readRDS(file = "Seurat_Objects/myeloid.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")

##================================================================================================================================
## Transfer velocyto data to our myeloid object
# Subset the velo object to contain only the myeloid cells
mye.refined.samples.velocyto.seurat <- mye.patients.seurat
mye.samples.velocyto.seurat <- subset(all.samples.velocyto.seurat, cells = colnames(mye.patients.seurat))

# Also subset vice versa because we might be missing some cells in the velo object
mye.refined.samples.velocyto.seurat <- subset(mye.refined.samples.velocyto.seurat, cells = colnames(mye.samples.velocyto.seurat))

##================================================================================================================================
## Match on features as well
DefaultAssay(mye.refined.samples.velocyto.seurat) <- "RNA"

# Subset the object
mye.samples.velocyto.seurat <- subset(mye.samples.velocyto.seurat, features = row.names(mye.refined.samples.velocyto.seurat))

##================================================================================================================================
## Transfer velocyto data to our main object
# Assays
mye.refined.samples.velocyto.seurat[["spliced"]]   <- CreateAssayObject(data = GetAssayData(mye.samples.velocyto.seurat, assay = "spliced"))
mye.refined.samples.velocyto.seurat[["unspliced"]] <- CreateAssayObject(data = GetAssayData(mye.samples.velocyto.seurat, assay = "unspliced"))
mye.refined.samples.velocyto.seurat[["ambiguous"]] <- CreateAssayObject(data = GetAssayData(mye.samples.velocyto.seurat, assay = "ambiguous"))

# Metadata
mye.refined.samples.velocyto.seurat <- AddMetaData(mye.refined.samples.velocyto.seurat, metadata = mye.samples.velocyto.seurat@meta.data)
head(mye.refined.samples.velocyto.seurat@meta.data)
tail(mye.refined.samples.velocyto.seurat@meta.data)

# Save the object
saveRDS(mye.refined.samples.velocyto.seurat, "Seurat_Objects/myeloid.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")
#mye.refined.samples.velocyto.seurat <- readRDS(file = "Seurat_Objects/myeloid.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")


##================================================================================================================================
## Transfer velocyto data to our macrophage object
# Subset the velo object to contain only the myeloid cells
mac.refined.samples.velocyto.seurat <- plaque.macrophages.seurat
mac.samples.velocyto.seurat <- subset(all.samples.velocyto.seurat, cells = colnames(plaque.macrophages.seurat))

# Also subset vice versa because we might be missing some cells in the velo object
mac.refined.samples.velocyto.seurat <- subset(mac.refined.samples.velocyto.seurat, cells = colnames(mac.samples.velocyto.seurat))

##================================================================================================================================
## Match on features as well
DefaultAssay(mac.refined.samples.velocyto.seurat) <- "RNA"

# Subset the object
mac.samples.velocyto.seurat <- subset(mac.samples.velocyto.seurat, features = row.names(mac.refined.samples.velocyto.seurat))

##================================================================================================================================
## Transfer velocyto data to our main object
# Assays
mac.refined.samples.velocyto.seurat[["spliced"]]   <- CreateAssayObject(data = GetAssayData(mac.samples.velocyto.seurat, assay = "spliced"))
mac.refined.samples.velocyto.seurat[["unspliced"]] <- CreateAssayObject(data = GetAssayData(mac.samples.velocyto.seurat, assay = "unspliced"))
mac.refined.samples.velocyto.seurat[["ambiguous"]] <- CreateAssayObject(data = GetAssayData(mac.samples.velocyto.seurat, assay = "ambiguous"))

# Metadata
mac.refined.samples.velocyto.seurat <- AddMetaData(mac.refined.samples.velocyto.seurat, metadata = mac.samples.velocyto.seurat@meta.data)
head(mac.refined.samples.velocyto.seurat@meta.data)
tail(mac.refined.samples.velocyto.seurat@meta.data)

# Save the object
saveRDS(mac.refined.samples.velocyto.seurat, "Seurat_Objects/macrophage.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")
#mac.refined.samples.velocyto.seurat <- readRDS(file = "Seurat_Objects/macropahge.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")


##================================================================================================================================
## Transfer velocyto data to our monocyte object
# Subset the velo object to contain only the myeloid cells
mono.refined.samples.velocyto.seurat <- pbmc.monocytes.seurat
mono.samples.velocyto.seurat <- subset(all.samples.velocyto.seurat, cells = colnames(pbmc.monocytes.seurat))

# Also subset vice versa because we might be missing some cells in the velo object
mono.refined.samples.velocyto.seurat <- subset(mono.refined.samples.velocyto.seurat, cells = colnames(mono.samples.velocyto.seurat))

##================================================================================================================================
## Match on features as well
DefaultAssay(mono.refined.samples.velocyto.seurat) <- "RNA"

# Subset the object
mono.samples.velocyto.seurat <- subset(mono.samples.velocyto.seurat, features = row.names(mono.refined.samples.velocyto.seurat))

##================================================================================================================================
## Transfer velocyto data to our main object
# Assays
mono.refined.samples.velocyto.seurat[["spliced"]]   <- CreateAssayObject(data = GetAssayData(mono.samples.velocyto.seurat, assay = "spliced"))
mono.refined.samples.velocyto.seurat[["unspliced"]] <- CreateAssayObject(data = GetAssayData(mono.samples.velocyto.seurat, assay = "unspliced"))
mono.refined.samples.velocyto.seurat[["ambiguous"]] <- CreateAssayObject(data = GetAssayData(mono.samples.velocyto.seurat, assay = "ambiguous"))

# Metadata
mono.refined.samples.velocyto.seurat <- AddMetaData(mono.refined.samples.velocyto.seurat, metadata = mono.samples.velocyto.seurat@meta.data)
head(mono.refined.samples.velocyto.seurat@meta.data)
tail(mono.refined.samples.velocyto.seurat@meta.data)

# Save the object
saveRDS(mono.refined.samples.velocyto.seurat, "Seurat_Objects/monocytes.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")
#mono.refined.samples.velocyto.seurat <- readRDS(file = "Seurat_Objects/monocytes.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")



##================================================================================================================================
## Clean up to save some memory
rm(all.samples.loom)
rm(all.samples.velocyto.seurat)
rm(mye.samples.velocyto.seurat)
rm(mac.samples.velocyto.seurat)
rm(mono.samples.velocyto.seurat)
rm(velo.cell.names)
rm(mye.cell.names)
rm(cell.names)


##================================================================================================================================
## Run velocyto natively on myeloid populations
DefaultAssay(mye.refined.samples.velocyto.seurat) <- "spliced"

# Run velocyto scoring
mye.refined.samples.velocyto.seurat <- RunVelocity(object = mye.refined.samples.velocyto.seurat, deltaT = 1, kCells = 25, fit.quantile = 0.02)

# Visualise
ident.colors       <- M.int_refined.pop.colors

# Get new refined idents
cell.colors        <- as.vector(Idents(integrated.mye.seurat))
names(cell.colors) <- names(Idents(integrated.mye.seurat))
cell.colors        <- cell.colors[names(cell.colors) %in% names(Idents(mye.refined.samples.velocyto.seurat))]

mono.cells         <- as.vector(Idents(mye.refined.samples.velocyto.seurat))[grep("SLAN", Idents(mye.refined.samples.velocyto.seurat))]
names(mono.cells)  <- names(Idents(mye.refined.samples.velocyto.seurat))[grep("SLAN", Idents(mye.refined.samples.velocyto.seurat))]
cell.colors        <- c(cell.colors, mono.cells)

Idents(mye.refined.samples.velocyto.seurat) <- cell.colors

# Set colors
for(theIdent in unique(Idents(mye.refined.samples.velocyto.seurat))){
  cell.colors[cell.colors == theIdent] <- ident.colors[theIdent]
}

pdf("velocyto_results/mye.cells.velo_embeddings_n200-40.pdf",width = 15, height = 15)
show.velocity.on.embedding.cor(emb                = Embeddings(object    = mye.refined.samples.velocyto.seurat, 
                                                               reduction = "umap"), 
                               vel                = Tool(object = mye.refined.samples.velocyto.seurat, 
                                                         slot = "RunVelocity"),
                               n                  = 200, 
                               scale              = "sqrt", 
                               cell.colors        = ac(x     = cell.colors, alpha = 0.5), 
                               cex                = 2,
                               arrow.scale        = 2, 
                               show.grid.flow     = TRUE, 
                               min.grid.cell.mass = 0.5,
                               grid.n             = 40, 
                               arrow.lwd          = 2, 
                               do.par             = FALSE, 
                               cell.border.alpha  = 0.1,
                               n.cores            = 6,
                               nPcs               = 30)
dev.off()

##================================================================================================================================
## Run velocyto natively on macrophage populations
DefaultAssay(mac.refined.samples.velocyto.seurat) <- "spliced"

# Run velocyto scoring
mac.refined.samples.velocyto.seurat <- RunVelocity(object = mac.refined.samples.velocyto.seurat, deltaT = 1, kCells = 25, fit.quantile = 0.02)

# Visualise
ident.colors       <- M.int_refined.pop.colors
# Get new refined idents
cell.colors        <- as.vector(Idents(integrated.mye.seurat))
names(cell.colors) <- names(Idents(integrated.mye.seurat))
cell.colors        <- cell.colors[names(cell.colors) %in% names(Idents(mac.refined.samples.velocyto.seurat))]

Idents(mac.refined.samples.velocyto.seurat) <- cell.colors

for(theIdent in unique(Idents(mac.refined.samples.velocyto.seurat))){
  cell.colors[cell.colors == theIdent] <- ident.colors[theIdent]
}

pdf("velocyto_results/mac.cells.velo_embeddings_n200-40.pdf",width = 15, height = 15)
show.velocity.on.embedding.cor(emb                = Embeddings(object    = mac.refined.samples.velocyto.seurat, 
                                                               reduction = "umap"), 
                               vel                = Tool(object = mac.refined.samples.velocyto.seurat, 
                                                         slot = "RunVelocity"),
                               n                  = 200, 
                               scale              = "sqrt", 
                               cell.colors        = ac(x     = cell.colors, 
                                                       alpha = 0.5), 
                               cex                = 2,
                               arrow.scale        = 2, 
                               show.grid.flow     = TRUE, 
                               min.grid.cell.mass = 0.5,
                               grid.n             = 40, 
                               arrow.lwd          = 2, 
                               do.par             = FALSE, 
                               cell.border.alpha  = 0.1,
                               n.cores            = 6,
                               nPcs               = 30)
dev.off()

##================================================================================================================================
## Run velocyto natively on monocyte populations
DefaultAssay(mono.refined.samples.velocyto.seurat) <- "spliced"

# Run velocyto scoring
mono.refined.samples.velocyto.seurat <- RunVelocity(object = mono.refined.samples.velocyto.seurat, deltaT = 1, kCells = 25, fit.quantile = 0.02)

# Visualise
ident.colors       <- M_refined.pop.colors
cell.colors        <- as.vector(Idents(mono.refined.samples.velocyto.seurat))
names(cell.colors) <- names(Idents(mono.refined.samples.velocyto.seurat))

for(theIdent in unique(Idents(mono.refined.samples.velocyto.seurat))){
  cell.colors[cell.colors == theIdent] <- ident.colors[theIdent]
}

pdf("velocyto_results/mono.cells.velo_embeddings_n200-40.pdf",width = 15, height = 15)
show.velocity.on.embedding.cor(emb                = Embeddings(object    = mono.refined.samples.velocyto.seurat, 
                                                               reduction = "umap"), 
                               vel                = Tool(object = mono.refined.samples.velocyto.seurat, 
                                                         slot = "RunVelocity"),
                               n                  = 200, 
                               scale              = "sqrt", 
                               cell.colors        = ac(x     = cell.colors, 
                                                       alpha = 0.5), 
                               cex                = 2,
                               arrow.scale        = 2, 
                               show.grid.flow     = TRUE, 
                               min.grid.cell.mass = 0.5,
                               grid.n             = 40, 
                               arrow.lwd          = 2, 
                               do.par             = FALSE, 
                               cell.border.alpha  = 0.1,
                               n.cores            = 6,
                               nPcs               = 30)
dev.off()



##================================================================================================================================
## Try re-running the embeddings
# Fix the ocunts slo for 'spliced' assay
tmp.assay <- GetAssayData(mye.refined.samples.velocyto.seurat, assay = "spliced")
mye.refined.samples.velocyto.seurat <- SetAssayData(mye.refined.samples.velocyto.seurat, assay = "spliced", slot = "counts", new.data = tmp.assay)
rm(tmp.assay)

# Run the normalisation and embedding
new_emb.mye.velo.seurat <- SCTransform(object = mye.refined.samples.velocyto.seurat, assay = "spliced")
new_emb.mye.velo.seurat <- RunPCA(object = new_emb.mye.velo.seurat, verbose = FALSE)
new_emb.mye.velo.seurat <- FindNeighbors(object = new_emb.mye.velo.seurat, dims = 1:20)
new_emb.mye.velo.seurat <- FindClusters(object = new_emb.mye.velo.seurat)
new_emb.mye.velo.seurat <- RunUMAP(object = new_emb.mye.velo.seurat, dims = 1:20)
new_emb.mye.velo.seurat <- RunVelocity(object = new_emb.mye.velo.seurat, deltaT = 1, kCells = 25, fit.quantile = 0.02)

# And have a look
ident.colors            <- (scales::hue_pal())(n = length(x = levels(x = new_emb.mye.velo.seurat)))
names(x = ident.colors) <- levels(x = new_emb.mye.velo.seurat)
cell.colors             <- ident.colors[Idents(object = new_emb.mye.velo.seurat)]
names(x = cell.colors)  <- colnames(x = new_emb.mye.velo.seurat)

show.velocity.on.embedding.cor(emb = Embeddings(object = new_emb.mye.velo.seurat, reduction = "umap"), 
                               vel = Tool(object = new_emb.mye.velo.seurat, slot = "RunVelocity"),
                               n = 200, 
                               scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, 
                               show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, 
                               grid.n = 40, 
                               arrow.lwd = 1, 
                               do.par = FALSE, 
                               cell.border.alpha = 0.1)

# Remove abject again to save space
rm(new_emb.mye.velo.seurat)