#======================================================================
## Save stuff as python objects
#======================================================================
# Setup virtual env
# For pyVIA
use_condaenv('ViaEnv')

# For cellrank adn decoupler
use_condaenv('cellrank', conda = "/Users/koenprange/mambaforge/bin/conda")
rm(py)
reticulate::py_config()

## Prep and convert the mac + mono seurat object (with velocyto data embedded)
# Add idents as 'final.pop.idents' metadata
final.pop.call.from_full.integrated.mye.velocyto.seurat <- AddMetaData(final.pop.call.from_full.integrated.mye.velocyto.seurat, metadata = Idents(final.pop.call.from_full.integrated.mye.velocyto.seurat), col.name = "final.pop.idents")

# Add archetypes
archetypes <- as.vector(Idents(final.pop.call.from_full.integrated.mye.velocyto.seurat))
archetypes[grep("Foamy", archetypes)]        <- "iLAM"
archetypes[grep("Lipid", archetypes)]        <- "LAM"
archetypes[grep("Resident", archetypes)]     <- "Resident-like"
archetypes[grep("Monocyte-derived", archetypes)] <- "Inflammatory"
archetypes[grep("Monocytes", archetypes)]    <- "Monocytes"

final.pop.call.from_full.integrated.mye.velocyto.seurat <- AddMetaData(final.pop.call.from_full.integrated.mye.velocyto.seurat, metadata = archetypes, col.name = "archetypes")


# Convert all metadata to character so we can read it properly in python
i <- sapply(final.pop.call.from_full.integrated.mye.velocyto.seurat@meta.data, is.factor)
final.pop.call.from_full.integrated.mye.velocyto.seurat@meta.data[i] <- lapply(final.pop.call.from_full.integrated.mye.velocyto.seurat@meta.data[i], as.character)

# Check the UMAP
DimPlot(final.pop.call.from_full.integrated.mye.velocyto.seurat, cols = archetype.colors, group.by = "archetypes") + NoLegend()

## Convert to h5seurat and then h5ad
# Set default assay
DefaultAssay(final.pop.call.from_full.integrated.mye.velocyto.seurat) <- "RNA"

# Manually populate median_umi slot of the SCTmodels. This slot form Seurat 4.1 is not set in our older object and UpdateSeurat() does not yet fix it...
slot(final.pop.call.from_full.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[1]], 'median_umi') = median(final.pop.call.from_full.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[1]]@cell.attributes$umi)
slot(final.pop.call.from_full.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[2]], 'median_umi') = median(final.pop.call.from_full.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[2]]@cell.attributes$umi)
slot(final.pop.call.from_full.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[3]], 'median_umi') = median(final.pop.call.from_full.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[3]]@cell.attributes$umi)
slot(final.pop.call.from_full.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[4]], 'median_umi') = median(final.pop.call.from_full.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[4]]@cell.attributes$umi)

head(row.names(final.pop.call.from_full.integrated.mye.velocyto.seurat))
"NOC2L" %in% head(row.names(final.pop.call.from_full.integrated.mye.velocyto.seurat))
"HSPA6" %in% row.names(final.pop.call.from_full.integrated.mye.velocyto.seurat)

# Convert the object first to h5Seurat
SaveH5Seurat(final.pop.call.from_full.integrated.mye.velocyto.seurat, filename = "Seurat_Objects/myeloid.velocyto.refined.h5Seurat", overwrite = T)

# And then to h5ad
Convert(source = "Seurat_Objects/myeloid.velocyto.refined.h5Seurat", dest= "h5ad", overwrite = T)

# Bridge the population color objects to the python session
pop_cols     <- as.list(M.int_refined.pop.colors)
py$pop_cols  <- r_to_py(pop_cols)
arch_cols    <- as.list(archetype.colors)
py$arch_cols <- r_to_py(arch_cols)

#======================================================================
## Prep and convert the mac seurat object (with velocyto data embedded)
# Add idents as 'final.pop.idents' metadata
final.pop.call.integrated.mye.velocyto.seurat <- AddMetaData(final.pop.call.integrated.mye.velocyto.seurat, metadata = Idents(final.pop.call.integrated.mye.velocyto.seurat), col.name = "final.pop.idents")

# Add archetypes
archetypes <- as.vector(Idents(final.pop.call.integrated.mye.velocyto.seurat))
archetypes[grep("Foamy", archetypes)]        <- "iLAM"
archetypes[grep("Lipid", archetypes)]        <- "LAM"
archetypes[grep("Resident", archetypes)]     <- "Resident-like"
archetypes[grep("Monocyte-derived", archetypes)] <- "Inflammatory"

final.pop.call.integrated.mye.velocyto.seurat <- AddMetaData(final.pop.call.integrated.mye.velocyto.seurat, metadata = archetypes, col.name = "archetypes")


# Convert all metadata to character so we can read it properly in python
i <- sapply(final.pop.call.integrated.mye.velocyto.seurat@meta.data, is.factor)
final.pop.call.integrated.mye.velocyto.seurat@meta.data[i] <- lapply(final.pop.call.integrated.mye.velocyto.seurat@meta.data[i], as.character)

# Check the UMAP
DimPlot(final.pop.call.integrated.mye.velocyto.seurat, cols = archetype.colors, group.by = "archetypes", pt.size = 3) + NoLegend()

## Convert to h5seurat and then h5ad
# Set default assay
DefaultAssay(final.pop.call.integrated.mye.velocyto.seurat) <- "RNA"

# Manually populate median_umi slot of the SCTmodels. This slot form Seurat 4.1 is not set in our older object and UpdateSeurat() does not yet fix it...
slot(final.pop.call.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[1]], 'median_umi') = median(final.pop.call.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[1]]@cell.attributes$umi)
slot(final.pop.call.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[2]], 'median_umi') = median(final.pop.call.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[2]]@cell.attributes$umi)
slot(final.pop.call.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[3]], 'median_umi') = median(final.pop.call.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[3]]@cell.attributes$umi)
slot(final.pop.call.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[4]], 'median_umi') = median(final.pop.call.integrated.mye.velocyto.seurat$SCT@SCTModel.list[[4]]@cell.attributes$umi)

tail(row.names(final.pop.call.integrated.mye.velocyto.seurat))
tail(row.names(final.pop.call.from_full.integrated.mye.velocyto.seurat))
tail(colnames(final.pop.call.integrated.mye.velocyto.seurat))
tail(colnames(final.pop.call.from_full.integrated.mye.velocyto.seurat))

"NOC2L" %in% head(row.names(final.pop.call.integrated.mye.velocyto.seurat))
"HSPA6" %in% row.names(final.pop.call.integrated.mye.velocyto.seurat)

# Keep only relevant metadata
final.pop.call.integrated.mye.velocyto.seurat@meta.data <- final.pop.call.integrated.mye.velocyto.seurat@meta.data[,colnames(final.pop.call.integrated.mye.velocyto.seurat@meta.data) %in% colnames(final.pop.call.from_full.integrated.mye.velocyto.seurat@meta.data)]

# Convert the object first to h5Seurat
SaveH5Seurat(final.pop.call.integrated.mye.velocyto.seurat, filename = "Seurat_Objects/mac.velocyto.refined.h5Seurat", overwrite = T, verbose = T)

# And then to h5ad
Convert(source = "Seurat_Objects/mac.velocyto.refined.h5Seurat", dest= "h5ad", overwrite = T)


# Create result dirs
dir.create("Cellrank", showWarnings = F)
dir.create("decoupler", showWarnings = F)
dir.create("Cellrank/macs_only", showWarnings = F)
dir.create("decoupler/macs_only", showWarnings = F)
dir.create("CellRank/Archetypes", showWarnings = F)
dir.create("CellRank/Archetypes_macs_only", showWarnings = F)

