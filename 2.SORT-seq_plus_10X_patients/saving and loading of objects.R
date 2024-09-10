 #=================================================================================================================================
## Overview of endpoint or branchpoint objects, saved to disk.
## Can be used to restart analysis at a midway or end point without having to re-run the whole thing.
##================================================================================================================================
dir.create("Seurat_Objects", showWarnings = F)

#=================================================================================================================================
## Objects to save (commented out so we don't accidentally overwrite our saved stuff)
##================================================================================================================================
## Gene ontology data
# saveRDS(human, "Seurat_Objects/human.mart.data.RDS")
# saveRDS(mouse, "Seurat_Objects/mouse.mart.data.RDS")
# 
# ## Seurat objects
# # List of 6 raw seurat objects, freshly loaded
# saveRDS(samples.seurat, "Seurat_Objects/main.seurat.raw.RDS")
# 
# # List of 3 seurat objects (per patient), SCtransformed, embedded, % mt etc. calculated
# saveRDS(patients.seurat, "Seurat_Objects/main.patients.seurat.prefilter.RDS")
# 
# # List of 3 seurat objects (per patients), SCtransformed, embedded, % mt etc. filtered
# saveRDS(patients.seurat, "Seurat_Objects/main.patients.seurat.filtered.RDS")
# 
# # Integrated seurat object, embedded
# saveRDS(patients.seurat, "Seurat_Objects/main.patients.seurat.filtered.integrated.RDS")
# 
# # Integrated seurat object, embedded, mapped to reference
# saveRDS(patients.seurat, "Seurat_Objects/main.patients.seurat.filtered.integrated.mapped.RDS")
#
# # Integrated seurat object, embedded, mapped to reference, populations refined
# saveRDS(refined.patients.seurat, "Seurat_Objects/main.patient.integrated.mapped.refined_populations.seurat.RDS")
#
# # Myeloid cells seurat object
# saveRDS(mye.patients.seurat, "Seurat_Objects/myeloid.cells.seurat.RDS")
# 
# # PBMC only monocytes eurat object
# saveRDS(pbmc.monocytes.seurat, "Seurat_Objects/pbmc.monocytes.seurat.RDS")
#
# # Plaque only macrophages seurat object
# saveRDS(plaque.macrophages.seurat, "Seurat_Objects/plaque.macrophages.seurat.RDS")
#
# # Integrated 43 patietns cel-seq and 10X myeloid cells seurat object
# saveRDS(integrated.mye.seurat, "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")
#
# # Full set 43 patients cel-seq and 10X cells seurat object
# saveRDS(integrated.full.seurat, "Seurat_Objects/full.43p_10X.integrated.seurat.RDS")
#
# # Full set 43 patients cel-seq and 10X cells cleaned and idents resolved seurat object
# saveRDS(integrated.full.seurat, "Seurat_Objects/full.43p_10X.integrated.cleaned.seurat.RDS")
#
# # From Full set 43 patients cel-seq and 10X cells cleaned and idents resolved myeloid  seurat object
# saveRDS(from_full.integrated.mye.seurat, "Seurat_Objects/from_full.43p_10X.integrated.cleaned.mye.seurat.RDS")
#
# # ssGSEA data ran on mye.integrated.seurat
# saveRDS(ssGSEA.integrated.mye.seurat, "Seurat_Objects/ssGSEA_data.43p_10X.integrated.RDS")
#
# # Myeloid cells refined Seurat object
# saveRDS(mye.patients.seurat, "Seurat_Objects/myeloid.cells.refined.seurat.RDS")
#
# # 43p seurat object, idents refined
# saveRDS(full.43p.seurat, file = "Seurat_Objects/full_43p.idents_refiened.RDS")
#

# ## Final population update
# saveRDS(final.pop.call.integrated.mye.seurat, file = "Seurat_Objects/final.pop.call.integrated.mye.seurat.Rds")
# saveRDS(from_full.integrated.mye.seurat, file = "Seurat_Objects/from_full.integrated.mye.seurat.Rds")
# saveRDS(final.pop.call.from_full.integrated.mye.seurat, file = "Seurat_Objects/final.pop.call.from_full.integrated.mye.seurat.Rds")
# saveRDS(final.pop.call.integrated.full.seurat, file = "Seurat_Objects/final.pop.call.integrated.full.seurat.Rds")
# saveRDS(final.pop.call.integrated.mye.velocyto.seurat, "Seurat_Objects/final.pop.call.mye.velo_plaque.integrated.mapped.refined_populations.seurat.RDS")
# saveRDS(final.pop.call.from_full.integrated.mye.velocyto.seurat, "Seurat_Objects/final.pop.call.mye.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")
#
# ## Marker genes
# # All markers from all patients and idents
# saveRDS(all.patients.seurat.markers, "Seurat_Objects/main.patients.markers.RDS")
# 
# # PBMC monocyte markers
# saveRDS(pbmc.monocytes.seurat.markers, "Seurat_Objects/PBMC.monocytes.markers.RDS")
#
# # Plaque macrophage markers
# saveRDS(plaque.macrophages.seurat.markers, "Seurat_Objects/plaque.macrophages.markers.RDS")
# 
# # Integrated myeloid markers
# saveRDS(integrated.mye.seurat.markers, file = "Seurat_Objects/integrated.mye.seurat.markers")
#
#
# ## Gene Ontologies
# 
# 
# ## Monocle objects
# 
# 
# ## Control samples
# saveRDS(ctrl.combined,         "Seurat_Objects/control.combined.Wang.seurat.RDS")
# saveRDS(ctrl_patient.combined, "Seurat_Objects/control_patient.combined.Wang.seurat.RDS")
# 
# ## Color Objects
# saveRDS(full_set.colors,          file = "Seurat_Objects/full_set.colors.RDS")
# saveRDS(M.int_refined.pop.colors, file = "Seurat_Objects/M.int_refined.pop.colors.RDS")
#
# ## Random odds and ends
# saveRDS(vars.to.regress,       file = "Seurat_Objects/vars.to.regress.Rds")
# saveRDS(Mye.markers.dgi,       file = "Seurat_Objects/mye_markers.DGI_object.RDS")
# saveRDS(exp.m.s.cpdb.relevant, file = "Seurat_Objects/cpdb_results.clean.RDS")
#
#
#=================================================================================================================================
## Objects to load
##================================================================================================================================
## Gene ontology data
ont_sets <- readList("cp_and_h.gmt.txt")
human    <- readRDS(file = "Seurat_Objects/human.mart.data.RDS")
mouse    <- readRDS(file = "Seurat_Objects/mouse.mart.data.RDS")

## Seurat objects
# List of 6 raw seurat objects, freshly loaded
samples.seurat <- readRDS(file = "Seurat_Objects/main.seurat.raw.RDS")

# List of 3 seurat objects (per patient), SCtransformed, embedded, % mt etc. calculated
patients.seurat <- readRDS(file = "Seurat_Objects/main.patients.seurat.prefilter.RDS")

# List of 3 seurat objects (per patient), SCtransformed, embedded, % mt etc. filtered
patients.seurat <- readRDS(file = "Seurat_Objects/main.patients.seurat.filtered.RDS")

# Integrated seurat object, embedded
patients.seurat <- readRDS(file = "Seurat_Objects/main.patients.seurat.integrated.RDS")

# Integrated seurat object, embedded, mapped to reference
patients.seurat <- readRDS(file = "Seurat_Objects/main.patients.seurat.filtered.mapped.RDS")

# Integrated seurat object, embedded, mapped to reference, populations refined
refined.patients.seurat <- readRDS(file = "Seurat_Objects/main.patient.integrated.mapped.refined_populations.seurat.RDS")

# Myeloid cells seurat object
mye.patients.seurat <- readRDS(file = "Seurat_Objects/myeloid.cells.seurat.RDS")

# PBMC only monocytes seurat object
pbmc.monocytes.seurat <- readRDS(file = "Seurat_Objects/pbmc.monocytes.seurat.RDS")

# Plaque only macrophages seurat object
plaque.macrophages.seurat <- readRDS(file = "Seurat_Objects/plaque.macrophages.seurat.RDS")

# Myeloid cells refined Seurat object
mye.patients.seurat <- readRDS(file = "Seurat_Objects/myeloid.cells.refined.seurat.RDS")

# Integrated 43 patietns cel-seq and 10X myeloid cells seurat object
integrated.mye.seurat <- readRDS(file = "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")

# Full set 43 patients cel-seq and 10X cells seurat object
integrated.full.seurat <- readRDS(file = "Seurat_Objects/full.43p_10X.integrated.seurat.RDS")

# Full set 43 patients cel-seq and 10X cells cleaned and idents resolved seurat object
integrated.full.seurat <- readRDS(file = "Seurat_Objects/full.43p_10X.integrated.cleaned.seurat.RDS")

# Full set 43 patients cel-seq and 10X cells cleaned and idents resolved to archetypes seurat object
integrated.full.seurat <- readRDS(file =  "Seurat_Objects/full.43p_10X.integrated.cleaned.archetypes.seurat.RDS")

# From Full set 43 patients cel-seq and 10X cells cleaned and idents resolved myeloid  seurat object
from_full.integrated.mye.seurat <- readRDS(file = "Seurat_Objects/from_full.43p_10X.integrated.cleaned.mye.seurat.RDS")

# ssGSEA data ran on mye.integrated.seurat
ssGSEA.integrated.mye.seurat <- readRDS(file = "Seurat_Objects/ssGSEA_data.43p_10X.integrated.RDS")

# 43p seurat obejct with refiend idents
full.43p.seurat <- readRDS(file = "Seurat_Objects/full_43p.idents_refiened.RDS")

# Final population update
final.pop.call.from_full.integrated.mye.seurat          <- readRDS(file = "Seurat_Objects/final.pop.call.from_full.integrated.mye.seurat.Rds")
final.pop.call.from_full.integrated.mye.velocyto.seurat <- readRDS(file = "Seurat_Objects/final.pop.call.mye.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")
final.pop.call.integrated.mye.velocyto.seurat           <- readRDS(file = "Seurat_Objects/final.pop.call.mye.velo_plaque.integrated.mapped.refined_populations.seurat.RDS")
final.pop.call.integrated.full.seurat                   <- readRDS(file = "Seurat_Objects/final.pop.call.integrated.full.seurat.Rds")


## Marker genes
# All markers from all patients and idents
all.patients.seurat.markers       <- readRDS(file = "Seurat_Objects/main.patients.markers.RDS")
pbmc.monocytes.seurat.markers     <- readRDS(file = "Seurat_Objects/PBMC.monocytes.markers.RDS")
plaque.macrophages.seurat.markers <- readRDS(file = "Seurat_Objects/plaque.macrophages.markers.RDS")
integrated.mye.seurat.markers     <- readRDS(file = "Seurat_Objects/integrated.mye.seurat.markers")

## Gene Ontologies


## Monocle objects


## Control samples
ctrl.combined         <- readRDS(file = "Seurat_Objects/control.combined.Wang.seurat.RDS")
ctrl_patient.combined <- readRDS(file = "Seurat_Objects/control_patient.combined.Wang.seurat.RDS")

## Random odds and ends
vars.to.regress       <- readRDS(file = "Seurat_Objects/vars.to.regress.Rds")
Mye.markers.dgi       <- readRDS(file = "Seurat_Objects/mye_markers.DGI_object.RDS")
Mye.type.markers.dgi  <- readRDS(file = "Seurat_Objects/mye_type.markers.DGI_object.RDS")
Mye.markers.dgi       <- readRDS(file = "Seurat_Objects/mye_markers_plus_LRI.DGI_object.RDS")
Mye.type.markers.dgi  <- readRDS(file = "Seurat_Objects/mye_type.markers_plus_LRI.DGI_object.RDS")
exp.m.s.cpdb.relevant <- readRDS(file = "Seurat_Objects/cpdb_results.clean.RDS")

## Color Objects
full_set.colors          <- readRDS(file = "Seurat_Objects/full_set.colors.RDS")
M.int_refined.pop.colors <- readRDS(file = "Seurat_Objects/M.int_refined.pop.colors.RDS")
