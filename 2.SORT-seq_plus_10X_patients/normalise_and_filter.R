#=========================================================================================================================
## Normalise the Patient objects
##========================================================================================================================
dir.create("norm_QC",   showWarnings = F)

## Merge the tissue types per patient
patients.seurat <- list()
patients.seurat[["P1"]] <- merge(samples.seurat[[1]], samples.seurat[[2]], add.cell.ids = c("PBMC", "plaque"), project = "P1")
patients.seurat[["P2"]] <- merge(samples.seurat[[3]], samples.seurat[[4]], add.cell.ids = c("PBMC", "plaque"), project = "P2")
patients.seurat[["P3"]] <- merge(samples.seurat[[5]], samples.seurat[[6]], add.cell.ids = c("PBMC", "plaque"), project = "P3")

##==========================================================================================================
## Clean up, normalize, and cluster until we are satisfied
##  BASE: run an unfiltered analysis as base
## Normalize
# run sctransform
patients.seurat <- lapply(patients.seurat, function(x)SCTransform(x, verbose = T))

# Run PCA and UMAP embedding
patients.seurat <- lapply(patients.seurat, function(x)RunPCA( x, verbose = T))
patients.seurat <- lapply(patients.seurat, function(x)RunUMAP(x, dims = 1:30, verbose = T))

# Show some stats
patients.seurat

## Clustering
set.seed(1)
patients.seurat <- lapply(patients.seurat, function(x)FindNeighbors(x, dims = 1:30, verbose = T, force.recalc = T))
patients.seurat <- lapply(patients.seurat, function(x)FindClusters( x, verbose = T, resolution = 1, random.seed = 666))

for(theSample in names(patients.seurat)){
  DimPlot(patients.seurat[[theSample]], label = TRUE, pt.size = 2, label.size = 10) + NoLegend() +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  
  ggsave(paste("norm_QC/", theSample, " Unfiltered Round UMAP clusters.pdf", sep = ""), width = 15, height = 15)

  ## Store mitochondrial percentage in object meta data
  patients.seurat[[theSample]] <- PercentageFeatureSet(patients.seurat[[theSample]], pattern = "^MT-", col.name = "percent.mt")
  VlnPlot(object = patients.seurat[[theSample]], features = c("percent.mt"), group.by = "Patient")
  ggsave(paste("norm_QC/", theSample, " Unfiltered Round percent.mt.pdf", sep = ""))
  
  VlnPlot(object = patients.seurat[[theSample]], group.by = "Patient", features = c("nFeature_RNA","nCount_RNA"), y.max = 15000) + 
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold"), 
          aspect.ratio     = 1
    )
  ggsave(paste("norm_QC/", theSample, " Unfiltered Round number of transcripts.pdf", sep = ""), width = 20, height = 20)
  
  ## Check apoptotic markers
  VlnPlot(object = patients.seurat[[theSample]], features = c("KCNQ1OT1","UGDH-AS1", "GHET1"), ncol = 3)+
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  
  ggsave(paste("norm_QC/", theSample, " Unfiltered Round Apoptotic markers and lncRNAs.pdf", sep = ""), width = 15, height = 15)
}

# Store apop gene percentage in object meta data
patients.seurat <- lapply(patients.seurat, function(x)PercentageFeatureSet(x, pattern = "KCNQ1OT1", col.name = "percent.KCNQ1OT1"))
patients.seurat <- lapply(patients.seurat, function(x)PercentageFeatureSet(x, pattern = "UGDH-AS1", col.name = "percent.UGDH.AS1"))
patients.seurat <- lapply(patients.seurat, function(x)PercentageFeatureSet(x, pattern = "GHET1",    col.name = "percent.GHET1"))

for(theSample in names(patients.seurat)){
  VlnPlot(object = patients.seurat[[theSample]], group.by = "Patient", features = c("percent.KCNQ1OT1", "percent.UGDH.AS1", "percent.GHET1"))
  ggsave(paste("norm_QC/", theSample, " Percentage Apoptotic markers.pdf", sep = ""))
}

# Save the object
saveRDS(patients.seurat, "Seurat_Objects/main.patients.seurat.prefilter.RDS")
#patients.seurat <- readRDS(file = "Seurat_Objects/main.patients.seurat.prefilter.RDS")


#==========================================================================================================
##=========================================================================================================
##  First round of clustering optimisation: remove cells with less than 200 or more than 5000 genes, mito reads > 5%, apop marker reads > 2%, and (long) noncoding (RNA) genes and other junk.
##  Also regress out cell cycle, % mt, % apop genes
# First, keep only cells with less then 5% mt reads, and less than 2% of apoptotic reads.
# Next, regress out % mt, % apop genes to make a fresh SCT assay
# Then, score cell cycle gene programs.
# Finally, regress out cell cycle and all the other variables from the remaining cells.

## Subset out the mito and apop cells
patients.seurat <- lapply(patients.seurat, function(x)subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & percent.KCNQ1OT1 < 2 & percent.UGDH.AS1 < 2 & percent.GHET1 < 2))
patients.seurat

## Kick some non-coding stuff
# Identify bad (long noncoding RNA) genes
kick.list <- vector()
for(theSample in names(patients.seurat)){
  kick.list <- c(kick.list, row.names(patients.seurat[[theSample]])[grep("^RNA|^RPL|^LOC|^LINC|^RPS|-AS1$|^EEF1A|^EEF1B|^FLJ|^MRPS|^MRPL", row.names(patients.seurat[[theSample]]))])
}

# Add some extra known ones by hand
kick.list <- c(kick.list, "KCNQ1OT1", "KCNA3", "GHET1", "SH3BGRL3", "NEAT1")

# We only need each gene once
kick.list <- unique(kick.list)

## Kick all the bad genes
for(theSample in names(patients.seurat)){
  # Invert the list to all genes we want to keep
  no.kick.list <- row.names(patients.seurat[[theSample]])[!row.names(patients.seurat[[theSample]]) %in% kick.list]
  
  # Add the ADT names or we lose the whole assay
  no.kick.list <- c(no.kick.list, row.names(patients.seurat[[theSample]]@assays$ADT))
  
  # Subset to keep this inverted list
  patients.seurat[[theSample]] <- subset(patients.seurat[[theSample]], features = no.kick.list)
}

# Set up regression variables
vars.to.regress <- c("percent.mt", "percent.KCNQ1OT1", "percent.UGDH.AS1", "percent.GHET1")

# Run sctransform to re-balance the SCT assay
patients.seurat <- lapply(patients.seurat, function(x)SCTransform(x, vars.to.regress = vars.to.regress, verbose = T))

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Score the cell cycle programs
patients.seurat <- lapply(patients.seurat, function(x)CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay = "SCT"))

patients.seurat <- lapply(patients.seurat, function(x)RunPCA(x, features = c(s.genes, g2m.genes), assay = "SCT"))
for(theSample in names(patients.seurat)){
  DimPlot(patients.seurat[[theSample]], group.by = "Phase") +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  ggsave(paste("norm_QC/", theSample, " First Round UMAP_cell_cylce.pdf", sep = ""))
}

## Normalize, scale, and recluster the subsetted and cell cycle scored object
# Set up regression variables
vars.to.regress   <- c(vars.to.regress, "S.Score", "G2M.Score")

# run sctransform
patients.seurat <- lapply(patients.seurat, function(x)SCTransform(x, vars.to.regress = vars.to.regress, verbose = T))

# Run PCA and UMAP embedding
patients.seurat <- lapply(patients.seurat, function(x)RunPCA(x, verbose = T))
patients.seurat <- lapply(patients.seurat, function(x)RunUMAP(x, dims = 1:30, verbose = T))

# Show some stats
patients.seurat

## Clustering
set.seed(1)
patients.seurat <- lapply(patients.seurat, function(x)FindNeighbors(x, dims = 1:30, verbose = T, force.recalc = T))
patients.seurat <- lapply(patients.seurat, function(x)FindClusters(x, verbose = T, resolution = 1, random.seed = 666))

for(theSample in names(patients.seurat)){
  DimPlot(patients.seurat[[theSample]], label = TRUE, pt.size = 2, label.size = 10, shuffle = T) + NoLegend() +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  ggsave(paste("norm_QC/", theSample, " First Round UMAP clusters.pdf", sep = ""))
}

## Have a look at the antibody capture data
for(theSample in names(patients.seurat)){
  DefaultAssay(patients.seurat[[theSample]]) <- "ADT"
}

# Normalise
patients.seurat <- lapply(patients.seurat, function(x)NormalizeData(x, normalization.method = "CLR", margin = 2, verbose = T))

# Plot all
  for(theSample in names(patients.seurat)){
  FeaturePlot(patients.seurat[[theSample]], row.names(patients.seurat[[theSample]]@assays$ADT), pt.size = 1, cols = c("lightgrey", "darkgreen"), order = T) +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
    ggsave(paste("norm_QC/", theSample, " First Round Antibody capture feature plots.pdf", sep = ""), width = 15, height = 15)
  
  ## Plot antibody overlaps
  # CD3 vs CD4
  FeaturePlot(patients.seurat[[theSample]], row.names(patients.seurat[[theSample]]@assays$ADT)[c(1,2)], pt.size = 1, blend = T, order = T) +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  ggsave(paste("norm_QC/", theSample, " First Round Antibody capture CD3 vs CD4.pdf", sep = ""), width = 20, height = 5)
  
  # CD3 vs CD8
  FeaturePlot(patients.seurat[[theSample]], row.names(patients.seurat[[theSample]]@assays$ADT)[c(1,3)], pt.size = 1, blend = T, order = T) +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  ggsave(paste("norm_QC/", theSample, " First Round Antibody capture CD3 vs CD8.pdf", sep = ""), width = 20, height = 5)
  
  # CD4 vs CD8
  FeaturePlot(patients.seurat[[theSample]], row.names(patients.seurat[[theSample]]@assays$ADT)[c(2,3)], pt.size = 1, blend = T, order = T) +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  ggsave(paste("norm_QC/", theSample, " First Round Antibody capture CD4 vs CD8.pdf", sep = ""), width = 20, height = 5)
  
  # CD3 vs CD14
  FeaturePlot(patients.seurat[[theSample]], row.names(patients.seurat[[theSample]]@assays$ADT)[c(1,4)], pt.size = 1, blend = T, order = T) +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  ggsave(paste("norm_QC/", theSample, " First Round Antibody capture CD3 vs CD14.pdf", sep = ""), width = 20, height = 5)
  
  # CD4 vs CD14
  FeaturePlot(patients.seurat[[theSample]], row.names(patients.seurat[[theSample]]@assays$ADT)[c(2,4)], pt.size = 1, blend = T, order = T) +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  ggsave(paste("norm_QC/", theSample, " First Round Antibody capture CD4 vs CD14.pdf", sep = ""), width = 20, height = 5)
  
  # CD8 vs CD14
  FeaturePlot(patients.seurat[[theSample]], row.names(patients.seurat[[theSample]]@assays$ADT)[c(3,4)], pt.size = 1, blend = T, order = T) +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  ggsave(paste("norm_QC/", theSample, " First Round Antibody capture CD8 vs CD14.pdf", sep = ""), width = 20, height = 5)
}

# Reset the default assay
for(theSample in names(patients.seurat)){
  DefaultAssay(patients.seurat[[theSample]]) <- "SCT"
}
## Save the object
saveRDS(patients.seurat, "Seurat_Objects/main.patients.seurat.filtered.RDS")
#patients.seurat <- readRDS(file = "Seurat_Objects/main.patients.seurat.filtered.RDS")

#=========================================================================================================================
## Clean up
##========================================================================================================================
rm(list = ls()[!ls() %in% c("patients.seurat", "vars.to.regress")])
source("code/functions.R")
save.image()
