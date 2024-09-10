##=========================================================================================================================
## Load gene expression and antibody capture data sets
##=========================================================================================================================
dir.create("Seurat_Objects", showWarnings = F)

# Fetch sample dirs
all.sample.dirs <- list.files(path = "~/Data/scRNA_of_Human_Plaques/10X_Data_3_plaques_and_PBMCs/raw/raw_files/", include.dirs = T, pattern = "f00", full.names = T)

##=========================================================================================================================
## Load  samples
##=========================================================================================================================
# Load raw data
samples.raw <- list()
for(theSample in all.sample.dirs){
  cat("Loading sample:", basename(theSample), "...\n")
  
  # Load from disk
  samples.raw[[basename(theSample)]] <- Read10X(data.dir = theSample)
  
  # Fix row names in antibody assay
  rownames(samples.raw[[basename(theSample)]][["Antibody Capture"]]) <- gsub(pattern     = "_[control_]*TotalSeqC", 
                                                                                  replacement = "", 
                                                                                  x           = rownames(samples.raw[[basename(theSample)]][["Antibody Capture"]]))
}

# Turn into seurat objects
samples.seurat <- list()
for(theSample in names(samples.raw)){
  cat("Processing:", theSample, "\n")
  samples.seurat[[theSample]] <- CreateSeuratObject(counts       = samples.raw[[theSample]][["Gene Expression"]],
                                                         project      = "scTCRseq", 
                                                         min.cells    = 3,
                                                         min.features = 200)
}

# Add the antibody assay data
for(theSample in names(samples.raw)){
  cat("Processing:", theSample, "\n")
  samples.seurat[[theSample]][["ADT"]] <- CreateAssayObject(counts = samples.raw[[theSample]][["Antibody Capture"]][, colnames(samples.seurat[[theSample]])])
} 

# Add metadata
for(theSample in names(samples.raw)){ # Patient
  samples.seurat[[theSample]] <- AddMetaData(object   = samples.seurat[[theSample]], 
                                             metadata = paste("P", floor(as.numeric(substr(theSample, 11, 11)) / 2), sep = ""), 
                                             col.name = "Patient")
}

for(theSample in names(samples.raw)[c(1,3,5)]){ # PBMC Samples
  samples.seurat[[theSample]] <- AddMetaData(object   = samples.seurat[[theSample]], 
                                                  metadata = "PBMC", 
                                                  col.name = "Tissue")
} 

for(theSample in names(samples.raw)[c(2,4,6)]){ # Plaque Samples
  samples.seurat[[theSample]] <- AddMetaData(object   = samples.seurat[[theSample]], 
                                             metadata = "plaque", 
                                             col.name = "Tissue")
} 

##=========================================================================================================================
## Clean up
##=========================================================================================================================
rm(list = ls()[!ls() %in% c("samples.seurat")])
save.image()

# Save the raw objects
saveRDS(samples.seurat,   "Seurat_Objects/main.seurat.raw.RDS")
#samples.seurat <- readRDS(file = "Seurat_Objects/main.seurat.raw.RDS")
