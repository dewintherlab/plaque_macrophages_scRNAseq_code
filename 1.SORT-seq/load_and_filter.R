#Load in the data
#Setup the matrixes and perform initial filtering
#-----------------------------------------------------------------------------------------
# Package install code
# install.packages("BiocManager")
# BiocManager::install(c("Seurat", "org.Hs.eg.db", "plyr","dplyr","glmGamPoi", "SeuratObject","escape","ggplot2","extrafont","pheatmap","reshape2","RColorBrewer", "SingleCellExperiment","genefilter","sctransform","viridis","ggpubr","ggsci","xlsx"))
# install.packages("devtools")


library(Seurat)
library(SeuratObject)
library(org.Hs.eg.db)
library(plyr)
library(dplyr)
library(ggplot2)
library(extrafont)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(SingleCellExperiment)
library(genefilter)
library(sctransform)
library(BiocParallel)
library(viridis)
library(glmGamPoi)
library(escape)
library(ggpubr)
library(ggsci)
library(Matrix)
library(xlsx)
source("code/functions.R")

#Get Arial font (skip if not on macos)
font_import(pattern = "Arial Unicode.ttf", prompt = F, recursive = F, paths = "/Library/Fonts/")
loadfonts()
loadfonts(device = "postscript")

#Set number of cells / plate.
n = 384

#-----------------------------------------------------------------------------------------
## prep meta data and files
# load meta data
meta.data <- read.table('raw_data/metadata.txt', header = T, skipNul = T, fill = T, stringsAsFactors = F)
# rename first column
names(meta.data)[1] <- "Patient"
# add ID column
meta.data$ID <- paste(meta.data$Patient, ".P",meta.data$Plate, sep = "")
# transform patient colum from integer to character
meta.data$Patient <- as.character(meta.data$Patient)

sum(duplicated(meta.data$ID))

# check files
all.files <- list.files(pattern="*.tsv", path = "raw_data/")
all.files
all.files %in% meta.data$File
all.files <- all.files[all.files %in% meta.data$File]
length(all.files) == dim(meta.data)[1]

#-----------------------------------------------------------------------------------------
## read in and process
seurat.object.list <- list()
unwanted.genes <- c("^ERCC", "^MALAT1$")

for (df in seq_along(all.files)) {
  cat(c("Current file: ",all.files[df],"\n"))
  current.df <- read.delim(paste("raw_data", all.files[df], sep = "/"))
  
  # fix gene names and rownames
  cat("Remove chrom info and fix row names...\n")
  current.df$GENEID <- sub("__.*", "" ,  current.df$GENEID)
  current.df <- current.df[c(!duplicated(current.df$GENEID)),]
  rownames(current.df) <- current.df$GENEID
  current.df <- current.df[,-1]
  
  # remove unwanted genes
  cat("Remove unwanted genes...\n")
  current.df <- current.df[grep(paste(unwanted.genes, collapse = "|"),rownames(current.df),invert=TRUE),]
  
  # match AE and plate number to file
  cat("Match metadata to file...\n")
  current.meta.data <- meta.data[meta.data$File == all.files[df],]
  colnames(current.df) <- paste(current.meta.data$ID, ".", 1:length(current.df), sep = "")
  
  # update gene names
  cat("Updating symbols...\n")
  current.df <- update.symbols(current.df, n)
  
  ## create object
  cat("Creating object...\n")
  current.object <- CreateSeuratObject(current.df, project = current.meta.data$ID)
  # add corresponding meta data
  cat("Adding metadata...\n")
  current.object@meta.data <- cbind(current.object@meta.data, current.meta.data)
  # add to list
  cat("Adding to list...\n")
  seurat.object.list[[df]] <- current.object
}
