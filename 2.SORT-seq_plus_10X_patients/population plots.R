#======================================================================
## Make pretty plots of the initital populations
#======================================================================
## 10X data 
## First need to rename the populations
celltypes.10X <- refined.patients.seurat$predicted.celltype.l2
unique(celltypes.10X)

# T cells
celltypes.10X[grep("CD4", celltypes.10X)]  <- "T Cells"
celltypes.10X[grep("CD8", celltypes.10X)]  <- "T Cells"
celltypes.10X[grep("Treg", celltypes.10X)] <- "T Cells"
celltypes.10X[grep("dnT", celltypes.10X)]  <- "T Cells"
celltypes.10X[grep("gdT", celltypes.10X)]  <- "T Cells"
celltypes.10X[grep("MAIT", celltypes.10X)] <- "T Cells"

# B cells
celltypes.10X[grep("B ", celltypes.10X)] <- "B Cells"

# NK Cells
celltypes.10X[grep("NK", celltypes.10X)] <- "NK Cells"

# Myeloid cells
celltypes.10X[grep("CD14", celltypes.10X)] <- "Myeloid Cells"
celltypes.10X[grep("CD16", celltypes.10X)] <- "Myeloid Cells"
celltypes.10X[grep("DC", celltypes.10X)]   <- "Myeloid Cells"

# Junk
celltypes.10X[grep("Eryth", celltypes.10X)]    <- "Erythroid Cells"
celltypes.10X[grep("ILC", celltypes.10X)]      <- "Innate Lymphoid Cells"
celltypes.10X[grep("Platelet", celltypes.10X)] <- "Platelets"

# Add back to the seurat object
refined.patients.seurat$celltypes.10X <- celltypes.10X

# Define some colors for these pops
unique(celltypes.10X)

# T
celltypes.10X.cols <- "orange3"

# B
celltypes.10X.cols <- c(celltypes.10X.cols, "purple")

# Plasmablast
celltypes.10X.cols <- c(celltypes.10X.cols, "azure2")

# NK
celltypes.10X.cols <- c(celltypes.10X.cols, "mediumslateblue")

# Myeloid
celltypes.10X.cols <- c(celltypes.10X.cols, "dodgerblue2")

# HSPC
celltypes.10X.cols <- c(celltypes.10X.cols, "azure3")

# Ery
celltypes.10X.cols <- c(celltypes.10X.cols, "red4")

# Doublet
celltypes.10X.cols <- c(celltypes.10X.cols, "azure")

# Platelets
celltypes.10X.cols <- c(celltypes.10X.cols, "azure4")

# ILC
celltypes.10X.cols <- c(celltypes.10X.cols, "bisque1")

# Name the pops
names(celltypes.10X.cols) <- unique(celltypes.10X)


## And draw the plot
customUMAP(object = refined.patients.seurat, 
           pt.size = 0.5,
           group.by = "celltypes.10X",
           title = "10X patients", cols = celltypes.10X.cols,
           shuffle = T, 
           legend.pos = "right", 
           seed = 666, plot.width = 10, plot.height = 10,
           file.name = "integration_QC/Integration UMAP 10X.pdf")


##====================================================================================
## CEL-seq data
## First need to rename the populations
celltypes.43p <- as.vector(Idents(full.43p.seurat))
unique(celltypes.43p)

# T cells
celltypes.43p[grep("T ", celltypes.43p)] <- "T Cells"
celltypes.43p[grep("FOXP3", celltypes.43p)]  <- "T Cells"

# B cells
celltypes.43p[grep("CD79", celltypes.43p)] <- "B Cells"

# NK Cells
celltypes.43p[grep("NK", celltypes.43p)] <- "NK Cells"

# Dendritic cells
celltypes.43p[grep("Dendritic", celltypes.43p)] <- "Dendritic Cells"

# Mastocytes
celltypes.43p[grep("KIT", celltypes.43p)] <- "Mastocytes"

# Macrophages
celltypes.43p[grep("CD68", celltypes.43p)] <- "Macrophages"

# SMC
celltypes.43p[grep("ACTA", celltypes.43p)] <- "Smooth Muscle Cells"

# Endothelial cells II
celltypes.43p[grep("Endothelial Cells II", celltypes.43p)] <- "Endothelial Cells II"

# Endothelial cells I
celltypes.43p[grep("Endothelial Cells I$", celltypes.43p)] <- "Endothelial Cells I"

# Add back to the seurat object
full.43p.seurat$celltypes.43p <- celltypes.43p

# Define some colors for these pops
unique(celltypes.43p)

# T
celltypes.43p.cols <- "orange3"

# Endothelial I
celltypes.43p.cols <- c(celltypes.43p.cols, "chartreuse4")

# NK
celltypes.43p.cols <- c(celltypes.43p.cols, "mediumslateblue")

# Macrophages
celltypes.43p.cols <- c(celltypes.43p.cols, "dodgerblue2")

# SMC
celltypes.43p.cols <- c(celltypes.43p.cols, "chartreuse2")

# Endothelial II
celltypes.43p.cols <- c(celltypes.43p.cols, "chartreuse3")

# DC
celltypes.43p.cols <- c(celltypes.43p.cols, "cadetblue2")

# B
celltypes.43p.cols <- c(celltypes.43p.cols, "purple")

# Mastocytes
celltypes.43p.cols <- c(celltypes.43p.cols, "dodgerblue")


# Name the pops
names(celltypes.43p.cols) <- unique(celltypes.43p)


## And draw the plot
customUMAP(object = full.43p.seurat, 
           pt.size = 0.5,
           group.by = "celltypes.43p",
           title = "43p patients", cols = celltypes.43p.cols,
           shuffle = T, 
           legend.pos = "right", 
           seed = 666, plot.width = 10, plot.height = 10,
           file.name = "integration_QC/Integration UMAP CEL-seq.pdf")


## Maker some TNFR plots
dir.create("cellphonedb_results/TNFRSF", showWarnings = F, recursive = T)

# Check which TNFRs we got
TNFRs <- row.names(from_full.integrated.mye.seurat)[grep("^TNFRSF", row.names(from_full.integrated.mye.seurat))]

# And plot them
for(theGene in TNFRs){
  stratifyByExpression(object    = from_full.integrated.mye.seurat, 
                       strat.by  = theGene, 
                       file.name = paste("cellphonedb_results/TNFRSF/", theGene, sep = ""), 
                       verbose   = F, 
                       onlyUMAP  = T)
}

# Plot ITGA4 in the same manner
stratifyByExpression(object    = from_full.integrated.mye.seurat, 
                     strat.by  = "ITGA4", 
                     file.name = paste("cellphonedb_results/TNFRSF/", "ITGA4", sep = ""), 
                     verbose   = F, 
                     onlyUMAP  = T)


full.43p.seurat$Phenotype
full.43p.seurat
refined.patients.seurat
integrated.mye.seurat

#=================================================================================================
## Some random plots
dir.create("various_plots", showWarnings = F)

bunchOfCustomPlots(object = from_full.integrated.mye.seurat, features = c("CSF1", "CSF1R"), Vln.draw.names = T, name = "various_plots/CSF1_CSF1R", Vln.color = M.int_refined.pop.colors )
bunchOfCustomPlots(object = from_full.integrated.mye.seurat, features = c("IL1B", "IL1R1"), Vln.draw.names = T, name = "various_plots/IL1B_IL1R1", Vln.color = M.int_refined.pop.colors )

# ANRIL
from_full.integrated.mye.seurat$RNA@counts[grep("CDKN", row.names(from_full.integrated.mye.seurat@assays$RNA@counts)),1:5, drop = F]
from_full.integrated.mye.seurat$RNA@counts[grep("ANRIL", row.names(from_full.integrated.mye.seurat@assays$RNA@data)),1:5, drop = F]

#bunchOfCustomPlots(object = from_full.integrated.mye.seurat, features = c("CDKN2B-AS"),         Vln.draw.names = T, name = "various_plots/ANRIL", Vln.color = M.int_refined.pop.colors)

#=================================================================================================
## Protein panel plots
dir.create("various_plots/AE_measured_proteins", showWarnings = F, recursive = T)

# Load the list of proteins measured in AE
AE.prot.list <- read.table(file = "raw_data/proteins_measured.txt", header = F, col.names = c("Protein", "Gene"))
AE.prot.list

# First plot globally to see if mac specific
for(theGene in AE.prot.list$Gene){
  stratifyByExpression(object    = integrated.full.seurat,
                       onlyUMAP = T,
                       strat.by  = theGene,
                       file.name = paste("various_plots/AE_measured_proteins/", AE.prot.list[which(AE.prot.list$Gene == theGene), "Protein"], "stratified.full_set"))
}

# Short list of genes potentially relevant to macs (deduced from global plots)
AE.prot.short.list <- AE.prot.list[c(1,2,3,12),]
AE.prot.short.list  

# Plot these in the myeloid clusters 
for(theGene in AE.prot.short.list$Gene){
  stratifyByExpression(object    = from_full.integrated.mye.seurat,
                       onlyUMAP  = F,
                       strat.by  = theGene,
                       file.name = paste("various_plots/AE_measured_proteins/", AE.prot.short.list[which(AE.prot.short.list$Gene == theGene), "Protein"], "stratified.myeloid_pops"))
}
  
# Boot MPO form the short list, mostly expressed in monocytes
AE.prot.short.list <- AE.prot.short.list[1:3,]

# Plot expression levels in the individual populations
bunchOfCustomPlots(object         = from_full.integrated.mye.seurat, 
                   features       = AE.prot.short.list$Gene, 
                   assay          = "RNA",
                   Vln.draw.names = F,
                   Vln.color      = M.int_refined.pop.colors, Vln.width = 15, Vln.height = 5,
                   ncol           = 3,
                   name           = "various_plots/AE_measured_proteins/shortlist.myeloid")

# Plot by cell cycle phase
bunchOfCustomPlots(object         = from_full.integrated.mye.seurat, 
                   features       = AE.prot.short.list$Gene, 
                   group.by       = "Phase", 
                   assay          = "RNA",
                   Vln.draw.names = T,
                   ncol           = 3,
                   name           = "various_plots/AE_measured_proteins/shortlist.myeloid.phase")
  
## Set up archetypes
archetypes <- as.vector(Idents(from_full.integrated.mye.seurat))
archetypes[grep("Monocytes", archetypes)]    <- "Monocytes"
archetypes[grep("Inflammatory", archetypes)] <- "Inflammatory"
archetypes[grep("Resident", archetypes)]     <- "Resident"
archetypes[grep("Foamy", archetypes)]        <- "Foamy"
archetypes[grep("Lipid", archetypes)]        <- "Resident"

# Sanity check
levels(as.factor(archetypes))

# Add to the seurat object
from_full.integrated.mye.seurat <- AddMetaData(from_full.integrated.mye.seurat,metadata = archetypes, col.name = "archetype")

# Set up colors
archetype.colors <- c("Resident" = "#528347", "Inflammatory" = "#87482F", "Foamy" = "#9476A3", "Monocytes" = "#6495ED")

# Plot by archetype
bunchOfCustomPlots(object         = from_full.integrated.mye.seurat, 
                   features       = AE.prot.short.list$Gene, 
                   group.by       = "archetype", 
                   assay          = "RNA",
                   Vln.color      = archetype.colors,
                   Vln.draw.names = T,
                   ncol           = 3,
                   name           = "various_plots/AE_measured_proteins/shortlist.myeloid.archetype")


#=================================================================================================
## OLINK panels available for AE patients
dir.create("various_plots/OLINK_panels", showWarnings = F, recursive = T)

## Get the genes in the panels
# Read the list of proteins in the panel
olink.panels.prot.list <- scan("raw_data/olink_protein_panels.txt", what = "list")
olink.panels.prot.list <- toupper(olink.panels.prot.list)

# Annotate with gene info
olink.panels.gene.list  <- AnnotationDbi::select(org.Hs.eg.db, keys = olink.panels.prot.list, keytype = "UNIPROT", columns = "SYMBOL", multiVals = "first")

# Fix IGLC2 which doesn't get assigned automatically for some reason
olink.panels.gene.list[grep("P0CG05", olink.panels.gene.list$UNIPROT),"SYMBOL"] <- "IGLC2"

# Clean up
olink.panel.gene.prot.key <- olink.panels.gene.list
olink.panels.gene.list    <- sort(olink.panels.gene.list$SYMBOL)
d                         <- duplicated(olink.panels.gene.list)
olink.panels.gene.list    <- olink.panels.gene.list[!d]

## Make some general plots
## Population heatmap
# Get average expression per cluster
avex.mye <- as.data.frame(AverageExpression(integrated.mye.seurat, assays = "RNA", slot = "data", features = olink.panels.gene.list))

# Keep only expressed genes
avex.mye <- avex.mye[!rowSums(avex.mye) == 0,]

# Prep the annotations
colnames(avex.mye)     <- gsub("RNA\\.", "", colnames(avex.mye))
hm.pop.anno            <- data.frame("Population" = colnames(avex.mye))
rownames(hm.pop.anno)  <- hm.pop.anno$Population
hm.pop.anno.col        <- M.int_refined.pop.colors[c(4,5,6,7,11,12,13,8,9,10)]
hm.pop.anno.col        <- list("Population" = hm.pop.anno.col)
names(hm.pop.anno.col$Population) <- hm.pop.anno$Population

# And plot
pheatmap(avex.mye, 
         show_colnames = F, 
         show_rownames = F, 
         scale = "row", 
         annotation_col = hm.pop.anno, 
         annotation_colors = hm.pop.anno.col,
         filename = "various_plots/OLINK_panels/population_heatmap.pdf", width = 10, height = 20)

set.seed(666)
pheatmap(avex.mye, 
         show_colnames = F, 
         show_rownames = T, 
         scale = "row", 
         kmeans_k = 8,
         annotation_col = hm.pop.anno, 
         annotation_colors = hm.pop.anno.col,
         filename = "various_plots/OLINK_panels/population_heatmap_kmeans.pdf", width = 10)

set.seed(666)
hm <- pheatmap(avex.mye,
               scale = "row", 
               kmeans_k = 8,
               silent = T)

# Extract per archetype
olink.inf.genes  <- names(hm$kmeans$cluster)[hm$kmeans$cluster == 1 | hm$kmeans$cluster == 2 | hm$kmeans$cluster == 4]
olink.foam.genes <- names(hm$kmeans$cluster)[hm$kmeans$cluster == 3 | hm$kmeans$cluster == 5]
olink.res.genes  <- names(hm$kmeans$cluster)[hm$kmeans$cluster == 6 | hm$kmeans$cluster == 7 | hm$kmeans$cluster == 8]

# Confirm in the seurat object
DoHeatmap(object = integrated.mye.seurat, features = c(olink.inf.genes, olink.foam.genes, olink.res.genes), label = F, raster = F, group.colors = M.int_refined.pop.colors)
ggsave("various_plots/OLINK_panels/population_heatmap.singlecell.pdf", width = 10)

# Let's take the top 10 per group
avex.mye.archetype <- data.frame("Inflammatory" = apply(avex.mye[,1:4], 1, mean), "Foamy" = apply(avex.mye[,8:10], 1, mean), "Resident" = apply(avex.mye[,5:7], 1, mean))

olink.inf.genes.top10  <- row.names(avex.mye.archetype[olink.inf.genes, ][order(avex.mye.archetype[olink.inf.genes, ]$Inflammatory, decreasing = T),])[1:10]
olink.foam.genes.top10 <- row.names(avex.mye.archetype[olink.foam.genes,][order(avex.mye.archetype[olink.foam.genes,]$Foamy,        decreasing = T),])[1:10]
olink.res.genes.top10  <- row.names(avex.mye.archetype[olink.res.genes, ][order(avex.mye.archetype[olink.res.genes, ]$Resident,     decreasing = T),])[1:10]

# Check the heatmap
DoHeatmap(object = integrated.mye.seurat, features = c(olink.inf.genes.top10, olink.foam.genes.top10, olink.res.genes.top10), label = F, raster = F, group.colors = M.int_refined.pop.colors)
ggsave("various_plots/OLINK_panels/population_heatmap_top10.singlecell.pdf", width = 10)

# Make the rest of the plots
bunchOfCustomPlots(object = integrated.mye.seurat, features = olink.inf.genes.top10,  assay = "RNA", Vln.draw.names = F, Vln.color = M.int_refined.pop.colors, ncol = 3, name = "various_plots/OLINK_panels/inf_pops_top10")
bunchOfCustomPlots(object = integrated.mye.seurat, features = olink.foam.genes.top10, assay = "RNA", Vln.draw.names = F, Vln.color = M.int_refined.pop.colors, ncol = 3, name = "various_plots/OLINK_panels/foam_pops_top10")
bunchOfCustomPlots(object = integrated.mye.seurat, features = olink.res.genes.top10,  assay = "RNA", Vln.draw.names = F, Vln.color = M.int_refined.pop.colors, ncol = 3, name = "various_plots/OLINK_panels/res_pops_top10")

## Set up archetypes
archetypes <- as.vector(Idents(integrated.mye.seurat))
archetypes[grep("Monocytes",    archetypes)] <- "Monocytes"
archetypes[grep("Inflammatory", archetypes)] <- "Inflammatory"
archetypes[grep("Resident",     archetypes)] <- "Resident"
archetypes[grep("Foamy",        archetypes)] <- "Foamy"
archetypes[grep("Lipid",        archetypes)] <- "Resident"

integrated.mye.seurat <- AddMetaData(integrated.mye.seurat, metadata = archetypes, col.name = "archetype")

bunchOfCustomPlots(object = integrated.mye.seurat, group.by = "archetype", features = olink.inf.genes.top10,  assay = "RNA", Vln.draw.names = T, Vln.color = archetype.colors, ncol = 3, name = "various_plots/OLINK_panels/inf_types_top10")
bunchOfCustomPlots(object = integrated.mye.seurat, group.by = "archetype", features = olink.foam.genes.top10, assay = "RNA", Vln.draw.names = T, Vln.color = archetype.colors, ncol = 3, name = "various_plots/OLINK_panels/foam_types_top10")
bunchOfCustomPlots(object = integrated.mye.seurat, group.by = "archetype", features = olink.res.genes.top10,  assay = "RNA", Vln.draw.names = T, Vln.color = archetype.colors, ncol = 3, name = "various_plots/OLINK_panels/res_types_top10")

## Check furhter down the chain for the inf group
dir.create("various_plots/OLINK_panels/inf_top20", showWarnings = F, recursive = T)
olink.inf.genes.top20  <- row.names(avex.mye.archetype[olink.inf.genes, ][order(avex.mye.archetype[olink.inf.genes, ]$Inflammatory, decreasing = T),])[11:20]
bunchOfCustomPlots(object = integrated.mye.seurat, features = olink.inf.genes.top20,  assay = "RNA", Vln.draw.names = F, Vln.color = M.int_refined.pop.colors, ncol = 3, name = "various_plots/OLINK_panels/inf_top20/inf_pops_top20")
bunchOfCustomPlots(object = from_full.integrated.mye.seurat, features = olink.inf.genes.top20,  assay = "RNA", Vln.draw.names = F, Vln.color = M.int_refined.pop.colors, ncol = 3, name = "various_plots/OLINK_panels/inf_top20/inf_pops_mono_top20")

dir.create("various_plots/OLINK_panels/inf_top30", showWarnings = F, recursive = T)
olink.inf.genes.top30  <- row.names(avex.mye.archetype[olink.inf.genes, ][order(avex.mye.archetype[olink.inf.genes, ]$Inflammatory, decreasing = T),])[21:30]
bunchOfCustomPlots(object = integrated.mye.seurat, features = olink.inf.genes.top30,  assay = "RNA", Vln.draw.names = F, Vln.color = M.int_refined.pop.colors, ncol = 3, name = "various_plots/OLINK_panels/inf_top30/inf_pops_top30")
bunchOfCustomPlots(object = from_full.integrated.mye.seurat, features = olink.inf.genes.top30,  assay = "RNA", Vln.draw.names = F, Vln.color = M.int_refined.pop.colors, ncol = 3, name = "various_plots/OLINK_panels/inf_top30/inf_pops_mono_top30")

dir.create("various_plots/OLINK_panels/inf_top40", showWarnings = F, recursive = T)
olink.inf.genes.top40  <- row.names(avex.mye.archetype[olink.inf.genes, ][order(avex.mye.archetype[olink.inf.genes, ]$Inflammatory, decreasing = T),])[31:40]
bunchOfCustomPlots(object = integrated.mye.seurat, features = olink.inf.genes.top40,  assay = "RNA", Vln.draw.names = F, Vln.color = M.int_refined.pop.colors, ncol = 3, name = "various_plots/OLINK_panels/inf_top40/inf_pops_top40")
bunchOfCustomPlots(object = from_full.integrated.mye.seurat, features = olink.inf.genes.top40,  assay = "RNA", Vln.draw.names = F, Vln.color = M.int_refined.pop.colors, ncol = 3, name = "various_plots/OLINK_panels/inf_top40/inf_pops_mono_top40")

dir.create("various_plots/OLINK_panels/inf_top50", showWarnings = F, recursive = T)
olink.inf.genes.top50  <- row.names(avex.mye.archetype[olink.inf.genes, ][order(avex.mye.archetype[olink.inf.genes, ]$Inflammatory, decreasing = T),])[41:50]
bunchOfCustomPlots(object = integrated.mye.seurat, features = olink.inf.genes.top50,  assay = "RNA", Vln.draw.names = F, Vln.color = M.int_refined.pop.colors, ncol = 3, name = "various_plots/OLINK_panels/inf_top50/inf_pops_top50")
bunchOfCustomPlots(object = from_full.integrated.mye.seurat, features = olink.inf.genes.top50,  assay = "RNA", Vln.draw.names = F, Vln.color = M.int_refined.pop.colors, ncol = 3, name = "various_plots/OLINK_panels/inf_top50/inf_pops_mono_top50")

#=================================================================================================
## method or CD45 selection bias in myeloid clusters

# Check CD45 distribution
bunchOfCustomPlots(object = integrated.mye.seurat, features = "PTPRC", group.by = "Method", assay = "RNA", name = "various_plots/CD45_by_method.myeloid.pdf")
customVln(object = integrated.mye.seurat, features = "PTPRC", splitPlot = T, split.by = "Method", assay = "RNA", name = "various_plots/CD45_by_method.split.myeloid.pdf")
integrated.mye.seurat <- stratifyByExpression(object = integrated.mye.seurat, strat.by = "PTPRC", file.name = "various_plots/CD45_stratified", return.object = T)
stratifyByExpression(object = integrated.mye.seurat, 
                     strat.by = "PTPRC", 
                     return.object = F, 
                     file.name = "various_plots/CD45_stratified_mac_genes", 
                     gene.groups = list("Mac genes" = c("CD14", "FCGR3A", "CD68", "TREM2", "FOLR2")))

# Check CD45 and method contribution to populaiton size
no_CD45.dist     <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$PTPRC_expr == "Zero",])])
method.dist      <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$Method == "CEL-seq",])])
total.dist       <- table(Idents(integrated.mye.seurat))
CD45.ratio       <- no_CD45.dist / total.dist
method.ratio     <- method.dist  / total.dist
method.ratio.10x <- 1 - method.ratio


# Plot the ratios
m           <- merge(data.frame(method.ratio.10x), data.frame(method.ratio), by = 1)
colnames(m) <- c("Population", "Method_10X", "Method_CEL-seq")
m           <- melt(m)
colnames(m) <- c("Population", "Type", "Ratio")

CEl.cells  <- sum(integrated.mye.seurat@meta.data$Method == "CEL-seq")
tenX.cells <- sum(integrated.mye.seurat@meta.data$Method == "10X")
CEL.ratio  <- CEl.cells  / (CEl.cells + tenX.cells)
tenX.ratio <- tenX.cells / (CEl.cells + tenX.cells)

ggplot(m, aes(x = Population, y = Ratio, group = Type, col = Type)) +
  geom_line() +
  geom_hline(yintercept = CEL.ratio, col = "purple", size = 0.25) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("various_plots/method_distribution.pdf")

# Look at CD45 split by method
no_CD45_and_method.dist     <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$Method == "CEL-seq" & integrated.mye.seurat@meta.data$PTPRC_expr == "Zero",])])
no_CD45_and_method10x.dist  <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$Method == "10X" & integrated.mye.seurat@meta.data$PTPRC_expr == "Zero",])])

# Normalise to method level
no_CD45_and_method.ratio    <- no_CD45_and_method.dist    / method.dist
no_CD45_and_method10x.ratio <- no_CD45_and_method10x.dist / (total.dist - method.dist)

# Total CD45 NULL ratio
CD45_NULL.ratio <- sum(no_CD45.dist) / (CEl.cells + tenX.cells)

# Plot the ratios
m           <- merge(data.frame(no_CD45_and_method.ratio), data.frame(no_CD45_and_method10x.ratio), by = 1)
colnames(m) <- c("Population", "CD45_NULL_CEL-seq", "CD45_NULL_10X")
m           <- melt(m)
colnames(m) <- c("Population", "Type", "Ratio")

ggplot(m, aes(x = Population, y = Ratio, group = Type, col = Type)) +
  geom_line() +
  geom_hline(yintercept = CD45_NULL.ratio, col = "purple", size = 0.25) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("various_plots/CD45_NULL_by_method.pdf")


## Check CD68 in the same way as a control
# Check CD68 distribution
bunchOfCustomPlots(object = integrated.mye.seurat, features = "CD68", group.by = "Method", assay = "RNA", name = "various_plots/CD68_by_method.myeloid.pdf")
customVln(object = integrated.mye.seurat, features = "CD68", splitPlot = T, split.by = "Method", assay = "RNA", name = "various_plots/CD68_by_method.split.myeloid.pdf")
integrated.mye.seurat <- stratifyByExpression(object = integrated.mye.seurat, strat.by = "CD68", file.name = "various_plots/CD68_stratified", return.object = T)

# Check CD45 and method contribution to populaiton size
no_CD68.dist     <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$CD68_expr == "Zero",])])
method.dist      <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$Method == "CEL-seq",])])
total.dist       <- table(Idents(integrated.mye.seurat))
CD68.ratio       <- no_CD68.dist / total.dist
method.ratio     <- method.dist  / total.dist
method.ratio.10x <- 1 - method.ratio

# Look at CD45 split by method
no_CD68_and_method.dist     <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$Method == "CEL-seq" & integrated.mye.seurat@meta.data$CD68_expr == "Zero",])])
no_CD68_and_method10x.dist  <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$Method == "10X" & integrated.mye.seurat@meta.data$CD68_expr == "Zero",])])

# Normalise to method level
no_CD68_and_method.ratio    <- no_CD68_and_method.dist    / method.dist
no_CD68_and_method10x.ratio <- no_CD68_and_method10x.dist / (total.dist - method.dist)

# Total CD45 NULL ratio
CD68_NULL.ratio <- sum(no_CD68.dist) / (CEl.cells + tenX.cells)

# Plot the ratios
m           <- merge(data.frame(no_CD68_and_method.ratio), data.frame(no_CD68_and_method10x.ratio), by = 1)
colnames(m) <- c("Population", "CD68_NULL_CEL-seq", "CD68_NULL_10X")
m           <- melt(m)
colnames(m) <- c("Population", "Type", "Ratio")

ggplot(m, aes(x = Population, y = Ratio, group = Type, col = Type)) +
  geom_line() +
  geom_hline(yintercept = CD68_NULL.ratio, col = "purple", size = 0.25) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("various_plots/CD68_NULL_by_method.pdf")


## Check CD14 in the same way as a control
# Check CD14 distribution
bunchOfCustomPlots(object = integrated.mye.seurat, features = "CD14", group.by = "Method", assay = "RNA", name = "various_plots/CD14_by_method.myeloid.pdf")
customVln(object = integrated.mye.seurat, features = "CD14", splitPlot = T, split.by = "Method", assay = "RNA", name = "various_plots/CD14_by_method.split.myeloid.pdf")
integrated.mye.seurat <- stratifyByExpression(object = integrated.mye.seurat, strat.by = "CD14", file.name = "various_plots/CD14_stratified", return.object = T)

# Check CD45 and method contribution to populaiton size
no_CD14.dist     <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$CD14_expr == "Zero",])])
method.dist      <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$Method == "CEL-seq",])])
total.dist       <- table(Idents(integrated.mye.seurat))
CD14.ratio       <- no_CD14.dist / total.dist
method.ratio     <- method.dist  / total.dist
method.ratio.10x <- 1 - method.ratio

# Look at CD45 split by method
no_CD14_and_method.dist     <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$Method == "CEL-seq" & integrated.mye.seurat@meta.data$CD14_expr == "Zero",])])
no_CD14_and_method10x.dist  <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$Method == "10X" & integrated.mye.seurat@meta.data$CD14_expr == "Zero",])])

# Normalise to method level
no_CD14_and_method.ratio    <- no_CD14_and_method.dist    / method.dist
no_CD14_and_method10x.ratio <- no_CD14_and_method10x.dist / (total.dist - method.dist)

# Total CD45 NULL ratio
CD14_NULL.ratio <- sum(no_CD14.dist) / (CEl.cells + tenX.cells)

# Plot the ratios
m           <- merge(data.frame(no_CD14_and_method.ratio), data.frame(no_CD14_and_method10x.ratio), by = 1)
colnames(m) <- c("Population", "CD14_NULL_CEL-seq", "CD14_NULL_10X")
m           <- melt(m)
colnames(m) <- c("Population", "Type", "Ratio")

ggplot(m, aes(x = Population, y = Ratio, group = Type, col = Type)) +
  geom_line() +
  geom_hline(yintercept = CD14_NULL.ratio, col = "purple", size = 0.25) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("various_plots/CD14_NULL_by_method.pdf")


## Check CD16 in the same way as a control
# Check CD16 distribution
bunchOfCustomPlots(object = integrated.mye.seurat, features = "FCGR3A", group.by = "Method", assay = "RNA", name = "various_plots/CD16_by_method.myeloid.pdf")
customVln(object = integrated.mye.seurat, features = "FCGR3A", splitPlot = T, split.by = "Method", assay = "RNA", name = "various_plots/CD16_by_method.split.myeloid.pdf")
integrated.mye.seurat <- stratifyByExpression(object = integrated.mye.seurat, strat.by = "FCGR3A", file.name = "various_plots/CD16_stratified", return.object = T)

# Check CD45 and method contribution to populaiton size
no_CD16.dist     <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$FCGR3A_expr == "Zero",])])
method.dist      <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$Method == "CEL-seq",])])
total.dist       <- table(Idents(integrated.mye.seurat))
CD16.ratio       <- no_CD16.dist / total.dist
method.ratio     <- method.dist  / total.dist
method.ratio.10x <- 1 - method.ratio

# Look at CD45 split by method
no_CD16_and_method.dist     <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$Method == "CEL-seq" & integrated.mye.seurat@meta.data$FCGR3A_expr == "Zero",])])
no_CD16_and_method10x.dist  <- table(Idents(integrated.mye.seurat)[row.names(integrated.mye.seurat@meta.data[integrated.mye.seurat@meta.data$Method == "10X" & integrated.mye.seurat@meta.data$FCGR3A_expr == "Zero",])])

# Normalise to method level
no_CD16_and_method.ratio    <- no_CD16_and_method.dist    / method.dist
no_CD16_and_method10x.ratio <- no_CD16_and_method10x.dist / (total.dist - method.dist)

# Total CD45 NULL ratio
CD16_NULL.ratio <- sum(no_CD16.dist) / (CEl.cells + tenX.cells)

# Plot the ratios
m           <- merge(data.frame(no_CD16_and_method.ratio), data.frame(no_CD16_and_method10x.ratio), by = 1)
colnames(m) <- c("Population", "CD16_NULL_CEL-seq", "CD16_NULL_10X")
m           <- melt(m)
colnames(m) <- c("Population", "Type", "Ratio")

ggplot(m, aes(x = Population, y = Ratio, group = Type, col = Type)) +
  geom_line() +
  geom_hline(yintercept = CD16_NULL.ratio, col = "purple", size = 0.25) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("various_plots/CD16_NULL_by_method.pdf")

#=================================================================================================
## Some Monaco targets
dir.create("various_plots/Monaco_targets", showWarnings = F, recursive = T)
bunchOfCustomPlots(object = from_full.integrated.mye.seurat, features = c("ZEB2", "CLEC4A", "IRF4", "IRF5"), assay = "RNA", ncol = 2, Vln.draw.names = F, name = "various_plots/Monaco_targets/Monaco_targets", Vln.color = M.int_refined.pop.colors )
for(theGene in c("ZEB2", "CLEC4A", "IRF4", "IRF5")){
  stratifyByExpression(object    = from_full.integrated.mye.seurat,
                       onlyUMAP  = F,
                       strat.by  = theGene,
                       file.name = paste("various_plots/Monaco_targets/",theGene, ".stratified.myeloid_pops.pdf", sep = ""))
}

# TREM1/TREM2 dichotomy
FeaturePlot(object = integrated.mye.seurat, features = c("TREM1", "TREM2"), pt.size = 2, order = T, blend = T)
plot.cor(object = integrated.mye.seurat, features = c("TREM2", "PLIN2"), cor.feature = "TREM1", file.name = "various_plots/Monaco_targets/TREM cor.pdf")

# Fix IRF5 'not found in assay'
stratifyByExpression(object    = from_full.integrated.mye.seurat,
                     assay     = "integrated",
                     onlyUMAP  = T,
                     strat.by  = "IRF5",
                     file.name = paste("various_plots/Monaco_targets/",theGene, ".stratified.myeloid_pops.pdf", sep = ""))
                     
bunchOfCustomPlots(object = from_full.integrated.mye.seurat, features = c("MKI67"), assay = "RNA", ncol = 2, Vln.draw.names = T, name = "various_plots/MKI67", Vln.color = archetype.colors )

# Split on phenotype
customVln(object = from_full.integrated.mye.seurat, features = c("ITGAX", "IRF5"), assay = "RNA", ncol = 2, splitPlot = T, split.by = "Phenotype", group.by = "archetype", draw.names = T, name = "various_plots/Monaco_targets/Monaco_targets_ITGAX_symp_split_archetype.pdf" )
customVln(object = from_full.integrated.mye.seurat, features = c("ITGAX", "IRF5"), assay = "RNA", ncol = 2, splitPlot = T, split.by = "Phenotype", draw.names = F, name = "various_plots/Monaco_targets/Monaco_targets_ITGAX_symp_split.pdf", stack = T )

## Do MKI67
stratifyByExpression(object    = integrated.mye.seurat,
                     assay     = "integrated",
                     onlyUMAP  = T,
                     strat.by  = "MKI67",
                     file.name = paste("various_plots/","MKI67", ".stratified.myeloid_pops.pdf", sep = ""))

FeaturePlot(object = from_full.integrated.mye.seurat, features = "IRF5", )  

bunchOfCustomPlots(object = from_full.integrated.mye.seurat, features = c("ZEB2", "CLEC4A", "IRF4", "IRF5"), assay = "RNA", ncol = 2, group.by = "archetype", Vln.draw.names = T, name = "various_plots/Monaco_targets/Monaco_targets_archetype", Vln.color = archetype.colors )

bunchOfCustomPlots(object = from_full.integrated.mye.seurat, features = "TREM1", assay = "RNA", name = "various_plots/Monaco_targets/TREM1", Vln.color = M.int_refined.pop.colors)
bunchOfCustomPlots(object = integrated.mye.seurat, features = "TREM1", assay = "RNA", name = "various_plots/Monaco_targets/TREM1 macs only", Vln.color = M.int_refined.pop.colors)


## Fibrotic genes Annette
# Create a sub dir
dir.create("various_plots/Annette")

# Load the requested genes
fib.genes <- c("Col1a1","Col1a2","Col3a1","Acta2", "Fn1","Sparc","Fbln2","Serpinh1")
fib.genes <- toupper(fib.genes)

# And plot
bunchOfCustomPlots(object          = integrated.mye.seurat, 
                   features        = fib.genes, 
                   assay           = "RNA", 
                   Vln.draw.names  = F, 
                   Vln.pt.size     = 1,
                   name            = "various_plots/Annette/fib_genes", 
                   feature.pt.size = 1.25,
                   Vln.color       = M.int_refined.pop.colors,
                   ncol            = 3,
                   Vln.width       = 20)

## Also check in the 18p data
# Load the seurat object
mye.18p.seurat <- readRDS(file = "raw_data/v3.all.seur.combined.M_clusters.Rds")

# Plot
bunchOfCustomPlots(object          = mye.18p.seurat,  
                   features        = fib.genes, 
                   assay           = "RNA",
                   idents          = c("Mye.0", "Mye.1", "Mye.2"), 
                   Vln.pt.size     = 1,
                   reduction       = "tsne",
                   Vln.draw.names  = T, 
                   name            = "various_plots/Annette/fib_genes.18_patients_from_paper", 
                   feature.pt.size = 1.25,
                   ncol            = 3,
                   Vln.width       = 20)
customUMAP(object     = mye.18p.seurat, 
           pt.size    = 4,
           label      = F, 
           title      = "18 patients macrophages", 
           reduction  = "tsne", 
           cells      = WhichCells(object = mye.18p.seurat, idents = c("Mye.0", "Mye.1", "Mye.2")), 
           legend.pos = "top", 
           file.name  = "various_plots/Annette/tSNE - 18 patients form paper.pdf")

# Clean it up again
rm(mye.18p.seurat)


## Plot monocytes panels
# Make the standard plots
unique(Idents(from_full.integrated.mye.seurat))
bunchOfCustomPlots(object         = from_full.integrated.mye.seurat, 
                   idents         = unique(Idents(from_full.integrated.mye.seurat))[1:3], 
                   features       = c("CD14", "LYZ", "S100A8", "CCR2", "ITGAX", "CD86", "FCGR3A", "MS4A7", "CDKN1C"), 
                   Vln.draw.names = F, 
                   name           = "various_plots/monocyte panels", 
                   ncol           = 3,
                   assay          = "RNA", 
                   Vln.color      = M.int_refined.pop.colors, Vln.pt.size = 0, dot.scale = 20, Vln.width = 16, Vln.height = 6
)

# Draw a heatmap of the markers
mono.sep.markers.top5 <- pbmc.monocytes.seurat.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

DoHeatmap(object = from_full.integrated.mye.seurat, 
           cells = WhichCells(object = from_full.integrated.mye.seurat, idents =  unique(Idents(from_full.integrated.mye.seurat))[1:3]),
        features = mono.sep.markers.top5$gene, assay = "RNA", label = F, raster = F, group.colors = M.int_refined.pop.colors)
ggsave("various_plots/mono marker heatmap.pdf")


## Plot druggable target
d.targets <- c("AKR1B1", "ALOX5", "ITGA4", "RARA")
customVln(final.pop.call.integrated.mye.seurat, features = d.targets, group.by = "archetype", cols = archetype.colors, assay = "RNA", ncol = 2,draw.names = F, name = "various_plots/drug targets violin.pdf")
customDot(final.pop.call.integrated.mye.seurat, features = d.targets, group.by = "archetype", assay = "RNA", name = "various_plots/drug targets dot.pdf")


## Plot some markers for the 43 populations
# Check what we have
DimPlot(full.43p.seurat, group.by = "celltypes.43p", cols = celltypes.43p.cols) + NoLegend()

# Define some cell type markers
celltype.markers <- c("CD3E", "CD8A", "CD4", "CD68", "CD1C", "NCAM1", "CD79A", "ACTA2", "CD34", "KIT")

# And make some plots
bunchOfCustomPlots(object          = full.43p.seurat, 
                   group.by        = "celltypes.43p", 
                   features        = celltype.markers, 
                   Vln.draw.names  = F, 
                   Vln.pt.size     = 0, 
                   feature.pt.size = 1.5,
                   name            = "various_plots/43p celltypes", 
                   assay           = "RNA", 
                   Vln.width       = 10,
                   Vln.height      = 15,
                   Vln.color       = celltypes.43p.cols,
                   ncol            = 2,
                   Dot.width       = 10, 
                   Dot.height      = 5, 
                   dot.scale       = 15
                )

## Plot pop markers
# Get pop markers
full.43p.seurat.celltype.idents         <- full.43p.seurat
Idents(full.43p.seurat.celltype.idents) <- full.43p.seurat$celltypes.43p
full.43p.seurat.celltype.markers        <- FindAllMarkers(object = full.43p.seurat.celltype.idents, assay = "RNA", min.pct = 0.25, min.diff.pct = 0.25)
full.43p.seurat.celltype.markers.top5   <- full.43p.seurat.celltype.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

# Plot a heatmap
DoHeatmap(object       = full.43p.seurat.celltype.idents, 
          features     = full.43p.seurat.celltype.markers.top5$gene, 
          group.colors = celltypes.43p.cols,
          assay        = "RNA",
          raster       = F, 
          label        = F)
ggsave(filename = "various_plots/43p celltypes - heatmap.pdf")


# Random genes
customVln(object = integrated.mye.seurat, features = "FABP4", splitPlot = F, assay = "RNA", name = "various_plots/FABP4.pdf")

# ERVs
bunchOfCustomPlots(object = final.pop.call.integrated.mye.seurat, 
                   features = c("HERVK", "LTR2B" ,"LTR23" ,"LTR80A", "MLT1E3", "LTR103", "MLT1", "THE1C"), 
                   group.by = "archetype",
                   assay = "RNA",
                   Vln.draw.names = T, 
                   feature.pt.size = 2, 
                   ncol = 2, 
                   name = "various_plots/ERV genes")
c("ERVK-6", "LTR2B" ,"LTR23" ,"LTR80A", "MLT1E3", "LTR103", "MLT1", "THE1C") %in% rownames(final.pop.call.integrated.mye.seurat)


## Check expression of markers available for spatial
# Load the markers
spatial.abs <- scan(file = "../Spatial/ab_list.txt", what = "list")

# Clean up the list
spatial.abs <- unique(spatial.abs)

# Update gene symbols
spatial.abs                   <- AnnotationDbi::select(org.Hs.eg.db, keys = spatial.abs, keytype = "ALIAS", columns = "SYMBOL", multiVals = "first")
spatial.abs.symbol            <- spatial.abs[!is.na(spatial.abs$SYMBOL),]
row.names(spatial.abs.symbol) <- NULL

# Write out unmatching entries to check by hand
spatial.abs.handpick <- spatial.abs[is.na(spatial.abs$SYMBOL),"ALIAS"]
write.table(file = "../Spatial/ab_list.no_symbol.txt", spatial.abs.handpick, quote = F, row.names = F, col.names = F)

# Retrieve the curated list and merge it
spatial.abs.handpick <- read.table(file = "../Spatial/ab_list.no_symbol.txt", header = F, col.names = c("ALIAS", "SYMBOL"))
spatial.abs          <- rbind(spatial.abs.symbol, spatial.abs.handpick)

# Fix FAP (wringly associated to CEL, the first hit in the ADbI query)
spatial.abs[spatial.abs$SYMBOL == "CEL", "SYMBOL"] <- "FAP"

# Keep overlap with the significant mac markers
spatial.abs.markers <- final.pop.call.integrated.mye.seurat.markers %>% filter(gene %in% spatial.abs$SYMBOL) %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 1)

# Add back ab aliases
spatial.abs.markers <- merge(spatial.abs.markers, spatial.abs, by.x = "gene", by.y = "SYMBOL")

# Keep only concatenated cluster info per gene
spatial.abs.markers <- spatial.abs.markers %>% group_by(gene) %>% arrange(gene) %>% dplyr::summarise(across(everything(), ~paste0(na.omit(.x), collapse = "; "))) %>% select(all_of(c("gene", "ALIAS", "cluster")))

# Plot the selection
bunchOfCustomPlots(object         = final.pop.call.integrated.mye.seurat, 
                   features       = spatial.abs.markers$gene, 
                   Vln.draw.names = F, 
                   Vln.width      = 25, 
                   Vln.height     = 25, 
                   name           = "../Spatial/Spatial ab selection", 
                   ncol           = 4,
                   Vln.color      = M.int_refined.pop.colors)

bunchOfCustomPlots(object         = integrated.full.seurat, 
                   features       = spatial.abs.markers$gene, 
                   Vln.draw.names = F,
                   Vln.width      = 30, 
                   Vln.height     = 30, 
                   name           = "../Spatial/Spatial ab selection all clusters", 
                   ncol           = 4,
                   Vln.color      = full_set.colors)
customUMAP(integrated.full.seurat, legend.pos = "right", file.name = "../Spatial/UMAP for pop color reference.pdf", cols = full_set.colors, plot.width = 20)

# Save to disk
write.xlsx(spatial.abs.markers, file = "../Spatial/ab_list.pop_markers.xlsx")

