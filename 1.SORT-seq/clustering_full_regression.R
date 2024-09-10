## Clean up, nomralize, and cluster until we are satisified
#==========================================================================================================
##=========================================================================================================
##  BASE: run an unfiltered analysis as base
## Normalize
# run sctransform
all.seur.combined <- SCTransform(all.seur.combined, verbose = T)

# Run PCA and UMAP embedding
all.seur.combined <- RunPCA(all.seur.combined, verbose = T)
all.seur.combined <- RunUMAP(all.seur.combined, dims = 1:30, verbose = T)

# Show some stats
all.seur.combined

## Clustering
set.seed(1)
all.seur.combined <- FindNeighbors(all.seur.combined, dims = 1:30, verbose = T, force.recalc = T)
all.seur.combined <- FindClusters(all.seur.combined, verbose = T, resolution = 1.25, random.seed = 666)
num.clus          <- length(unique(Idents(all.seur.combined)))
DimPlot(all.seur.combined, label = TRUE, pt.size = 2, label.size = 10) + NoLegend() + theme(panel.background = element_blank(),
                                                                                            axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/Unfiltered Round UMAP clusters.pdf", width = 15, height = 15)

DimPlot(all.seur.combined, label = F, pt.size = 2, label.size = 10, group.by = "ID") + theme(panel.background = element_blank(),
                                                                                             axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/Unfiltered Round UMAP ID.pdf", width = 15, height = 15)
DimPlot(all.seur.combined, label = F, pt.size = 2, label.size = 10, group.by = "Sex") + theme(panel.background = element_blank(),
                                                                                              axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/Unfiltered Round UMAP Sex.pdf", width = 15, height = 15)
DimPlot(all.seur.combined, label = F, pt.size = 2, label.size = 10, group.by = "Patient") + theme(panel.background = element_blank(),
                                                                                                  axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/Unfiltered Round UMAP Patient.pdf", width = 15, height = 15)
DimPlot(all.seur.combined, label = F, pt.size = 2, label.size = 10, group.by = "Batch") + theme(panel.background = element_blank(),
                                                                                                axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/Unfiltered Round UMAP Batch.pdf", width = 15, height = 15)

## Scale all genes from RNA assay for diff gene calling and visualisation
# Set default to RNA
DefaultAssay(all.seur.combined) <- "RNA"

# Normalize
all.seur.combined <- NormalizeData(all.seur.combined, normalization.method = "LogNormalize")

#Scale
all.seur.combined <- ScaleData(all.seur.combined, features = row.names(all.seur.combined))


## Check apoptotic markers
VlnPlot(object = all.seur.combined, features = c("KCNQ1OT1","UGDH-AS1", "GHET1", "ROR1-AS1"), ncol = 2) + 
  theme(panel.background = element_blank(),
        text             = element_text(size = 16, face = "bold"),
        axis.text        = element_text(size = 14, face = "plain",angle=45,hjust=1)
  )

ggsave("results/Unfiltered Round Apoptotic markers and lncRNAs.pdf", width = 15, height = 15)

## Define marker genes per cluster
all.seur.combined.markers <- FindAllMarkers(object = all.seur.combined, only.pos = TRUE, min.pct = 0.25, 
                                            thresh = 0.25)

#Save the top9 markers per cluster
sep.markers <- all.seur.combined.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

#plot the top9 markers per cluster
for(i in 0:(num.clus-1)){
  VlnPlot(object = all.seur.combined,
          features = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"]))
  ) + theme(panel.background = element_blank(),
            text             = element_text(size = 16, face = "bold"),
            axis.text        = element_text(size = 14, face = "plain")
  )
  ggsave(paste("results/Unfiltered Round cluster",i, " markers violin.pdf", sep = ""), width = 15, height = 15)
}

VlnPlot(object = all.seur.combined, features = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                                 "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                                 "CD79A", "ACTA2", "CD34", "KIT", "OLR1" ), ncol = 5) + theme(panel.background = element_blank(),
                                                                                                              text             = element_text(size = 16, face = "bold"),
                                                                                                              axis.text        = element_text(size = 14, face = "plain")
                                                 )
ggsave("results/Unfiltered Round cluster defining genes violin.pdf", width = 25, height = 15)

# Store the object
unfiltered.all.seur.combined <- all.seur.combined

#==========================================================================================================
##=========================================================================================================
##  Clustering optimisation: remove cells with mito reads >25%, apop marker reads > 2%, and (long) noncoding (RNA)genes and other junk.
##  Also regress out cell cycle, % mt, % apop genes, and correlation with lncRNAs.
# First, keep only cells with less then 25% mt reads, and less than 2% of apoptotic reads.
# Next, regress out % mt, % apop genes, and correlation to lncRNAs to make a fresh SCT assay
# Then, score cell cycle gene programs.
# Penultimately, regress out cell cycle and all the other variables from the remaining cells.
# Finally: turns out we need to tweak the resolution to split the NK clusters properly
## Subset out the mito and apop cells
all.seur.combined <- subset(all.seur.combined, subset = percent.mt < 25 & percent.KCNQ1OT1 < 2 & percent.UGDH.AS1 < 2 & percent.GHET1 < 2 & percent.ROR1.AS1 < 2)
all.seur.combined
VlnPlot(object = all.seur.combined, features = c("percent.mt"))
ggsave("results/First Round percent.mt.filtered.pdf")

## Identify bad (long noncoding RNA) genes
# first get the ones we can define in bulk from the marker gene lists
sig.markers <- list()
for (theclus in unique(Idents(all.seur.combined))){
  sig.markers[[theclus]]  <- subset(all.seur.combined.markers, p_val_adj < 0.1 & cluster == theclus)$gene
}

kick.list <- list()
for (theclus in unique(Idents(all.seur.combined))){
  kick.list[[theclus]]  <- sig.markers[[theclus]][grep("^RNA|^RPL|^LOC|^LINC|^RPS|-AS1$|^EEF1A|^EEF1B|^FLJ|^MRPS|^MRPL", sig.markers[[theclus]])]
}

# Add some extra known ones by hand
for (theclus in unique(Idents(all.seur.combined))){
  kick.list[[theclus]]  <- c(kick.list[[theclus]], "KCNQ1OT1", "KCNA3", "GHET1", "SH3BGRL3", "NEAT1")
}

# Add a correlation score per cluster
for (theclus in unique(Idents(all.seur.combined))){
  all.seur.combined <- AddModuleScore(all.seur.combined, features = list(bad_genes=kick.list[[theclus]]),  seed = 666, name = paste("bad_lnc_", theclus, sep = ""))
}

# Kick all the bad genes
kick.list <- row.names(all.seur.combined)[grep("^RNA|^RPL|^LOC|^LINC|^RPS|-AS1$|^EEF1A|^EEF1B|^FLJ|^MRPS|^MRPL", row.names(all.seur.combined))]
kick.list <- c(kick.list, c("KCNQ1OT1", "KCNA3", "GHET1", "SH3BGRL3", "NEAT1"))

no.kick.list <- row.names(all.seur.combined)[!row.names(all.seur.combined) %in% kick.list]
all.seur.combined <- subset(all.seur.combined, features = no.kick.list)

# Set up regression variables
vars.to.regress   <- c("percent.mt", "percent.KCNQ1OT1", "percent.UGDH.AS1", "percent.GHET1", "percent.ROR1.AS1")
vars.to.regress   <- c(vars.to.regress, colnames(all.seur.combined@meta.data)[grep("bad", colnames(all.seur.combined@meta.data))])

# run sctransform to re-balance the SCT assay
all.seur.combined <- SCTransform(all.seur.combined, vars.to.regress = vars.to.regress, verbose = T)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. We can segregate this list into markers of G2/M phase and markers of S phase
s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Score the cell cycle programs
all.seur.combined <- CellCycleScoring(all.seur.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay = "SCT")

# Visualize the distribution of cell cycle markers
RidgePlot(all.seur.combined, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
ggsave("results/First Round cell_cylce_Ridgeplot.pdf")

all.seur.combined <- RunPCA(all.seur.combined, features = c(s.genes, g2m.genes), assay = "SCT")
DimPlot(all.seur.combined)
ggsave("results/First Round UMAP_cell_cylce.pdf")

## Normalize, scale, and recluster the subsetted and cell cycle scored object
# Set up regression variables
vars.to.regress   <- c(vars.to.regress, "S.Score", "G2M.Score")

# run sctransform
all.seur.combined <- SCTransform(all.seur.combined, vars.to.regress = vars.to.regress, verbose = T)

# Run PCA and UMAP embedding
all.seur.combined <- RunPCA(all.seur.combined, verbose = T)
all.seur.combined <- RunUMAP(all.seur.combined, dims = 1:30, verbose = T)

# Show some stats
all.seur.combined

## Clustering
set.seed(1)
all.seur.combined <- FindNeighbors(all.seur.combined, dims = 1:30, verbose = T, force.recalc = T)
all.seur.combined <- FindClusters(all.seur.combined, verbose = T, resolution = 1.4, random.seed = 666)
num.clus          <- length(unique(Idents(all.seur.combined)))
DimPlot(all.seur.combined, label = TRUE, pt.size = 2, label.size = 10) + NoLegend() + theme(panel.background = element_blank(),
                                                                                            axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/Filtered Round UMAP clusters.pdf", width = 15, height = 15)

DimPlot(all.seur.combined, label = F, pt.size = 2, label.size = 10, group.by = "ID") + theme(panel.background = element_blank(),
                                                                                             axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/Filtered Round UMAP ID.pdf", width = 15, height = 15)
DimPlot(all.seur.combined, label = F, pt.size = 2, label.size = 10, group.by = "Sex") + theme(panel.background = element_blank(),
                                                                                              axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/Filtered Round UMAP Sex.pdf", width = 15, height = 15)
DimPlot(all.seur.combined, label = F, pt.size = 2, label.size = 10, group.by = "Patient") + theme(panel.background = element_blank(),
                                                                                                  axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/Filtered Round UMAP Patient.pdf", width = 15, height = 15)

DimPlot(all.seur.combined, label = F, pt.size = 2, label.size = 10, group.by = "Batch") + theme(panel.background = element_blank(),
                                                                                                axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/Filtered Round UMAP Batch.pdf", width = 15, height = 15)

# Quick check for the NK clusters
VlnPlot(all.seur.combined, features = c("CD3E", "CD4", "CD8A", "CD14", "CD68", "NCAM1"))
ggsave("results/Filtered Round NK check SCT.pdf", width = 15, height = 15)

## Scale all genes from RNA assay for diff gene calling and visualisation
# Set default to RNA
DefaultAssay(all.seur.combined) <- "RNA"

# Normalize
all.seur.combined <- NormalizeData(all.seur.combined, normalization.method = "LogNormalize")

#Scale
all.seur.combined <- ScaleData(all.seur.combined, features = row.names(all.seur.combined))


## Check apoptotic markers
VlnPlot(object = all.seur.combined, features = c("KCNQ1OT1","UGDH-AS1", "GHET1", "ROR1-AS1"), ncol = 2) + 
  theme(panel.background = element_blank(),
        text             = element_text(size = 16, face = "bold"),
        axis.text        = element_text(size = 14, face = "plain",angle=45,hjust=1)
  )

ggsave("results/Filtered Round Apoptotic markers and lncRNAs.pdf", width = 15, height = 15)

## Define marker genes per cluster
all.seur.combined.markers <- FindAllMarkers(object = all.seur.combined, only.pos = TRUE, min.pct = 0.25, 
                                            thresh = 0.25)

#Save the top9 markers per cluster
sep.markers <- all.seur.combined.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

#plot the top9 markers per cluster
for(i in 0:(num.clus-1)){
  VlnPlot(object = all.seur.combined,
          features = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"]))
  ) + theme(panel.background = element_blank(),
            text             = element_text(size = 16, face = "bold"),
            axis.text        = element_text(size = 14, face = "plain")
  )
  ggsave(paste("results/Filtered Round cluster",i, " markers violin.pdf", sep = ""), width = 15, height = 15)
}

VlnPlot(object = all.seur.combined, features = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                                 "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                                 "CD79A", "ACTA2", "CD34", "KIT", "OLR1" ), ncol = 5) + theme(panel.background = element_blank(),
                                                                                                              text             = element_text(size = 16, face = "bold"),
                                                                                                              axis.text        = element_text(size = 14, face = "plain")
                                                 )
ggsave("results/Filtered Round cluster defining genes violin.pdf", width = 25, height = 15)

# Quick check for the NK clusters
VlnPlot(all.seur.combined, features = c("CD3E", "CD4", "CD8A", "CD14", "CD68", "NCAM1"))
ggsave("results/Filtered Round NK check.pdf", width = 15, height = 15)

# Store the object
round.one.all.seur.combined <- all.seur.combined


##==========================================================================================================
##==========================================================================================================
## Inspect the final clusters
#Get cluster colors for subsequent plots
t                      <- DimPlot(all.seur.combined)
tbuild                 <- ggplot_build(t)
tdata                  <- tbuild$data[[1]] 
tdata                  <- tdata[order(tdata$group), ]
cluster_colours        <- unique(tdata$colour)
names(cluster_colours) <- 0:(length(levels(Idents(all.seur.combined)))-1)
# Run this next line only afer establishing cluster identities
#names(cluster_colours) <- new.ident

#plot the top9 markers per cluster
for(i in 0:(num.clus-1)){
  FeaturePlot(object = all.seur.combined,
              features = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])),
              cols = c("grey", "blue"), 
              reduction = "umap", pt.size = 1
  ) + theme(panel.background = element_blank(),
            text             = element_text(size = 16, face = "bold"),
            axis.text        = element_text(size = 14, face = "plain")
  )
  ggsave(paste("results/cluster",i, " markers.pdf", sep = ""), width = 15, height = 15)
  VlnPlot(object = all.seur.combined,
          features = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"]))
  ) + theme(panel.background = element_blank(),
            text             = element_text(size = 16, face = "bold"),
            axis.text        = element_text(size = 14, face = "plain")
  )
  ggsave(paste("results/cluster",i, " markers violin.pdf", sep = ""), width = 15, height = 15)
}

#Look at some genes
FeaturePlot(object = all.seur.combined, pt.size = 1, features = c("CD3E", "CD4", "CD8A", "CD14", 
                                                                  "CD68", "CD79A", "MYH11", "CD34", "KIT"), min.cutoff = "q9", cols = c("lightgrey", 
                                                                                                                                        "blue")) + theme(panel.background = element_blank(),
                                                                                                                                                         text             = element_text(size = 16, face = "bold"),
                                                                                                                                                         axis.text        = element_text(size = 14, face = "plain")
                                                                                                                                        )
ggsave("results/cluster defining genes.pdf", width = 15, height = 15)

VlnPlot(object = all.seur.combined, features = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                                 "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                                 "CD79A", "ACTA2", "CD34", "KIT", "IL7R" ), ncol = 5) + theme(panel.background = element_blank(),
                                                                                                              text             = element_text(size = 16, face = "bold"),
                                                                                                              axis.text        = element_text(size = 14, face = "plain")
                                                 )
ggsave("results/cluster defining genes extended violin.pdf", width = 25, height = 15)

VlnPlot(object = all.seur.combined, features = c("CD3E", "CD4", "CD8A", "CD14", 
                                                 "CD68", "CD79A", "ACTA2", "CD34", "KIT")) + theme(panel.background = element_blank(),
                                                                                                   text             = element_text(size = 16, face = "bold"),
                                                                                                   axis.text        = element_text(size = 14, face = "plain")
                                                 )
ggsave("results/cluster defining genes violin.pdf", width = 15, height = 15)

#Plot top markers in a heatmap
top10 <- all.seur.combined.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(object = all.seur.combined, features = top10$gene) + theme(panel.background = element_blank(),
                                                                     text             = element_text(size = 16, face = "bold"),
                                                                     axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/clusters.pdf", width = 15, height =15)


##==========================================================================================================
##==========================================================================================================
## Determine cluster cell types by expression of marker genes

#Rename clusters and replot UMAP
new.ident <-  c(             "CD3+ T Cells I",                    "CD3+ T Cells II",                      "CD3+ T Cells III",                               #0,  1,  2
                             "CD3+ T Cells IV",                   "CD3+ T Cells V",                       "CD3+ T Cells VI",                                #3,  4,  5
                             "CD3+ T Cells VII",                  "ACTA2+ Smooth Muscle Cells",           "CD68+IL18+TLR4+TREM2+ Resident macrophages",     #6,  7,  8
                             "CD34+ Endothelial Cells I",         "CD3+CD56+ NK Cells I",                 "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", #9,  10, 11
                             "CD68+ABCA1+OLR1+TREM2+ Foam Cells", "CD34+ Endothelial Cells II",           "CD68+CD1C+ Dendritic Cells",                     #12, 13, 14
                             "CD3+CD56+ NK Cells II",             "CD79A+ Class-switched Memory B Cells", "FOXP3+ T Cells",                                 #15, 16, 17
                             "CD68+KIT+ Mast Cells" ,             "CD68+CD4+ Monocytes",                  "CD79+ Plasma B Cells")                           #18, 19, 20
length(new.ident)
length(unique(new.ident))
names(new.ident)  <- 0:20
new.ident
safetysafe.all.seur.combined <- all.seur.combined
all.seur.combined <- RenameIdents(all.seur.combined, new.ident)
DimPlot(all.seur.combined, label = T, pt.size = 2, label.size = 10) + theme(legend.position  = "none",
                                                                            panel.background = element_blank(),
                                                                            text             = element_text(size = 14, face = "bold"),
                                                                            axis.text        = element_text(size = 14, face = "plain")
)
ggsave("results/UMAP clusters labels.pdf", width = 15, height = 15)

# Save the object
saveRDS(all.seur.combined, "Seurat Objects/main.seurat.filtered.clustered.RDS")
