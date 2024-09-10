#=========================================================================================================================
## Map the refined plaque and PBMC clusters back to the main populations
##========================================================================================================================
dir.create("myeloid_backmapping_results", showWarnings = F)
## First map back to the integrated myeloid cells

# Get refined idents per cell
mye.idents        <- c(as.vector(Idents(pbmc.monocytes.seurat)),
                       as.vector(Idents(plaque.macrophages.seurat))
)
names(mye.idents) <- c(names(Idents(pbmc.monocytes.seurat)),
                       names(Idents(plaque.macrophages.seurat))
)

# Get original idents per cell
old.mye.ident         <- as.vector(Idents(mye.patients.seurat))
names(old.mye.ident)  <- names(Idents(mye.patients.seurat))

# Discard cells that were discarded in the refined populations
mye.patients.seurat <- subset(mye.patients.seurat, cells = names(mye.idents))

# Rerun the clustering to get the cells together properly
DefaultAssay(mye.patients.seurat) <- "integrated"
mye.patients.seurat <- RunPCA(mye.patients.seurat, verbose = T, seed.use = 666)
mye.patients.seurat <- RunUMAP(mye.patients.seurat, dims = 1:30, verbose = T)

# Set the new (refined) idents
Idents(mye.patients.seurat) <- mye.idents


## Save the object
saveRDS(mye.patients.seurat, "Seurat_Objects/myeloid.cells.refined.seurat.RDS")
#mye.patients.seurat <- readRDS(file = "Seurat_Objects/myeloid.cells.refined.seurat.RDS")


## Define useful pop colors
M_refined.pops <- levels(Idents(mye.patients.seurat))

# Monocytes
M_refined.pop.colors        <- colorSpacer(startcolor = "cyan3", middlecolors = "dodgerblue", endcolor = "blue4", steps = 3, return.colors = T)
names(M_refined.pop.colors) <- M_refined.pops[grep(" Monocytes", M_refined.pops)][c(1,3,2)]

# Inflammatory Mo-Macs
tmp.names                   <- names(M_refined.pop.colors)
M_refined.pop.colors        <- c(M_refined.pop.colors, colorSpacer(startcolor = "firebrick4", endcolor = "goldenrod3", steps = 4, return.colors = T))
names(M_refined.pop.colors) <- c(tmp.names, M_refined.pops[grep("derived Macro", M_refined.pops)])

# Foam Cells
tmp.names                   <- names(M_refined.pop.colors)
M_refined.pop.colors        <- c(M_refined.pop.colors, colorSpacer(startcolor = "plum1", endcolor = "darkslateblue", steps = 2, return.colors = T))
names(M_refined.pop.colors) <- c(tmp.names, M_refined.pops[grep("Foamy Mac", M_refined.pops)])

# Res Macs
tmp.names                   <- names(M_refined.pop.colors)
M_refined.pop.colors        <- c(M_refined.pop.colors, colorSpacer(startcolor = "chartreuse1", endcolor = "chartreuse4", steps = 3, return.colors = T))
names(M_refined.pop.colors) <- c(tmp.names, M_refined.pops[grep("Resident", M_refined.pops)])

length(M_refined.pop.colors)

## And plot
# First plot all
customUMAP(object     = mye.patients.seurat,
           label      = F, 
           pt.size    = 3,
           label.size = 10, 
           reduction  = "umap", 
           shuffle    = T,
           cols       = M_refined.pop.colors, file.name = "myeloid_backmapping_results/Myeloid UMAP clusters labels.pdf",
           plot.width = 25, 
           plot.height = 15
          )

# Highlight each cluster
dir.create("myeloid_backmapping_results/Myeloid cells cluster highlights", showWarnings = F, recursive = T)
for(theIdent in levels(Idents(mye.patients.seurat))){
  customUMAP(object     = mye.patients.seurat,
             label      = F, 
             pt.size    = 3,
             label.size = 10, 
             reduction  = "umap", 
             shuffle    = T,
             cols       = M_refined.pop.colors, file.name = paste("myeloid_backmapping_results/Myeloid cells cluster highlights/", theIdent, " UMAP.pdf", sep = ""),
             plot.width = 20, 
             plot.height = 15,
             sizes.highlight = 4,
             cells.highlight = names(mye.idents)[which(mye.idents == theIdent)])
}
 
#=========================================================================================================================
## Now map back to the full complement of cells
##========================================================================================================================
# Get original idents per cell
old.ident         <- as.vector(Idents(patients.seurat)) 
names(old.ident)  <- names(Idents(patients.seurat))

# Discard (myeloid!) cells that were discarded in the refined populations
not.mye.cells           <- row.names(patients.seurat@meta.data[grep("^Mono$", invert = T, patients.seurat$predicted.celltype.l1),])
refined.patients.seurat <- subset(patients.seurat, cells = c(not.mye.cells, names(mye.idents)))

# Set the new idents
Idents(refined.patients.seurat, cells = names(mye.idents)) <- mye.idents

## Save the object
saveRDS(refined.patients.seurat, "Seurat_Objects/main.patient.integrated.mapped.refined_populations.seurat.RDS")
#refined.patients.seurat <- readRDS(file = "Seurat_Objects/main.patient.integrated.mapped.refined_populations.seurat.RDS")


## And plot
# First plot all
customUMAP(object     = refined.patients.seurat,
           label      = F, 
           pt.size    = 3,
           label.size = 10, 
           reduction  = "umap", 
           shuffle    = T,
           cols       = M_refined.pop.colors, file.name = "myeloid_backmapping_results/Full set UMAP clusters labels.pdf",
           plot.width = 25, 
           plot.height = 15
)

customUMAP(object     = refined.patients.seurat,
           label      = F, 
           pt.size    = 3,
           label.size = 10, 
           reduction  = "umap", legend.pos = "right",
           shuffle    = T,
           cols       = M_refined.pop.colors, file.name = "myeloid_backmapping_results/Full set UMAP clusters labels legend.pdf",
           plot.width = 25, 
           plot.height = 15
)

# Highlight each cluster
dir.create("myeloid_backmapping_results/Full set cluster highlights", showWarnings = F, recursive = T)
for(theIdent in levels(Idents(refined.patients.seurat))){
  customUMAP(object     = refined.patients.seurat,
             label      = F, 
             pt.size    = 3,
             label.size = 10, 
             reduction  = "umap", 
             shuffle    = T,
             cols       = M_refined.pop.colors, file.name = paste("myeloid_backmapping_results/Full set cluster highlights/", theIdent, " UMAP.pdf", sep = ""),
             plot.width = 20, 
             plot.height = 15,
             sizes.highlight = 4,
             cells.highlight = names(mye.idents)[which(mye.idents == theIdent)])
}



#=========================================================================================================================
## Explore the link between moncoytes and plaque macs
##========================================================================================================================
# Plot som emono and mac genes
bunchOfCustomPlots(features  = c("CD14","FCGR3A", "CD68", 'CD86', "SELL", "CSF1R", "ITGAM", "ITGAX", "PECAM1", "CCR2", "PADI4", "SECISBP2L"), 
                   object    = refined.patients.seurat,
                   ncol      = 3, Vln.stack = T,
                   name      = "myeloid_backmapping_results/macro and mono genes",
                   Vln.width = 15, Vln.height = 15)




##==================================================================================
##==================================================================================
## Use new 10X - CEL-seq integrated populations
# Get new refined idents
int.idents <- Idents(integrated.mye.seurat)

# Add monocytes
int.idents <- c(int.idents, Idents(mye.patients.seurat)[grep("Mono", Idents(mye.patients.seurat))])

mye.patients.seurat <- AddMetaData(mye.patients.seurat, metadata = int.idents[names(int.idents) %in% colnames(mye.patients.seurat)], col.name = "method.int.idents")

## Save the object
saveRDS(mye.patients.seurat, "Seurat_Objects/myeloid.cells.refined.seurat.RDS")
#mye.patients.seurat <- readRDS(file = "Seurat_Objects/myeloid.cells.refined.seurat.RDS")


# Plot
customUMAP(object     = mye.patients.seurat, 
           group.by   = "method.int.idents", 
           cols       = M.int_refined.pop.colors, 
           title      = "43p refined populations", 
           legend.pos = "right", 
           plot.width = 20, shuffle = F,
           file.name  = "myeloid_backmapping_results/UMAP 43p integrated populations.pdf" )


## Also check the pre-merged integreated clusters 2, 4, and 6
int.idents <- integrated.mye.seurat$method.int.idents.pre.merge
mye.patients.seurat <- AddMetaData(mye.patients.seurat, metadata = int.idents[names(int.idents) %in% colnames(mye.patients.seurat)], col.name = "numbered.method.int.idents")

# Plot
customUMAP(object     = mye.patients.seurat, 
           group.by   = "numbered.method.int.idents", 
           title      = "43p numbered populations", 
           legend.pos = "right", 
           plot.width = 15,
           file.name  = "myeloid_backmapping_results/UMAP 43p integrated numbered populations.pdf" )




