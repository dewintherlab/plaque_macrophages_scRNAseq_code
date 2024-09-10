#======================================================================
## Analyze GM_CSF signatures in th emac populations
#======================================================================
dir.create("csfmac", showWarnings = F)

## Construct signatures using GSE135491 data of human monocyte derived macrophages.
## Load and preprocess data
# Read RPKM file
csf.macs.RPKM <- read.table(file = "raw_data/GSE135491_Processeddatafile_MacrophagesProject_GM_vs_M_CSF.clean.txt", header = T, row.names = 1, sep = "\t")

# Remove (for us irrelevant) TAM samples
csf.macs.RPKM <- csf.macs.RPKM[,1:10]

# Remove (for us irrelevant) Ribosomal genes and small or nc RNAs
csf.macs.RPKM <- csf.macs.RPKM[grep("RPL|SNOR|RPS|LINC|RNA|ERV|RNU|RACK|EEF|EIF|RMR|RP5|RP11|RN7|RP3|AC00|MIR|RPP|GAS6|DLEU|CTA-|CTB-|AP0|RNF|orf|XXbac|CTD|RP4|Metazoa|KCNQ|GS1|XXyac|AC127|AC10|LA16|CTC-|KB-|AE00|AC0",csf.macs.RPKM[,"gene_name"], invert = T),]

# Split off gene info for later
csf.macs.genes <- csf.macs.RPKM[,1,drop = F]
csf.macs.RPKM  <- csf.macs.RPKM[,-1]

# Define Group indices
MCSF  <- c(1, 2, 3, 7, 8, 9)
GMCSF <- c(4, 5, 6)

# Remove lowly expressed genes (keep median RPKM > 1 in at least one group)
keep           <- apply(csf.macs.RPKM[,MCSF],1,median) > 1 | apply(csf.macs.RPKM[,GMCSF],1,median) > 1
csf.macs.RPKM  <- csf.macs.RPKM[ keep,]
csf.macs.genes <- csf.macs.genes[keep,,drop = F]

# Remove NAs
indx           <- apply(csf.macs.RPKM, 1, function(x) any(is.na(x)))
csf.macs.RPKM  <- csf.macs.RPKM[ -!indx,]
csf.macs.genes <- csf.macs.genes[-!indx,,drop = F]

# Construct metadata table
colData      <- data.frame(row.names = colnames(csf.macs.RPKM), Type = c("M-CSF", "M-CSF", "M-CSF", "GM-CSF", "GM-CSF", "GM-CSF", "M-CSF", "M-CSF", "M-CSF"))
colData$Type <- factor(colData$Type, levels = c("M-CSF", "GM-CSF"))

# Define annotations
annColors <- list("Type" = c("M-CSF" = "dodgerblue", "GM-CSF" = "goldenrod"))

## Inspect the data
# Extract top variable genes
topVarGenes <- head(order(-rowVars(as.matrix(csf.macs.RPKM))),5000)

# Heatmap of top variable genes
pheatmap(mat = csf.macs.RPKM[topVarGenes,], scale = "row", cluster_cols = T, cluster_rows = T, annotation_col = colData, annotation_colors = annColors, show_colnames = T, show_rownames = F)

# Looks like CA2, CB2, and CB3 are outliers, let's remove them
csf.macs.RPKM <- csf.macs.RPKM[,c(1,3,4,5,6,7)]
colData       <- colData[c(1,3,4,5,6,7),, drop = F]

# Extract top variable genes
topVarGenes <- head(order(-rowVars(as.matrix(csf.macs.RPKM))),1000)

# Heatmap of top variable genes
pheatmap(mat = csf.macs.RPKM[topVarGenes,], scale = "row", cluster_cols = T, cluster_rows = T, annotation_col = colData, annotation_colors = annColors, show_colnames = T, show_rownames = F, filename = "csfmac/topvargenes.pdf")

# Define Group indices
MCSF  <- c(1, 2, 6)
GMCSF <- c(3, 4, 5)

# Remove lowly expressed genes (keep median RPKM > 1 in at least one group)
keep           <- apply(csf.macs.RPKM[,MCSF],1,median) > 1 | apply(csf.macs.RPKM[,GMCSF],1,median) > 1
csf.macs.RPKM  <- csf.macs.RPKM[ keep,]
csf.macs.genes <- csf.macs.genes[keep,,drop = F]

# Remove NAs
indx           <- apply(csf.macs.RPKM, 1, function(x) any(is.na(x)))
csf.macs.RPKM  <- csf.macs.RPKM[ -!indx,]
csf.macs.genes <- csf.macs.genes[-!indx,,drop = F]

# Kmeans heatmap of all genes
pheatmap(mat = csf.macs.RPKM, scale = "row", cluster_cols = T, cluster_rows = T, annotation_col = colData, annotation_colors = annColors, show_colnames = F, show_rownames = F, kmeans_k = 2)

## Conclusion: shit data or nothing much changes?
## Let's manually extract the largest differences
# Get type averages
csf.macs.meds <- data.frame(row.names = row.names(csf.macs.RPKM), "MCSF" = apply(csf.macs.RPKM[,MCSF],1,mean), "GMCSF" = apply(csf.macs.RPKM[,GMCSF],1,mean))

# Get the order of expression for both types
mcsf.high  <- order(csf.macs.meds$MCSF, decreasing = T)
gmcsf.high <- order(csf.macs.meds$GMCSF, decreasing = T)

# Plot them
pheatmap(csf.macs.meds[mcsf.high[ 1:100],],cluster_rows = F, cluster_cols = F, show_rownames = F)
pheatmap(csf.macs.meds[gmcsf.high[1:100],],cluster_rows = F, cluster_cols = F, show_rownames = F)

csf.macs.genes[mcsf.high[1:10],]
csf.macs.genes[gmcsf.high[1:10],]

# Get the diff
csf.macs.meds$log2FC <- log2(csf.macs.meds$MCSF / csf.macs.meds$GMCSF)

# Get the order of expression for both types
mcsf.high  <- order(csf.macs.meds$log2FC, decreasing = T)
gmcsf.high <- order(csf.macs.meds$log2FC, decreasing = F)

# Plot them
pheatmap(csf.macs.meds[mcsf.high[ 1:50],1:2],cluster_rows = T, cluster_cols = F, show_rownames = F)
pheatmap(csf.macs.meds[gmcsf.high[1:50],1:2],cluster_rows = T, cluster_cols = F, show_rownames = F)

## Zoom in to the top |log2FC|> 2
# Fetch the top 
mcsf.high.macs  <- csf.macs.meds[mcsf.high[ 1:sum(csf.macs.meds$log2FC > 2)],]
gmcsf.high.macs <- csf.macs.meds[gmcsf.high[1:sum(csf.macs.meds$log2FC < -2)],]

# Order the top
mcsf.high.macs  <- mcsf.high.macs[ order( mcsf.high.macs$MCSF,  decreasing = T), ]
gmcsf.high.macs <- gmcsf.high.macs[order(gmcsf.high.macs$GMCSF, decreasing = T), ]

# Keep only the expressed top (RPKM > 1)
mcsf.high.macs  <- mcsf.high.macs[  mcsf.high.macs$MCSF  > 1, ]
gmcsf.high.macs <- gmcsf.high.macs[gmcsf.high.macs$GMCSF > 1, ]

# Retrieve the gene names
mcsf.high.genes <- csf.macs.genes[row.names(csf.macs.genes) %in% row.names(mcsf.high.macs),, drop = F]
mcsf.high.genes <- mcsf.high.genes[row.names(mcsf.high.macs),]

gmcsf.high.genes <- csf.macs.genes[row.names(csf.macs.genes) %in% row.names(gmcsf.high.macs),, drop = F]
gmcsf.high.genes <- gmcsf.high.genes[row.names(gmcsf.high.macs),]


# Plot the top
pheatmap(mcsf.high.macs[,1:2], color = colorSpacer(endcolor = "dodgerblue", startcolor = "grey", steps = length(mcsf.high.genes), return.colors = T), cellwidth = 30, cellheight = 10, border_color = NA, cluster_rows = F, cluster_cols = F, show_rownames = T, angle_col = 45, labels_row = mcsf.high.genes, filename = "csfmac/mcsftopFC.rpkm.pdf")
pheatmap(gmcsf.high.macs[,1:2], color = colorSpacer(endcolor = "goldenrod", startcolor = "grey", steps = length(gmcsf.high.genes), return.colors = T), cellwidth = 30, cellheight = 10, border_color = NA, cluster_rows = F, cluster_cols = F, show_rownames = T, angle_col = 45, labels_row = gmcsf.high.genes, filename = "csfmac/gmcsftopFC.rpkm.pdf")

# Fix the FC
mcsf.high.macs$log2FC[!is.finite( mcsf.high.macs$log2FC)] <- max( mcsf.high.macs$log2FC[is.finite( mcsf.high.macs$log2FC)]) + 1

gmcsf.high.macs$log2FC[!is.finite(gmcsf.high.macs$log2FC)] <- min(gmcsf.high.macs$log2FC[is.finite(gmcsf.high.macs$log2FC)]) - 1

# Invert FC for ease
gmcsf.high.macs$log2FC <- abs(gmcsf.high.macs$log2FC)

# Plot the FC
pheatmap(mcsf.high.macs$log2FC, 
         color = colorSpacer(endcolor = "black",middlecolors = "red", startcolor = "grey", steps = length(mcsf.high.genes), return.colors = T),
         cellwidth = 30, cellheight = 10, border_color = NA, cluster_rows = F, cluster_cols = F, show_rownames = T, angle_col = 45, 
         labels_row = mcsf.high.genes,
         legend_breaks = c(min(mcsf.high.macs$log2FC[is.finite(mcsf.high.macs$log2FC)]), 4, 6, 8, 10, 12, 14, max(mcsf.high.macs$log2FC[is.finite(mcsf.high.macs$log2FC)])),
         legend_labels = c(2,4,6,8,10,12,14, "Inf"), filename = "csfmac/mcsftopFC.log2FC.pdf")

pheatmap(gmcsf.high.macs$log2FC, 
         color = colorSpacer(endcolor = "black",middlecolors = "red", startcolor = "grey", steps = length(gmcsf.high.genes), return.colors = T),
         cellwidth = 30, cellheight = 10, border_color = NA, cluster_rows = F, cluster_cols = F, show_rownames = T, angle_col = 45, 
         labels_row = gmcsf.high.genes,
         legend_breaks = c(min(gmcsf.high.macs$log2FC[is.finite(gmcsf.high.macs$log2FC)]), 4, 6, 8, 10, 12, 14, max(gmcsf.high.macs$log2FC[is.finite(gmcsf.high.macs$log2FC)])),
         legend_labels = c(2,4,6,8,10,12,14, "Inf"), filename = "csfmac/gmcsftopFC.log2FC.pdf")


##==================================================
## Check for enrichment in the single cell data
# Make genesets for use with AUCell
# Fetch the top 
mcsf.high.macs  <- csf.macs.meds[mcsf.high[ 1:sum(csf.macs.meds$log2FC > 2)],]
gmcsf.high.macs <- csf.macs.meds[gmcsf.high[1:sum(csf.macs.meds$log2FC < -2)],]

# Order the top
mcsf.high.macs  <- mcsf.high.macs[ order( mcsf.high.macs$MCSF,  decreasing = T), ]
gmcsf.high.macs <- gmcsf.high.macs[order(gmcsf.high.macs$GMCSF, decreasing = T), ]

# Retrieve the gene names
mcsf.high.genes <- csf.macs.genes[row.names(csf.macs.genes) %in% row.names(mcsf.high.macs),, drop = F]
mcsf.high.genes <- mcsf.high.genes[row.names(mcsf.high.macs),]

gmcsf.high.genes <- csf.macs.genes[row.names(csf.macs.genes) %in% row.names(gmcsf.high.macs),, drop = F]
gmcsf.high.genes <- gmcsf.high.genes[row.names(gmcsf.high.macs),]

geneSets <- list("MCSF" = mcsf.high.genes, "GMCSF" = gmcsf.high.genes)

# Retrieve the expression matrix from the seurat object
exprMat <- GetAssayData(final.pop.call.integrated.mye.velocyto.seurat, "counts")
DimPlot(final.pop.call.integrated.mye.velocyto.seurat) + theme(legend.position = "none")
mcsf.high.genes %in% row.names(exprMat)

# Calculate enrichment scores
cells_AUC <- AUCell_run(exprMat, geneSets, BiocParallel::MulticoreParam(5))

# Optimize tresholds
set.seed(333)
par(mfrow=c(2,1)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
warningMsg[which(warningMsg!="")]

cells_assignment$MCSF$aucThr$selected
cells_assignment$GMCSF$aucThr$selected
length(cells_assignment$MCSF$assignment)
length(cells_assignment$GMCSF$assignment)


##==================================================
## Try a microarray set directly form GEO
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE99056", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13497", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "111XXX000XXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("MCSF","GMCSF"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","SPOT_ID","GB_ACC","SEQUENCE"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE99056", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE99056", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 3, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=3", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE99056")

## Extract sig genes
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.01, lfc=2)

mcsf.up  <- row.names(dT@.Data)[dT@.Data == 1]
gmcsf.up <- row.names(dT@.Data)[dT@.Data == -1]

mcsf.up.genes  <- select(x = hgug4112a.db, keys = mcsf.up,  columns = "SYMBOL", keytype = "PROBEID")$SYMBOL
gmcsf.up.genes <- select(x = hgug4112a.db, keys = gmcsf.up, columns = "SYMBOL", keytype = "PROBEID")$SYMBOL

mcsf.up.genes  <- mcsf.up.genes[ !is.na(mcsf.up.genes)]
gmcsf.up.genes <- gmcsf.up.genes[!is.na(gmcsf.up.genes)]

##==================================================
## Check for enrichment in the single cell data
# Make genesets for use with AUCell
geneSets <- list("MCSF" = mcsf.up.genes, "GMCSF" = gmcsf.up.genes)


# Retrieve the expression matrix from the seurat object
exprMat <- GetAssayData(final.pop.call.integrated.mye.velocyto.seurat, assay = "RNA", layer = "counts")
DimPlot(final.pop.call.integrated.mye.velocyto.seurat) + theme(legend.position = "none")
mcsf.up.genes %in% row.names(exprMat)

# Calculate enrichment scores
cells_AUC <- AUCell_run(exprMat, geneSets, BiocParallel::MulticoreParam(5), aucMaxRank = 20)

# Optimize tresholds
set.seed(333)
par(mfrow=c(2,1)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
warningMsg[which(warningMsg!="")]

cells_assignment$MCSF$aucThr$thresholds
cells_assignment$MCSF$aucThr$selected

cells_assignment$GMCSF$aucThr$thresholds
cells_assignment$GMCSF$aucThr$selected
length(cells_assignment$MCSF$assignment)
length(cells_assignment$GMCSF$assignment)

# Add to seurat metadata
act.sig <- data.frame(row.names = colnames(final.pop.call.integrated.mye.velocyto.seurat), "act.sig" = rep("none", length(colnames(final.pop.call.integrated.mye.velocyto.seurat))))
act.sig[cells_assignment$MCSF$assignment,  "act.sig"] <- "M-CSF"
act.sig[cells_assignment$GMCSF$assignment, "act.sig"] <- "GM-CSF"

archetype.integrated.mye.velocyto.seurat <- AddMetaData(archetype.integrated.mye.velocyto.seurat, metadata = act.sig, col.name = "act.sig")
DimPlot(archetype.integrated.mye.velocyto.seurat, group.by = "act.sig")
customUMAP(object = archetype.integrated.mye.velocyto.seurat, group.by = "act.sig",   pt.size = 3, shuffle = T, legend.pos = "top", seed = 666, file.name = "csfmac/umap.pdf", cols = c("M-CSF" = "dodgerblue", "GM-CSF" = "orangered", "none" = "grey"))
customUMAP(object = archetype.integrated.mye.velocyto.seurat, pt.size = 3, shuffle = T, legend.pos = "top", seed = 666, file.name = "csfmac/umap archatypes.pdf", cols = archetype.colors)


## Construct signatures using GSE135491 data of human monocyte derived macrophages after LPS stimulation.
## Load and preprocess data
# Read RPKM file
lps.macs.counts <- read.xlsx2(file = "raw_data/GSE227737_counts.xlsx", sheetIndex = 1)

head(lps.macs.counts)

# Remove (for us irrelevant) Ribosomal genes and small or nc RNAs
lps.macs.counts <- lps.macs.counts[grep("RPL|SNOR|RPS|LINC|RNA|ERV|RNU|RACK|EEF|EIF|RMR|RP5|RP11|RN7|RP3|AC00|MIR|RPP|GAS6|DLEU|CTA-|CTB-|AP0|RNF|orf|XXbac|CTD|RP4|Metazoa|KCNQ|GS1|XXyac|AC127|AC10|LA16|CTC-|KB-|AE00|AC0",lps.macs.counts[,"gene_name"], invert = T),]

# Split off gene info for later
lps.macs.genes <- lps.macs.counts[,1,drop = F]
lps.macs.counts  <- lps.macs.counts[,-1]

# Define Group indices
Mlps  <- c(1, 2, 3, 7, 8, 9)
GMlps <- c(4, 5, 6)

# Remove lowly expressed genes (keep median RPKM > 1 in at least one group)
keep           <- apply(lps.macs.counts[,Mlps],1,median) > 1 | apply(lps.macs.counts[,GMlps],1,median) > 1
lps.macs.counts  <- lps.macs.counts[ keep,]
lps.macs.genes <- lps.macs.genes[keep,,drop = F]

# Remove NAs
indx           <- apply(lps.macs.counts, 1, function(x) any(is.na(x)))
lps.macs.counts  <- lps.macs.counts[ -!indx,]
lps.macs.genes <- lps.macs.genes[-!indx,,drop = F]

# Construct metadata table
colData      <- data.frame(row.names = colnames(lps.macs.counts), Type = c("M-lps", "M-lps", "M-lps", "GM-lps", "GM-lps", "GM-lps", "M-lps", "M-lps", "M-lps"))
colData$Type <- factor(colData$Type, levels = c("M-lps", "GM-lps"))

# Define annotations
annColors <- list("Type" = c("M-lps" = "dodgerblue", "GM-lps" = "goldenrod"))

## Inspect the data
# Extract top variable genes
topVarGenes <- head(order(-rowVars(as.matrix(lps.macs.counts))),5000)

# Heatmap of top variable genes
pheatmap(mat = lps.macs.counts[topVarGenes,], scale = "row", cluster_cols = T, cluster_rows = T, annotation_col = colData, annotation_colors = annColors, show_colnames = T, show_rownames = F)

# Looks like CA2, CB2, and CB3 are outliers, let's remove them
lps.macs.counts <- lps.macs.counts[,c(1,3,4,5,6,7)]
colData       <- colData[c(1,3,4,5,6,7),, drop = F]

# Extract top variable genes
topVarGenes <- head(order(-rowVars(as.matrix(lps.macs.counts))),1000)

# Heatmap of top variable genes
pheatmap(mat = lps.macs.counts[topVarGenes,], scale = "row", cluster_cols = T, cluster_rows = T, annotation_col = colData, annotation_colors = annColors, show_colnames = T, show_rownames = F, filename = "lpsmac/topvargenes.pdf")

# Define Group indices
Mlps  <- c(1, 2, 6)
GMlps <- c(3, 4, 5)

# Remove lowly expressed genes (keep median RPKM > 1 in at least one group)
keep           <- apply(lps.macs.counts[,Mlps],1,median) > 1 | apply(lps.macs.counts[,GMlps],1,median) > 1
lps.macs.counts  <- lps.macs.counts[ keep,]
lps.macs.genes <- lps.macs.genes[keep,,drop = F]

# Remove NAs
indx           <- apply(lps.macs.counts, 1, function(x) any(is.na(x)))
lps.macs.counts  <- lps.macs.counts[ -!indx,]
lps.macs.genes <- lps.macs.genes[-!indx,,drop = F]

# Kmeans heatmap of all genes
pheatmap(mat = lps.macs.counts, scale = "row", cluster_cols = T, cluster_rows = T, annotation_col = colData, annotation_colors = annColors, show_colnames = F, show_rownames = F, kmeans_k = 2)

## Conclusion: shit data or nothing much changes?
## Let's manually extract the largest differences
# Get type averages
lps.macs.meds <- data.frame(row.names = row.names(lps.macs.counts), "Mlps" = apply(lps.macs.counts[,Mlps],1,mean), "GMlps" = apply(lps.macs.counts[,GMlps],1,mean))

# Get the order of expression for both types
mlps.high  <- order(lps.macs.meds$Mlps, decreasing = T)
gmlps.high <- order(lps.macs.meds$GMlps, decreasing = T)

# Plot them
pheatmap(lps.macs.meds[mlps.high[ 1:100],],cluster_rows = F, cluster_cols = F, show_rownames = F)
pheatmap(lps.macs.meds[gmlps.high[1:100],],cluster_rows = F, cluster_cols = F, show_rownames = F)

lps.macs.genes[mlps.high[1:10],]
lps.macs.genes[gmlps.high[1:10],]

# Get the diff
lps.macs.meds$log2FC <- log2(lps.macs.meds$Mlps / lps.macs.meds$GMlps)

# Get the order of expression for both types
mlps.high  <- order(lps.macs.meds$log2FC, decreasing = T)
gmlps.high <- order(lps.macs.meds$log2FC, decreasing = F)

# Plot them
pheatmap(lps.macs.meds[mlps.high[ 1:50],1:2],cluster_rows = T, cluster_cols = F, show_rownames = F)
pheatmap(lps.macs.meds[gmlps.high[1:50],1:2],cluster_rows = T, cluster_cols = F, show_rownames = F)

## Zoom in to the top |log2FC|> 2
# Fetch the top 
mlps.high.macs  <- lps.macs.meds[mlps.high[ 1:sum(lps.macs.meds$log2FC > 2)],]
gmlps.high.macs <- lps.macs.meds[gmlps.high[1:sum(lps.macs.meds$log2FC < -2)],]

# Order the top
mlps.high.macs  <- mlps.high.macs[ order( mlps.high.macs$Mlps,  decreasing = T), ]
gmlps.high.macs <- gmlps.high.macs[order(gmlps.high.macs$GMlps, decreasing = T), ]

# Keep only the expressed top (RPKM > 1)
mlps.high.macs  <- mlps.high.macs[  mlps.high.macs$Mlps  > 1, ]
gmlps.high.macs <- gmlps.high.macs[gmlps.high.macs$GMlps > 1, ]

# Retrieve the gene names
mlps.high.genes <- lps.macs.genes[row.names(lps.macs.genes) %in% row.names(mlps.high.macs),, drop = F]
mlps.high.genes <- mlps.high.genes[row.names(mlps.high.macs),]

gmlps.high.genes <- lps.macs.genes[row.names(lps.macs.genes) %in% row.names(gmlps.high.macs),, drop = F]
gmlps.high.genes <- gmlps.high.genes[row.names(gmlps.high.macs),]


# Plot the top
pheatmap(mlps.high.macs[,1:2], color = colorSpacer(endcolor = "dodgerblue", startcolor = "grey", steps = length(mlps.high.genes), return.colors = T), cellwidth = 30, cellheight = 10, border_color = NA, cluster_rows = F, cluster_cols = F, show_rownames = T, angle_col = 45, labels_row = mlps.high.genes, filename = "lpsmac/mlpstopFC.rpkm.pdf")
pheatmap(gmlps.high.macs[,1:2], color = colorSpacer(endcolor = "goldenrod", startcolor = "grey", steps = length(gmlps.high.genes), return.colors = T), cellwidth = 30, cellheight = 10, border_color = NA, cluster_rows = F, cluster_cols = F, show_rownames = T, angle_col = 45, labels_row = gmlps.high.genes, filename = "lpsmac/gmlpstopFC.rpkm.pdf")

# Fix the FC
mlps.high.macs$log2FC[!is.finite( mlps.high.macs$log2FC)] <- max( mlps.high.macs$log2FC[is.finite( mlps.high.macs$log2FC)]) + 1

gmlps.high.macs$log2FC[!is.finite(gmlps.high.macs$log2FC)] <- min(gmlps.high.macs$log2FC[is.finite(gmlps.high.macs$log2FC)]) - 1

# Invert FC for ease
gmlps.high.macs$log2FC <- abs(gmlps.high.macs$log2FC)

# Plot the FC
pheatmap(mlps.high.macs$log2FC, 
         color = colorSpacer(endcolor = "black",middlecolors = "red", startcolor = "grey", steps = length(mlps.high.genes), return.colors = T),
         cellwidth = 30, cellheight = 10, border_color = NA, cluster_rows = F, cluster_cols = F, show_rownames = T, angle_col = 45, 
         labels_row = mlps.high.genes,
         legend_breaks = c(min(mlps.high.macs$log2FC[is.finite(mlps.high.macs$log2FC)]), 4, 6, 8, 10, 12, 14, max(mlps.high.macs$log2FC[is.finite(mlps.high.macs$log2FC)])),
         legend_labels = c(2,4,6,8,10,12,14, "Inf"), filename = "lpsmac/mlpstopFC.log2FC.pdf")

pheatmap(gmlps.high.macs$log2FC, 
         color = colorSpacer(endcolor = "black",middlecolors = "red", startcolor = "grey", steps = length(gmlps.high.genes), return.colors = T),
         cellwidth = 30, cellheight = 10, border_color = NA, cluster_rows = F, cluster_cols = F, show_rownames = T, angle_col = 45, 
         labels_row = gmlps.high.genes,
         legend_breaks = c(min(gmlps.high.macs$log2FC[is.finite(gmlps.high.macs$log2FC)]), 4, 6, 8, 10, 12, 14, max(gmlps.high.macs$log2FC[is.finite(gmlps.high.macs$log2FC)])),
         legend_labels = c(2,4,6,8,10,12,14, "Inf"), filename = "lpsmac/gmlpstopFC.log2FC.pdf")

