#-----------------------------------------------------------------------------------------
#Setup data

#Make directory for results
dir.create("results", showWarnings = FALSE)
dir.create("Seurat Objects", showWarnings = FALSE)

# merge seurat objects from list
all.seur.combined <- merge(seurat.object.list[[1]], y = seurat.object.list[2:length(seurat.object.list)])
all.seur.combined

# HOUSEKEEPING: Remove the list objects to save memory
rm(current.df)
rm(current.meta.data)
rm(current.object)
rm(seurat.object.list)

## Filter cells with too low and too high gene counts
# pick thresholds
VlnPlot(object = all.seur.combined, features = c("nFeature_RNA","nCount_RNA"), y.max = 10000) + theme(aspect.ratio = 1/4)
ggsave("results/number of transcripts.pdf", width = 20, height = 20)

all.seur.combined <- subset(all.seur.combined, subset = nFeature_RNA > 500 & nFeature_RNA < 10000)
all.seur.combined

## Filter and regress out the MT reads and apoptotic genes
# store mitochondrial percentage in object meta data
all.seur.combined <- PercentageFeatureSet(all.seur.combined, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(object = all.seur.combined, features = c("percent.mt"))
ggsave("results/percent.mt.pdf")

# store apoptotic gene percentage in object meta data
all.seur.combined <- PercentageFeatureSet(all.seur.combined, pattern = "KCNQ1OT1", col.name = "percent.KCNQ1OT1")
all.seur.combined <- PercentageFeatureSet(all.seur.combined, pattern = "UGDH-AS1", col.name = "percent.UGDH.AS1")
all.seur.combined <- PercentageFeatureSet(all.seur.combined, pattern = "GHET1", col.name = "percent.GHET1")
all.seur.combined <- PercentageFeatureSet(all.seur.combined, pattern = "ROR1-AS1", col.name = "percent.ROR1.AS1")
VlnPlot(object = all.seur.combined, features = c("percent.KCNQ1OT1"))
ggsave("results/percent.KCN.pdf")

# Save the object
saveRDS(all.seur.combined, "Seurat Objects/main.seurat.prefilter.RDS")