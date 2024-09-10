#=================================================================================================================================
## DecoupleR R version code -- Slower than Python, but there is a sh stopping bug in the python version as of 19-09-2023. 
##================================================================================================================================
## TF Network activity
# Download TF network
net <- get_collectri(organism='human', split_complexes=FALSE)
net

# Run ulm
acts <- run_ulm(mat=as.matrix(GetAssayData(final.pop.call.integrated.mye.velocyto.seurat, assay = "RNA")), net=net, .source='source', .target='target',.mor='mor', minsize = 5)
acts

# Extract ulm and store it in tfsulm in the SEurat object
final.pop.call.integrated.mye.velocyto.seurat[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = final.pop.call.integrated.mye.velocyto.seurat) <- "tfsulm"

# Scale the data
final.pop.call.integrated.mye.velocyto.seurat <- ScaleData(final.pop.call.integrated.mye.velocyto.seurat)
final.pop.call.integrated.mye.velocyto.seurat@assays$tfsulm@data <- final.pop.call.integrated.mye.velocyto.seurat@assays$tfsulm@scale.data

# Plot an example
FeaturePlot(final.pop.call.integrated.mye.velocyto.seurat, features = c("EGR1"), pt.size = 3) +
  scale_colour_gradient2(low = 'Darkblue', mid = 'orangered', high = 'yellow') +
  ggtitle('EGR1 activity')

## Plot a heatmap of top diff TFs
n_tfs <- 20

# Extract activities from object as a long dataframe
df <- t(as.matrix(final.pop.call.integrated.mye.velocyto.seurat@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = final.pop.call.integrated.mye.velocyto.seurat$archetype) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  dplyr::summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  group_by(source) %>%
  dplyr::summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Alt, get the highest activty per cluster
tfs <- df %>% group_by(cluster) %>% top_n(5, mean) %>% pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks, cellwidth = 15, cellheight = 20, filename = "decoupler/TF activity.pdf") 

for(i in tfs){
  FeaturePlot(final.pop.call.integrated.mye.velocyto.seurat, features = i, pt.size = 3) +
    scale_colour_gradient2(low = 'blue4', mid = 'grey90', high = 'red') +
    ggtitle(paste(i, "activity", sep = " "))
  ggsave(paste("decoupler/", i, " activity UMAP.pdf"))
}

## Pathway activity
# Download pathways
net <- get_progeny(organism = 'human', top = 500)

# Run ulm
acts <- run_ulm(mat=as.matrix(GetAssayData(final.pop.call.integrated.mye.velocyto.seurat, assay = "RNA")), net=net, .source='source', .target='target',.mor='weight', minsize = 5)
acts

# Extract mlm and store it in pathwaysmlm in the SEurat object
final.pop.call.integrated.mye.velocyto.seurat[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = final.pop.call.integrated.mye.velocyto.seurat) <- "pathwaysmlm"

# Scale the data
final.pop.call.integrated.mye.velocyto.seurat <- ScaleData(final.pop.call.integrated.mye.velocyto.seurat)
final.pop.call.integrated.mye.velocyto.seurat@assays$pathwaysmlm@data <- final.pop.call.integrated.mye.velocyto.seurat@assays$pathwaysmlm@scale.data

# Plot an example
FeaturePlot(final.pop.call.integrated.mye.velocyto.seurat, features = c("NFkB"), pt.size = 3) +
  scale_colour_gradient2(low = 'blue4', mid = 'grey90', high = 'red') +
  ggtitle('NFkB activity')

## Plot a heatmap of pathways
# Extract activities from object as a long dataframe
df <- t(as.matrix(final.pop.call.integrated.mye.velocyto.seurat@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = final.pop.call.integrated.mye.velocyto.seurat$archetype) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  dplyr::summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks, cellwidth = 15, cellheight = 20, filename = "decoupler/pathway activity.pdf") 

for(i in Features(final.pop.call.integrated.mye.velocyto.seurat)){
  FeaturePlot(final.pop.call.integrated.mye.velocyto.seurat, features = i, pt.size = 3) +
    scale_colour_gradient2(low = 'blue4', mid = 'grey90', high = 'red') +
    ggtitle(paste(i, "activity", sep = " "))
  ggsave(paste("decoupler/", i, " activity UMAP.pdf"))
}

