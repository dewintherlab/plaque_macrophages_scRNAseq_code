#=========================================================================================================================
## Subset the refined macs and monos form the full integrated set
##========================================================================================================================
dir.create("from_full_myeloid_43p_celseq_integration/")

# Subset the object
from_full.integrated.mye.seurat <- subset(integrated.full.seurat, cells = WhichCells(integrated.full.seurat, idents = names(M.int_refined.pop.colors)))
from_full.integrated.mye.seurat
Idents(from_full.integrated.mye.seurat) <- factor(x = Idents(from_full.integrated.mye.seurat), levels = names(M.int_refined.pop.colors))

# Plot it! (to verify)
customUMAP(object = from_full.integrated.mye.seurat, cols = M.int_refined.pop.colors, 
           file.name = "from_full_myeloid_43p_celseq_integration/Initial Population UMAP.pdf", plot.width = 15, legend.pos = "right")

# Rescale
from_full.integrated.mye.seurat <- ScaleData(from_full.integrated.mye.seurat, verbose = F)
from_full.integrated.mye.seurat <- RunPCA(from_full.integrated.mye.seurat, verbose = F)
from_full.integrated.mye.seurat <- RunUMAP(from_full.integrated.mye.seurat, verbose = F, reduction = "pca", dims = 1:30)

# Plot it! (to verify)
customUMAP(object = from_full.integrated.mye.seurat, cols = M.int_refined.pop.colors, 
           file.name = "from_full_myeloid_43p_celseq_integration/Population UMAP.pdf", plot.width = 20, legend.pos = "right", pt.size = 3, shuffle = T)
customUMAP(object = from_full.integrated.mye.seurat, cols = M.int_refined.pop.colors, 
           file.name = "from_full_myeloid_43p_celseq_integration/Population UMAP no legend.pdf", plot.width = 10, legend.pos = "none", pt.size = 3, shuffle = T)


## Run monocle
##========================================================================================================================
## Import the seurat object
#  Convert to monocle cell data set object using SeuratWrappers
from_full.mye.int.monocle <- as.cell_data_set(from_full.integrated.mye.seurat, assay = "RNA")

## Build the single-cell trajectories
# Cluster cells with monocle
from_full.mye.int.monocle <- cluster_cells(from_full.mye.int.monocle, reduction_method = "UMAP", resolution = 0.01)

p1 <- plot_cells(from_full.mye.int.monocle, show_trajectory_graph = FALSE, cell_size = 2)
p2 <- plot_cells(from_full.mye.int.monocle, show_trajectory_graph = FALSE, cell_size = 2, color_cells_by = "partition")
p3 <- plot_cells(from_full.mye.int.monocle, show_trajectory_graph = FALSE, cell_size = 2, color_cells_by = "ident", label_cell_groups = F) + 
  scale_colour_manual(breaks =  names(M.int_refined.pop.colors), values = M.int_refined.pop.colors)
wrap_plots(p1, p2, p3)
ggsave(file = "from_full_myeloid_43p_celseq_integration/monocle_clusters.UMAP.pdf", width = 30, height = 10)

# Determine trajectories
from_full.mye.int.monocle <- learn_graph(from_full.mye.int.monocle, verbose = T, use_partition = T, learn_graph_control = list(orthogonal_proj_tip = T,
                                                                                                                               rann.k = NULL, 
                                                                                                                               prune_graph = T,
                                                                                                                               minimal_branch_len = 8,
                                                                                                                               geodesic_distance_ratio = 1/3,
                                                                                                                               euclidean_distance_ratio = 3))

# plot by ident
plot_cells(from_full.mye.int.monocle, 
           label_groups_by_cluster = T, 
           label_leaves = T, 
           label_branch_points = T,
           color_cells_by = "ident", trajectory_graph_color = "red", trajectory_graph_segment_size = 1.5,
           label_cell_groups = F, label_roots = T, cell_size = 3) + 
  scale_colour_manual(breaks =  names(M.int_refined.pop.colors), values = M.int_refined.pop.colors)
ggsave(file = "from_full_myeloid_43p_celseq_integration/monocle_trajectory.idents.pdf", width = 20, height = 10)

# plot by monocle cluster
plot_cells(from_full.mye.int.monocle, 
           label_groups_by_cluster = T, 
           label_leaves = T, 
           label_branch_points = T,
           label_cell_groups = F, label_roots = T, cell_size = 1.5)
ggsave(file = "from_full_myeloid_43p_celseq_integration/monocle_trajectory.clusters.pdf", width = 20, height = 10)

## Plot by pseudotime
# Set cells from classical monocyte idents to root
clas.mono.cells <- WhichCells(from_full.mye.int.seurat, 
                              idents = unique(Idents(from_full.mye.int.seurat))[grep("Clas", unique(Idents(from_full.mye.int.seurat)))]
)
from_full.mye.int.monocle <- order_cells(from_full.mye.int.monocle, root_cells = clas.mono.cells)

plot_cells(from_full.mye.int.monocle, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_roots = F,
           label_branch_points = FALSE, cell_size = 2, show_trajectory_graph = FALSE)
ggsave(file = "from_full_myeloid_43p_celseq_integration/monocle_pseudotime.clas_mon_root.pdf")

# Set principle node as root
plot_cells(from_full.mye.int.monocle, 
           label_principal_points = T)
ggsave(file = "from_full_myeloid_43p_celseq_integration/monocle_trajectory_princple_points.pdf", limitsize = F)

from_full.mye.int.monocle <- order_cells(from_full.mye.int.monocle, root_pr_nodes = "Y_274")

plot_cells(from_full.mye.int.monocle, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_roots = F,
           label_branch_points = FALSE, cell_size = 1.5, show_trajectory_graph = F)
ggsave(file = "from_full_myeloid_43p_celseq_integration/monocle_pseudotime.pr_node_root.pdf")
