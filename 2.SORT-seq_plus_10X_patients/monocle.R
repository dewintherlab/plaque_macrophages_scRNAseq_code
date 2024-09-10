#=========================================================================================================================
## Calculate cell trajectories using monocle 3
##========================================================================================================================
dir.create("monocle_results", showWarnings = F)

## All Myeloid cells
##========================================================================================================================
## Import the seurat object
#  Convert to monocle cell data set object using SeuratWrappers
mye.patients.monocle <- as.cell_data_set(mye.patients.seurat, assay = "RNA", group.by = "method.int.idents")

## Build the single-cell trajectories
# Cluster cells with monocle
mye.patients.monocle <- cluster_cells(mye.patients.monocle, reduction_method = "UMAP", resolution = 0.01)

p1 <- plot_cells(mye.patients.monocle, show_trajectory_graph = FALSE, cell_size = 2)
p2 <- plot_cells(mye.patients.monocle, show_trajectory_graph = FALSE, cell_size = 2, color_cells_by = "partition")
p3 <- plot_cells(mye.patients.monocle, show_trajectory_graph = FALSE, cell_size = 2, color_cells_by = "method.int.idents", label_cell_groups = F) + 
  scale_colour_manual(breaks =  names(M.int_refined.pop.colors), values = M.int_refined.pop.colors)
wrap_plots(p1, p2, p3)
ggsave(file = "monocle_results/monocle_clusters.UMAP.pdf", width = 30, height = 10)

# Determine trajectories
mye.patients.monocle <- learn_graph(mye.patients.monocle, verbose = T, use_partition = F, learn_graph_control = list(orthogonal_proj_tip = T,
                                                                                                                     rann.k = NULL, 
                                                                                                                     prune_graph = T,
                                                                                                                     minimal_branch_len = 8,
                                                                                                                     geodesic_distance_ratio = 1/3,
                                                                                                                     euclidean_distance_ratio = 3))

# plot by ident
plot_cells(mye.patients.monocle, 
           label_groups_by_cluster = T, 
           label_leaves = T, 
           label_branch_points = T,
           color_cells_by = "method.int.idents",
           label_cell_groups = F, label_roots = T, cell_size = 4, trajectory_graph_color = "red", trajectory_graph_segment_size = 2) + 
  scale_colour_manual(breaks =  names(M.int_refined.pop.colors), values = M.int_refined.pop.colors)
ggsave(file = "monocle_results/monocle_trajectory.idents.pdf", width = 20, height = 10)

# plot by monocle cluster
plot_cells(mye.patients.monocle, 
           label_groups_by_cluster = T, 
           label_leaves = T, 
           label_branch_points = T,
           label_cell_groups = F, label_roots = T, cell_size = 1.5)
ggsave(file = "monocle_results/monocle_trajectory.clusters.pdf", width = 20, height = 10)

## Plot by pseudotime
# Set cells from classical monocyte idents to root
clas.mono.cells <- WhichCells(mye.patients.seurat, 
                              idents = unique(Idents(mye.patients.seurat))[grep("Clas", unique(Idents(mye.patients.seurat)))]
)
mye.patients.monocle <- order_cells(mye.patients.monocle, root_cells = clas.mono.cells)

plot_cells(mye.patients.monocle, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_roots = F,
           label_branch_points = FALSE, cell_size = 2, show_trajectory_graph = FALSE)
ggsave(file = "monocle_results/monocle_pseudotime.clas_mon_root.pdf")

# Set principle node as root
plot_cells(mye.patients.monocle, 
           label_principal_points = T)
ggsave(file = "monocle_results/monocle_trajectory_princple_points.pdf", limitsize = F)

mye.patients.monocle <- order_cells(mye.patients.monocle, root_pr_nodes = "Y_1")

plot_cells(mye.patients.monocle, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_roots = F,
           label_branch_points = FALSE, cell_size = 2, show_trajectory_graph = F)
ggsave(file = "monocle_results/monocle_pseudotime.pr_node_root.pdf")


## Monocle diff expression analyses by module
##========================================================================================================================
# First define genes focally expressed in the UMAP
pr_graph_test_res  <- graph_test(cds            = mye.patients.monocle, 
                                 neighbor_graph = "knn", 
                                 cores          = 11)
pr_deg_ids         <- row.names(subset(pr_graph_test_res, q_value < 0.05))

# Group genes into modules based on expression
# Run preprocess_cds to populate some fields we need
mye.patients.monocle <- preprocess_cds(cds = mye.patients.monocle)
gene_module_df       <- find_gene_modules(cds = mye.patients.monocle[pr_deg_ids,], 
                                          resolution = 1e-2)

# Aggregate gene expression per module for visualisation of specificity
cell_group_df      <- tibble(cell      = row.names(colData(mye.patients.monocle)), 
                            cell_group = as.vector(mye.patients.monocle@colData$ident))
agg_mat            <- aggregate_gene_expression(mye.patients.monocle, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ",    row.names(agg_mat))
colnames(agg_mat)  <- colnames(agg_mat)

pheatmap(mat               = agg_mat, 
         cluster_rows      = TRUE, 
         cluster_cols      = TRUE,
         scale             = "column", 
         clustering_method = "ward.D2",
         cellwidth         = 15, 
         cellheight        = 7,
         fontsize          = 6,
         filename          = "monocle_results/module_gene_expression.pdf"  
         )

# Plot genes in the modules
dir.create("monocle_results/module", showWarnings = F, recursive = T)
for(i in levels(gene_module_df$module)){
  cat("Working on module", i, "of", length(levels(gene_module_df$module)), "...\n")
  module_genes <- subset(gene_module_df, module == i)$id
  DoHeatmap(object       = mye.patients.seurat, 
            features     = module_genes, 
            cells        = row.names(colData(mye.patients.monocle)), 
            group.colors = M.int_refined.pop.colors, 
            group.bar    = T, 
            label        = F, 
            raster       = F) + theme(legend.position = "none")
  ggsave(filename = paste("monocle_results/module/Module ", i, " genes in Myeloid clusters.pdf", sep = ""), width = 20, height = (3 + (0.1 * length(module_genes))), limitsize = F)
}

# Save genes in the modules
write.table(x = gene_module_df, file = "monocle_results/module/module_genes.txt", quote = F, row.names = F, sep = "\t")

# Plot pathways for genes in the modules
ont.M.modules <- list()
for(i in levels(gene_module_df$module)){
  cat("Working on module", i, "of", length(levels(gene_module_df$module)), "...\n")
  df                <- data.frame(gene = subset(gene_module_df, module == i)$id, avg_log2FC = 0)
  ont.M.modules[[i]] <- get_ontology(res          = df, 
                                    name         = i, 
                                    outdir       = "monocle_results/module/", 
                                    full_GSEA    = F, 
                                    return.data  = T, 
                                    plot.top.n   = 20, 
                                    universe     = mye.patients.seurat, 
                                    volcano.plot = F)
}
  
  
## Monocle diff expression analyses by trajectory
##========================================================================================================================
# Calculate genes changing over pseudotime
pr_test_res <- graph_test(cds            = mye.patients.monocle, 
                          neighbor_graph = "principal_graph",
                          cores          = 11
                          )
pr_test_res <- pr_test_res[order(pr_test_res$q_value),]
pr_deg_ids  <- row.names(subset(pr_test_res, q_value < 0.05))

# Plot some top hits
pr_test_res[1:50,]
rowData(mye.patients.monocle)$gene_name       <- row.names(mye.patients.monocle)
rowData(mye.patients.monocle)$gene_short_name <- rowData(mye.patients.monocle)$gene_name
plot_cells(mye.patients.monocle, genes=c("C1QA", "S100A8", "FCER1G", "FCGR3A", "IL1B", "CXCL2", "EGR1", "CD14", "CSF1R"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

# Define modules
mye.patients.monocle <- preprocess_cds(cds = mye.patients.monocle)
gene_module_df  <- find_gene_modules(cds = mye.patients.monocle[pr_deg_ids,], 
                                     resolution = c(10^seq(-6,-1)))

# Aggregate gene expression per module for visualisation of specificity
cell_group_df      <- tibble(cell      = row.names(colData(mye.patients.monocle)), 
                             cell_group = as.vector(mye.patients.monocle@colData$ident))
agg_mat            <- aggregate_gene_expression(mye.patients.monocle, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ",    row.names(agg_mat))
colnames(agg_mat)  <- colnames(agg_mat)

pheatmap(mat               = agg_mat, 
         cluster_rows      = TRUE, 
         cluster_cols      = TRUE,
         scale             = "column", 
         clustering_method = "ward.D2",
         cellwidth         = 15, 
         cellheight        = 7,
         fontsize          = 6,
         filename          = "monocle_results/trajectory_gene_expression.pdf"  
)

# Plot genes in the trajectory
dir.create("monocle_results/trajectory", showWarnings = F, recursive = T)
for(i in levels(gene_module_df$module)){
  cat("Working on module", i, "of", length(levels(gene_module_df$module)), "...\n")
  module_genes <- subset(gene_module_df, module == i)$id
  DoHeatmap(object       = mye.patients.seurat, 
            features     = module_genes, 
            cells        = row.names(colData(mye.patients.monocle)), 
            group.colors = M.int_refined.pop.colors, 
            group.bar    = T, 
            label        = F, 
            raster       = F) + theme(legend.position = "none")
  ggsave(filename = paste("monocle_results/trajectory/Module ", i, " genes in Myeloid clusters.pdf", sep = ""), width = 20, height = (3 + (0.1 * length(module_genes))), limitsize = F)
}

# Save genes in the trajectory
write.table(x = gene_module_df, file = "monocle_results/trajectory/module_genes.txt", quote = F, row.names = F, sep = "\t")

# Plot pathways for genes in the trajectories
ont.M.modules <- list()
for(i in levels(gene_module_df$module)){
  cat("Working on module", i, "of", length(levels(gene_module_df$module)), "...\n")
  df                <- data.frame(gene = subset(gene_module_df, module == i)$id, avg_log2FC = 0)
  ont.M.modules[[i]] <- get_ontology(res          = df, 
                                     name         = i, 
                                     outdir       = "monocle_results/trajectory/", 
                                     full_GSEA    = F, 
                                     return.data  = T, 
                                     plot.top.n   = 20, 
                                     universe     = mye.patients.seurat, 
                                     volcano.plot = F)
}

