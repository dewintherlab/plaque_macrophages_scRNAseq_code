#==========================================================================================
## Run cellrank trajectory analyses
#==========================================================================================

## Import modules
import sys

if "google.colab" in sys.modules:
  !pip install -q git+https://github.com/theislab/cellrank@dev
  !pip install python-igraph

import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
from matplotlib import pyplot as plt

scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2

import warnings

warnings.simplefilter("ignore", category=UserWarning)
warnings.simplefilter("ignore", category=FutureWarning)
warnings.simplefilter("ignore", category=DeprecationWarning)

## Load data mac + mono
adata = sc.read_h5ad("Seurat_Objects/myeloid.velocyto.refined.h5ad")
adata.obs['final.pop.idents'] = adata.obs['final.pop.idents'].astype('category')
scv.pl.proportions(adata, groupby='final.pop.idents',save="CellRank/Splice Ratio.pdf")
adata

## Set vars
true_label=adata.obs.loc[:,"final.pop.idents"]
random_seed=666


# Plot UMAP to check
sc.pl.umap(adata, color="final.pop.idents", palette=pop_cols, return_fig=True)
#plt.show()
plt.savefig("CellRank/UMAP.pdf", bbox_inches='tight')

#Pre-process
scv.pp.filter_and_normalize(adata, n_top_genes=2000, subset_highly_variable=False, min_shared_counts=20)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=50)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

# Run scvelo
scv.tl.recover_dynamics(adata, n_jobs=5)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

# Visualize
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_fontsize=12, title="", density=3, alpha= 1, size=3, smooth=0.8,n_neighbors=200, min_mass=1, show=False, color='final.pop.idents', palette=pop_cols, legend_loc='right margin')
plt.savefig("CellRank/velocyto.png", bbox_inches='tight')

# Find terminal states
cr.tl.terminal_states(adata, weight_connectivities=0.2, cluster_key='final.pop.idents')
cr.pl.terminal_states(adata,save="CellRank/terminal states.pdf")
cr.pl.terminal_states(adata, discrete = True, save="CellRank/auto terminal states discrete.pdf")
cr.pl.terminal_states(adata)

# Find initial states
cr.tl.initial_states(adata, cluster_key="final.pop.idents")
cr.pl.initial_states(adata, discrete=True, save="CellRank/initial states discrete.pdf")
cr.pl.initial_states(adata, discrete=False, save="CellRank/initial states.pdf")

# Compute lineages and plot fate maps
cr.tl.lineages(adata)
cr.pl.lineages(adata, same_plot=False)


#==========================================================================================
## Run customized cellrank in order to set initial and terminal states properley
#==========================================================================================
## Import kernels and compute the trnasition matrixes
from cellrank.tl.kernels import VelocityKernel
from cellrank.tl.kernels import ConnectivityKernel

vk = VelocityKernel(adata).compute_transition_matrix()
ck = ConnectivityKernel(adata).compute_transition_matrix()

## combine kernels
combined_kernel = 0.8 * vk + 0.2 * ck

## Add the GPCCA estimator
from cellrank.tl.estimators import GPCCA

g = GPCCA(combined_kernel)

## Compute matrix decomposition
g.compute_schur(n_components=20)
g.plot_spectrum()
plt.show()

# Compute macrostates
g.compute_macrostates(cluster_key="final.pop.idents", n_states=4)
g.plot_macrostates(discrete=True)
plt.savefig("CellRank/macrostates discrete.pdf", bbox_inches='tight')
g.plot_macrostates()
plt.savefig("CellRank/macrostates.pdf", bbox_inches='tight')

# Set terminal states
g.set_terminal_states_from_macrostates(names = ["CD14+CD16+CD64-SLAN+ Non-Classical Monocytes", "CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages"], cluster_key="final.pop.idents")
g.plot_terminal_states(discrete=True)
plt.savefig("CellRank/terminal states discrete.pdf", bbox_inches='tight')
g.plot_terminal_states()
plt.savefig("CellRank/terminal states.pdf", bbox_inches='tight')

# Set initial states
g._set_initial_states_from_macrostates(names = ["CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages", "CD14+CD16-CD64+SLAN- Classical Monocytes"])

# Compute fate probalities
g.compute_absorption_probabilities()

# Compute lineage drivers
g.compute_lineage_drivers()

g.plot_lineage_drivers(lineage="CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages", n_genes=9, ncols=3)
plt.show()

g.plot_lineage_drivers(lineage="CD14+CD16+CD64-SLAN+ Non-Classical Monocytes", n_genes=9, ncols=3)
plt.show()

## Now visualize using 'normal' cr methods
# Compute fate maps
cr.tl.lineages(adata)
cr.pl.lineages(adata, same_plot=False)
plt.savefig("CellRank/lineages separate.pdf", bbox_inches='tight')
cr.pl.lineages(adata, same_plot=True)
plt.savefig("CellRank/lineages combined.pdf", bbox_inches='tight')

# Compute Latent time
scv.tl.recover_latent_time(adata, root_key="initial_states_probs", end_key="terminal_states_probs")
scv.tl.paga(
    adata,
    groups="final.pop.idents",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="velocity_pseudotime")
    
cr.pl.cluster_fates(
    adata,
    mode="paga_pie",
    cluster_key="final.pop.idents",
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    legend_loc="top left out",
    node_size_scale=5,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA"
)
plt.show()
plt.savefig("CellRank/cell fate directed PAGA.pdf", bbox_inches='tight')

cr.pl.cluster_fates(
    adata,
    mode="violin",
    cluster_key="final.pop.idents",
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    legend_loc="top left out",
    node_size_scale=5,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA"
)
plt.show()
plt.savefig("CellRank/cell fate violin.pdf", bbox_inches='tight')

cr.pl.cluster_fates(
    adata,
    mode="heatmap",
    cluster_key="final.pop.idents",
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    title="directed PAGA"
)
plt.show()
plt.savefig("CellRank/cell fate heatmap.pdf", bbox_inches='tight')

cr.pl.cluster_fates(adata, mode="bar", cluster_key='final.pop.idents', ncols = 3, figsize=(15,25))
plt.show()
plt.savefig("CellRank/cell fate barplot.pdf", bbox_inches='tight')

# Find lineage drivers
cr.tl.lineage_drivers(adata)

# Plot top 5 per lineage
cr.pl.lineage_drivers(adata, lineage="CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages", n_genes=9, ncols=3)
plt.show()
plt.savefig("CellRank/res lineage drivers.pdf", bbox_inches='tight')



cr.pl.lineage_drivers(adata, lineage="CD14+CD16+CD64-SLAN+ Non-Classical Monocytes", n_genes=9, ncols=3)
plt.show()
plt.savefig("CellRank/non clas lineage drivers.pdf", bbox_inches='tight')

cr.pl.lineage_drivers(adata, lineage="CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages", n_genes=9, ncols=3)
plt.show()
plt.savefig("CellRank/foam lineage drivers.pdf", bbox_inches='tight')

## Gene expression trends
# From Clas mons
root_idx = np.where(adata.obs["initial_states"] == "CD14+CD16-CD64+SLAN- Classical Monocytes")[0][0]
adata.uns["iroot"] = root_idx
sc.tl.dpt(adata)

scv.pl.scatter(
    adata,
    color=["final.pop.idents", root_idx, "latent_time", "dpt_pseudotime"],
    fontsize=16,
    cmap="viridis",
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    title=["clusters", "root cell", "latent time", "dpt pseudotime"],save="CellRank/clas pseudotime.pdf"
)

model = cr.ul.models.GAM(adata)
cr.pl.gene_trends(
    adata,
    model=model,
    data_key="X",
    genes=["PLIN2", "TIMP1", "ABCA1", "FOLR2", "CD9", "PLTP", "FCGR3A", "CDKN1C", "IFITM1"],
    ncols=3,
    time_key="latent_time",
    legend_loc="none",
    same_plot=True,
    hide_cells=True,
    figsize=(15, 15),
    n_test_points=200,
)
plt.show()
plt.savefig("CellRank/clas gene expression trends some genes latent time.pdf", bbox_inches='tight')

cr.pl.gene_trends(
    adata,
    model=model,
    data_key="X",
    genes=["PLIN2", "TIMP1", "ABCA1", "FOLR2", "CD9", "PLTP", "FCGR3A", "CDKN1C", "IFITM1"],
    ncols=3,
    time_key="dpt_pseudotime",
    legend_loc="none",
    same_plot=True,
    hide_cells=True,
    figsize=(15, 15),
    n_test_points=200,
)
plt.show()
plt.savefig("CellRank/clas gene expression trends some genes dpt.pdf", bbox_inches='tight')

cr.pl.heatmap(
    adata,
    model,
    genes=adata.varm['terminal_lineage_drivers']["CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages_corr"].sort_values(ascending=False).index[:100],
    show_absorption_probabilities=True,
    lineages="CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages",
    n_jobs=1,
    backend="loky",
)
plt.show()
plt.savefig("CellRank/clas foam gene expression trends heatmap.pdf", bbox_inches='tight')

cr.pl.heatmap(
    adata,
    model,
    genes=adata.varm['terminal_lineage_drivers']["CD14+CD16+CD64-SLAN+ Non-Classical Monocytes_corr"].sort_values(ascending=False).index[:100],
    show_absorption_probabilities=True,
    lineages="CD14+CD16+CD64-SLAN+ Non-Classical Monocytes",
    n_jobs=1,
    backend="loky",
)
plt.show()
plt.savefig("CellRank/clas non clas gene expression trends heatmap.pdf", bbox_inches='tight')

# From res
root_idx = np.where(adata.obs["initial_states"] == "CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages")[0][0]
adata.uns["iroot"] = root_idx
sc.tl.dpt(adata)

scv.pl.scatter(
    adata,
    color=["final.pop.idents", root_idx, "latent_time", "dpt_pseudotime"],
    fontsize=16,
    cmap="viridis",
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    title=["clusters", "root cell", "latent time", "dpt pseudotime"],save="CellRank/res pseudotime.pdf"
)

model = cr.ul.models.GAM(adata)
cr.pl.gene_trends(
    adata,
    model=model,
    data_key="X",
    genes=["PLIN2", "TIMP1", "ABCA1", "FOLR2", "CD9", "PLTP", "FCGR3A", "CDKN1C", "IFITM1"],
    ncols=3,
    time_key="dpt_pseudotime",
    legend_loc="none",
    same_plot=True,
    hide_cells=True,
    figsize=(15, 15),
    n_test_points=200,
)
plt.show()
plt.savefig("CellRank/res gene expression trends some genes dpt.pdf", bbox_inches='tight')

cr.pl.heatmap(
    adata,
    model,
    genes=adata.varm['terminal_lineage_drivers']["CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages_corr"].sort_values(ascending=False).index[:100],
    show_absorption_probabilities=True,
    lineages="CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages",
    n_jobs=1,
    backend="loky",
)
plt.show()
plt.savefig("CellRank/res foam gene expression trends heatmap.pdf", bbox_inches='tight')

# PLot some random walks
vk.plot_random_walks(
    50,
    start_ixs={"final.pop.idents": "CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages"},
    max_iter=100,
    successive_hits=5,
    show_progress_bar=False,
    cmap="viridis",linealpha=0.1,
    seed=420,
    ixs_legend_loc="best"
)
plt.show()
plt.savefig("CellRank/random walks res start.pdf", bbox_inches='tight')

vk.plot_random_walks(
    50,
    start_ixs={"final.pop.idents": "CD14+CD16-CD64+SLAN- Classical Monocytes"},
    max_iter=100,
    successive_hits=5,
    show_progress_bar=False,
    cmap="viridis",linealpha=0.1,
    seed=420,
    ixs_legend_loc="best"
)
plt.show()
plt.savefig("CellRank/random walks clas mon start.pdf", bbox_inches='tight')

vk.plot_random_walks(
    50,
    start_ixs={"final.pop.idents": "CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages"},
    stop_ixs={"final.pop.idents": "CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages"},
    max_iter=100,
    successive_hits=5,
    show_progress_bar=False,
    cmap="viridis",linealpha=0.1,
    seed=420,
    ixs_legend_loc="best"
)
plt.show()
plt.savefig("CellRank/random walks res start foam stop.pdf", bbox_inches='tight')

vk.plot_random_walks(
    50,
    start_ixs={"final.pop.idents": "CD14+CD16-CD64+SLAN- Classical Monocytes"},
    stop_ixs={"final.pop.idents": "CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages"},
     max_iter=200,
    successive_hits=5,
    show_progress_bar=False,
    cmap="viridis",linealpha=0.1,
    seed=420,
    ixs_legend_loc="best"
)
plt.show()
plt.savefig("CellRank/random walks clas mon start foam stop.pdf", bbox_inches='tight')


#==========================================================================================
## Run cellrank with cyto trace kernel
#==========================================================================================
from cellrank.tl.kernels import CytoTRACEKernel
ctk = CytoTRACEKernel(adata)

scv.pl.scatter(
    adata,
    c=["ct_pseudotime", "final.pop.idents"],
    legend_loc="right",
    color_map="gnuplot2",save="CellRank/cytotrace pseudotime.pdf"
)
sc.pl.violin(adata, keys=["ct_pseudotime"], groupby="final.pop.idents", rotation=90)

# Compute transition matrix
ctk.compute_transition_matrix(threshold_scheme="soft", nu=0.5)

# Visualize
ctk.compute_projection()
scv.pl.velocity_embedding_stream(adata, color="final.pop.idents", vkey="T_fwd", legend_loc="right")
plt.savefig("CellRank/cytotrace velo.png", bbox_inches='tight')

#==========================================================================================
## Run cellrank without monocytes
#==========================================================================================
# Remove mono cells
adata = adata[adata.obs['final.pop.idents'].str.find('CD64') == -1].copy()

# Plot UMAP to check
sc.pl.umap(adata, color="final.pop.idents", palette=pop_cols, return_fig=True)
#plt.show()
plt.savefig("CellRank/UMAP.pdf", bbox_inches='tight')

#Pre-process
scv.pp.filter_and_normalize(adata, n_top_genes=2000, subset_highly_variable=False, min_shared_counts=20)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=50)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

# Run scvelo
scv.tl.recover_dynamics(adata, n_jobs=5)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

# Visualize
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_fontsize=12, title="", density=3, alpha= 1, size=3, smooth=0.8,n_neighbors=200, min_mass=1, show=False, color='final.pop.idents', palette=pop_cols, legend_loc='right margin')
plt.savefig("CellRank/velocyto.png", bbox_inches='tight')

# Find terminal states
cr.tl.terminal_states(adata, weight_connectivities=0.2, cluster_key='final.pop.idents', n_states = 2)
cr.pl.terminal_states(adata,save="CellRank/terminal states.pdf")
cr.pl.terminal_states(adata, discrete = True, save="CellRank/auto terminal states discrete.pdf")
cr.pl.terminal_states(adata)

# Find initial states
cr.tl.initial_states(adata, cluster_key="final.pop.idents")
cr.pl.initial_states(adata, discrete=True, save="CellRank/initial states discrete.pdf")
cr.pl.initial_states(adata, discrete=False, save="CellRank/initial states.pdf")

# Compute lineages and plot fate maps
cr.tl.lineages(adata)
cr.pl.lineages(adata, same_plot=False)


#==========================================================================================
## Run customized cellrank in order to set initial and terminal states properley
#==========================================================================================
## Import kernels and compute the trnasition matrixes
from cellrank.tl.kernels import VelocityKernel
from cellrank.tl.kernels import ConnectivityKernel

vk = VelocityKernel(adata).compute_transition_matrix()
ck = ConnectivityKernel(adata).compute_transition_matrix()

## combine kernels
combined_kernel = 0.8 * vk + 0.2 * ck

## Add the GPCCA estimator
from cellrank.tl.estimators import GPCCA

g = GPCCA(combined_kernel)

## Compute matrix decomposition
g.compute_schur(n_components=20)
g.plot_spectrum()
plt.show()

# Compute macrostates
g.compute_macrostates(cluster_key="final.pop.idents", n_states = 3)
g.plot_macrostates(discrete=True, save = "CellRank/macrostates discrete.pdf")
g.plot_macrostates(save = "CellRank/macrostates.pdf")

# Set terminal states
g.set_terminal_states_from_macrostates(names = ["CD14+TREM2+FOLR2-ABCG+ Lipid Associated Macrophages", "CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages"], cluster_key="final.pop.idents")
g.plot_terminal_states(discrete=True, save = "CellRank/terminal states discrete.pdf")
g.plot_terminal_states(save = "CellRank/terminal states.pdf")

# Set initial states
g._set_initial_states_from_macrostates(names = ["CD14+-IL1B+SELL+CD16+ Migrating Inflammatory Monocyte-derived Macrophages"])

# Compute fate probalities
g.compute_absorption_probabilities()

# Compute lineage drivers
g.compute_lineage_drivers()

g.plot_lineage_drivers(lineage="CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages", n_genes=9, ncols=3)
plt.show()

g.plot_lineage_drivers(lineage="CD14+TREM2+FOLR2-ABCG+ Lipid Associated Macrophages", n_genes=9, ncols=3)
plt.show()

## Now visualize using 'normal' cr methods
# Compute fate maps
cr.tl.transition_matrix(adata)
cr.tl.lineages(adata)
cr.pl.lineages(adata, same_plot=False, save = "CellRank/lineages separate.pdf")
cr.pl.lineages(adata, same_plot=True, save = "CellRank/lineages combined.pdf")

# Compute Latent time
scv.tl.recover_latent_time(adata, root_key="initial_states_probs", end_key="terminal_states_probs")
scv.tl.paga(
    adata,
    groups="final.pop.idents",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="velocity_pseudotime")
    
cr.pl.cluster_fates(
    adata,
    mode="paga_pie",
    cluster_key="final.pop.idents",
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    legend_loc="top left out",
    node_size_scale=5,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA"
)
plt.show()

cr.pl.cluster_fates(
    adata,
    mode="violin",
    cluster_key="final.pop.idents",
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    legend_loc="top left out",
    node_size_scale=5,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA"
)
plt.show()
plt.savefig("CellRank/cell fate violin.pdf", bbox_inches='tight')

cr.pl.cluster_fates(
    adata,
    mode="heatmap",
    cluster_key="final.pop.idents",
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    title="directed PAGA"
)
plt.show()
plt.savefig("CellRank/cell fate heatmap.pdf", bbox_inches='tight')

cr.pl.cluster_fates(adata, mode="bar", cluster_key='final.pop.idents', ncols = 3, figsize=(15,25))
plt.show()
plt.savefig("CellRank/cell fate barplot.pdf", bbox_inches='tight')

# Find lineage drivers
cr.tl.lineage_drivers(adata)

# Plot top 5 per lineage
cr.pl.lineage_drivers(adata, lineage="CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages", n_genes=9, ncols=3)
plt.show()



cr.pl.lineage_drivers(adata, lineage="CD14+TREM2+FOLR2-ABCG+ Lipid Associated Macrophages", n_genes=9, ncols=3)
plt.show()

# From inf
root_idx = np.where(adata.obs["initial_states"] == "CD14+-IL1B+SELL+CD16+ Migrating Inflammatory Monocyte-derived Macrophages")[0][0]
adata.uns["iroot"] = root_idx
sc.tl.dpt(adata)

scv.pl.scatter(
    adata,
    color=["final.pop.idents", root_idx, "latent_time", "dpt_pseudotime"],
    fontsize=16,
    cmap="viridis",
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    title=["clusters", "root cell", "latent time", "dpt pseudotime"],save="CellRank/res pseudotime.pdf"
)


cr.pl.heatmap(
    adata,
    model,
    genes=adata.varm['terminal_lineage_drivers']["CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages_corr"].sort_values(ascending=False).index[:100],
    show_absorption_probabilities=True,
    lineages="CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages",
    n_jobs=1,
    backend="loky",
)
plt.show()
plt.savefig("CellRank/res foam gene expression trends heatmap.pdf", bbox_inches='tight')

# PLot some random walks
vk.plot_random_walks(
    25,
    start_ixs={"final.pop.idents": "CD14+-IL1B+SELL+CD16+ Migrating Inflammatory Monocyte-derived Macrophages"},
    max_iter=50,
    successive_hits=5,
    show_progress_bar=False,
    cmap="viridis",linealpha=0.1,
    seed=420,
    ixs_legend_loc="best"
)
plt.show()

vk.plot_random_walks(
    25,
    start_ixs={"final.pop.idents": "CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages"},
    max_iter=50,
    successive_hits=5,
    show_progress_bar=False,
    cmap="viridis",linealpha=0.1,
    seed=420,
    ixs_legend_loc="best"
)
plt.show()

vk.plot_random_walks(
    25,
    start_ixs={"final.pop.idents": "CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages"},
    max_iter=50,
    successive_hits=5,
    show_progress_bar=False,
    cmap="viridis",linealpha=0.1,
    seed=420,
    ixs_legend_loc="best"
)
plt.show()


#=====================================================================================
#=====================================================================================
## Load data mac + mono
adata = sc.read_h5ad("Seurat_Objects/mac.velocyto.refined.h5ad")
adata.obs['final.pop.idents'] = adata.obs['final.pop.idents'].astype('category')
scv.pl.proportions(adata, groupby='final.pop.idents',save="CellRank/macs_only/Splice Ratio.pdf")
adata

## Set vars
true_label=adata.obs.loc[:,"final.pop.idents"]
random_seed=666


# Plot UMAP to check
sc.pl.umap(adata, color="final.pop.idents", palette=pop_cols, return_fig=True)
#plt.show()
plt.savefig("CellRank/macs_only/UMAP.pdf", bbox_inches='tight')

#Pre-process
scv.pp.filter_and_normalize(adata, n_top_genes=2000, subset_highly_variable=False, min_shared_counts=20)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=50)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

# Run scvelo
scv.tl.recover_dynamics(adata, n_jobs=5)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

# Visualize
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_fontsize=12, title="", density=3, alpha= 1, size=3, smooth=0.8,n_neighbors=200, min_mass=1, show=False, color='final.pop.idents', palette=pop_cols, legend_loc='right margin')
plt.savefig("CellRank/macs_only/velocyto.png", bbox_inches='tight')

# Find terminal states
cr.tl.terminal_states(adata, weight_connectivities=0.2, cluster_key='final.pop.idents')
cr.pl.terminal_states(adata,save="CellRank/macs_only/terminal states.pdf")
cr.pl.terminal_states(adata, discrete = True, save="CellRank/macs_only/auto terminal states discrete.pdf")
cr.pl.terminal_states(adata)

# Find initial states
cr.tl.initial_states(adata, cluster_key="final.pop.idents")
cr.pl.initial_states(adata, discrete=True, save="CellRank/macs_only/initial states discrete.pdf")
cr.pl.initial_states(adata, discrete=False, save="CellRank/macs_only/initial states.pdf")

# Compute lineages and plot fate maps
cr.tl.lineages(adata)
cr.pl.lineages(adata, same_plot=False)


#==========================================================================================
## Run customized cellrank in order to set initial and terminal states properley
#==========================================================================================
## Import kernels and compute the trnasition matrixes
from cellrank.tl.kernels import VelocityKernel
from cellrank.tl.kernels import ConnectivityKernel

vk = VelocityKernel(adata).compute_transition_matrix()
ck = ConnectivityKernel(adata).compute_transition_matrix()

## combine kernels
combined_kernel = 0.8 * vk + 0.2 * ck

## Add the GPCCA estimator
from cellrank.tl.estimators import GPCCA

g = GPCCA(combined_kernel)

## Compute matrix decomposition
g.compute_schur(n_components=20)
g.plot_spectrum()
plt.show()

# Compute macrostates
g.compute_macrostates(cluster_key="final.pop.idents", n_states=4)
g.plot_macrostates(discrete=True)
plt.savefig("CellRank/macs_only/macrostates discrete.pdf", bbox_inches='tight')
g.plot_macrostates()
plt.savefig("CellRank/macs_only/macrostates.pdf", bbox_inches='tight')

# Set terminal states
g.set_terminal_states_from_macrostates(names = ["CD14+TREM2-OLR1+ABCA+ Foamy Macrophages"], cluster_key="final.pop.idents")
g.plot_terminal_states(discrete=True)
plt.savefig("CellRank/macs_only/terminal states discrete.pdf", bbox_inches='tight')
g.plot_terminal_states()
plt.savefig("CellRank/macs_only/terminal states.pdf", bbox_inches='tight')

# Set initial states
g._set_initial_states_from_macrostates(names = ["CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages", "CD14+-IL1B+SELL+CD16+ Migrating Inflammatory Monocyte-derived Macrophages"])

# Compute fate probalities
g.compute_absorption_probabilities()

# Compute lineage drivers
g.compute_lineage_drivers()

g.plot_lineage_drivers(lineage="CD14+TREM2-OLR1+ABCA+ Foamy Macrophages", n_genes=9, ncols=3)
plt.show()

g.plot_lineage_drivers(lineage="CD14+CD16+CD64-SLAN+ Non-Classical Monocytes", n_genes=9, ncols=3)
plt.show()

## Now visualize using 'normal' cr methods
# Compute fate maps
cr.tl.lineages(adata)
cr.pl.lineages(adata, same_plot=False)
plt.savefig("CellRank/macs_only/lineages separate.pdf", bbox_inches='tight')
cr.pl.lineages(adata, same_plot=True)
plt.savefig("CellRank/macs_only/lineages combined.pdf", bbox_inches='tight')

# Compute Latent time
scv.tl.recover_latent_time(adata, root_key="initial_states_probs", end_key="terminal_states_probs")
scv.tl.paga(
    adata,
    groups="final.pop.idents",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="velocity_pseudotime")
    
cr.pl.cluster_fates(
    adata,
    mode="paga_pie",
    cluster_key="final.pop.idents",
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    legend_loc="top left out",
    node_size_scale=5,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA"
)
plt.show()
plt.savefig("CellRank/macs_only/cell fate directed PAGA.pdf", bbox_inches='tight')

cr.pl.cluster_fates(
    adata,
    mode="violin",
    cluster_key="final.pop.idents",
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    legend_loc="top left out",
    node_size_scale=5,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA"
)
plt.show()
plt.savefig("CellRank/macs_only/cell fate violin.pdf", bbox_inches='tight')

cr.pl.cluster_fates(
    adata,
    mode="heatmap",
    cluster_key="final.pop.idents",
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    title="directed PAGA"
)
plt.show()
plt.savefig("CellRank/macs_only/cell fate heatmap.pdf", bbox_inches='tight')

cr.pl.cluster_fates(adata, mode="bar", cluster_key='final.pop.idents', ncols = 3, figsize=(15,25))
plt.show()
plt.savefig("CellRank/macs_only/cell fate barplot.pdf", bbox_inches='tight')

# Find lineage drivers
cr.tl.lineage_drivers(adata)

# Plot top 5 per lineage
cr.pl.lineage_drivers(adata, lineage="CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages", n_genes=9, ncols=3)
plt.show()
plt.savefig("CellRank/macs_only/res lineage drivers.pdf", bbox_inches='tight')



cr.pl.lineage_drivers(adata, lineage="CD14+CD16+CD64-SLAN+ Non-Classical Monocytes", n_genes=9, ncols=3)
plt.show()
plt.savefig("CellRank/macs_only/non clas lineage drivers.pdf", bbox_inches='tight')

cr.pl.lineage_drivers(adata, lineage="CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages", n_genes=9, ncols=3)
plt.show()
plt.savefig("CellRank/macs_only/foam lineage drivers.pdf", bbox_inches='tight')

## Gene expression trends
# From Clas mons
root_idx = np.where(adata.obs["initial_states"] == "CD14+-IL1B+SELL+CD16+ Migrating Inflammatory Monocyte-derived Macrophages")[0][0]
adata.uns["iroot"] = root_idx
sc.tl.dpt(adata)

scv.pl.scatter(
    adata,
    color=["final.pop.idents", root_idx, "latent_time", "dpt_pseudotime"],
    fontsize=16,
    cmap="viridis",
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    title=["clusters", "root cell", "latent time", "dpt pseudotime"],save="CellRank/macs_only/clas pseudotime.pdf"
)

model = cr.ul.models.GAM(adata)
cr.pl.gene_trends(
    adata,
    model=model,
    data_key="X",
    genes=["PLIN2", "TIMP1", "ABCA1", "FOLR2", "CD9", "PLTP", "FCGR3A", "CDKN1C", "IFITM1"],
    ncols=3,
    time_key="latent_time",
    legend_loc="none",
    same_plot=True,
    hide_cells=True,
    figsize=(15, 15),
    n_test_points=200,
)
plt.show()
plt.savefig("CellRank/macs_only/clas gene expression trends some genes latent time.pdf", bbox_inches='tight')

cr.pl.gene_trends(
    adata,
    model=model,
    data_key="X",
    genes=["PLIN2", "TIMP1", "ABCA1", "FOLR2", "CD9", "PLTP", "FCGR3A", "CDKN1C", "IFITM1"],
    ncols=3,
    time_key="dpt_pseudotime",
    legend_loc="none",
    same_plot=True,
    hide_cells=True,
    figsize=(15, 15),
    n_test_points=200,
)
plt.show()
plt.savefig("CellRank/macs_only/clas gene expression trends some genes dpt.pdf", bbox_inches='tight')

cr.pl.heatmap(
    adata,
    model,
    genes=adata.varm['terminal_lineage_drivers']["CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages_corr"].sort_values(ascending=False).index[:100],
    show_absorption_probabilities=True,
    lineages="CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages",
    n_jobs=1,
    backend="loky",
)
plt.show()
plt.savefig("CellRank/macs_only/clas foam gene expression trends heatmap.pdf", bbox_inches='tight')

cr.pl.heatmap(
    adata,
    model,
    genes=adata.varm['terminal_lineage_drivers']["CD14+CD16+CD64-SLAN+ Non-Classical Monocytes_corr"].sort_values(ascending=False).index[:100],
    show_absorption_probabilities=True,
    lineages="CD14+CD16+CD64-SLAN+ Non-Classical Monocytes",
    n_jobs=1,
    backend="loky",
)
plt.show()
plt.savefig("CellRank/macs_only/clas non clas gene expression trends heatmap.pdf", bbox_inches='tight')

# From res
root_idx = np.where(adata.obs["initial_states"] == "CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages")[0][0]
adata.uns["iroot"] = root_idx
sc.tl.dpt(adata)

scv.pl.scatter(
    adata,
    color=["final.pop.idents", root_idx, "latent_time", "dpt_pseudotime"],
    fontsize=16,
    cmap="viridis",
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    title=["clusters", "root cell", "latent time", "dpt pseudotime"],save="CellRank/macs_only/res pseudotime.pdf"
)

model = cr.ul.models.GAM(adata)
cr.pl.gene_trends(
    adata,
    model=model,
    data_key="X",
    genes=["PLIN2", "TIMP1", "ABCA1", "FOLR2", "CD9", "PLTP", "FCGR3A", "CDKN1C", "IFITM1"],
    ncols=3,
    time_key="dpt_pseudotime",
    legend_loc="none",
    same_plot=True,
    hide_cells=True,
    figsize=(15, 15),
    n_test_points=200,
)
plt.show()
plt.savefig("CellRank/macs_only/res gene expression trends some genes dpt.pdf", bbox_inches='tight')

cr.pl.heatmap(
    adata,
    model,
    genes=adata.varm['terminal_lineage_drivers']["CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages_corr"].sort_values(ascending=False).index[:100],
    show_absorption_probabilities=True,
    lineages="CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages",
    n_jobs=1,
    backend="loky",
)
plt.show()
plt.savefig("CellRank/macs_only/res foam gene expression trends heatmap.pdf", bbox_inches='tight')

# PLot some random walks
vk.plot_random_walks(
    25,
    start_ixs={"final.pop.idents": "CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages"},
    max_iter=50,
    successive_hits=5,
    show_progress_bar=False,
    cmap="viridis",linealpha=0.1,
    seed=420,
    ixs_legend_loc="best"
)
plt.show()

vk.plot_random_walks(
    25,
    start_ixs={"final.pop.idents": "CD14+-IL1B+SELL+CD16+ Migrating Inflammatory Monocyte-derived Macrophages"},
    max_iter=50,
    successive_hits=10,
    show_progress_bar=False,
    cmap="viridis",linealpha=0.1,
    seed=420,
    ixs_legend_loc="best"
)
plt.show()


##===================================================================================================
## Run on archetypes
##===================================================================================================
## Load data mac + mono
adata = sc.read_h5ad("Seurat_Objects/myeloid.velocyto.refined.h5ad")
adata.obs['archetypes'] = adata.obs['archetypes'].astype('category')
scv.pl.proportions(adata, groupby='archetypes',save="CellRank/Archetypes/Splice Ratio.pdf")
adata

## Set vars
true_label=adata.obs.loc[:,"archetypes"]
random_seed=666


# Plot UMAP to check
sc.pl.umap(adata, color="archetypes", palette=arch_cols, return_fig=True)
#plt.show()
plt.savefig("CellRank/Archetypes/UMAP.pdf", bbox_inches='tight')

#Pre-process
scv.pp.filter_and_normalize(adata, n_top_genes=2000, subset_highly_variable=False, min_shared_counts=20)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=50)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

# Run scvelo
scv.tl.recover_dynamics(adata, n_jobs=5)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

# Visualize
scv.pl.velocity_embedding_stream(adata, basis="umap", legend_fontsize=12, title="", density=3, alpha= 1, size=3, smooth=0.8,n_neighbors=200, min_mass=1, show=False, color='archetypes', palette=arch_cols, legend_loc='right margin')
plt.savefig("CellRank/Archetypes/velocyto.png", bbox_inches='tight')

# Find terminal states
cr.tl.terminal_states(adata, weight_connectivities=0.2, cluster_key='archetypes')
cr.pl.terminal_states(adata,save="CellRank/Archetypes/terminal states.pdf")
cr.pl.terminal_states(adata, discrete = True, save="CellRank/Archetypes/auto terminal states discrete.pdf")
cr.pl.terminal_states(adata)

# Find initial states
cr.tl.initial_states(adata, cluster_key="archetypes")
cr.pl.initial_states(adata, discrete=True, save="CellRank/Archetypes/initial states discrete.pdf")
cr.pl.initial_states(adata, discrete=False, save="CellRank/Archetypes/initial states.pdf")

# Compute lineages and plot fate maps
cr.tl.lineages(adata)
cr.pl.lineages(adata, same_plot=False)


#==========================================================================================
## Run customized cellrank in order to set initial and terminal states properley
#==========================================================================================
## Import kernels and compute the trnasition matrixes
from cellrank.tl.kernels import VelocityKernel
from cellrank.tl.kernels import ConnectivityKernel

vk = VelocityKernel(adata).compute_transition_matrix()
ck = ConnectivityKernel(adata).compute_transition_matrix()

## combine kernels
combined_kernel = 0.8 * vk + 0.2 * ck

## Add the GPCCA estimator
from cellrank.tl.estimators import GPCCA

g = GPCCA(combined_kernel)

## Compute matrix decomposition
g.compute_schur(n_components=20)
g.plot_spectrum()
plt.show()

# Compute macrostates
g.compute_macrostates(cluster_key="archetypes", n_states=4)
g.plot_macrostates(discrete=True)
g.plot_macrostates()

# Set terminal states
g.set_terminal_states_from_macrostates(names = ["Foamy", "Monocytes_1"], cluster_key="archetypes")
g.plot_terminal_states(discrete=True)
g.plot_terminal_states()

# Set initial states
g._set_initial_states_from_macrostates(names = ["Monocytes_2", "Resident"])

# Compute fate probalities
g.compute_absorption_probabilities()

# Compute lineage drivers
g.compute_lineage_drivers()

g.plot_lineage_drivers(lineage="Foamy", n_genes=9, ncols=3)
plt.show()

g.plot_lineage_drivers(lineage="Monocytes_1", n_genes=9, ncols=3)
plt.show()

## Now visualize using 'normal' cr methods
# Compute fate maps
cr.tl.lineages(adata)
cr.pl.lineages(adata, same_plot=False)
cr.pl.lineages(adata, same_plot=True)

# Compute Latent time
scv.tl.recover_latent_time(adata, root_key="initial_states_probs", end_key="terminal_states_probs")
scv.tl.paga(
    adata,
    groups="archetypes",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="velocity_pseudotime")
    
cr.pl.cluster_fates(
    adata,
    mode="paga_pie",
    cluster_key="archetypes",
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    legend_loc="top left out",
    node_size_scale=5,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA"
)
plt.show()
plt.savefig("CellRank/Archetypes/cell fate directed PAGA.pdf", bbox_inches='tight')

cr.pl.cluster_fates(
    adata,
    mode="violin",
    cluster_key="archetypes",
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    legend_loc="top left out",
    node_size_scale=5,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA"
)
plt.show()
plt.savefig("CellRank/Archetypes/cell fate violin.pdf", bbox_inches='tight')

cr.pl.cluster_fates(
    adata,
    mode="heatmap",
    cluster_key="archetypes",
    basis="umap",
    legend_kwargs={"loc": "top right out"},
    title="directed PAGA"
)
plt.show()
plt.savefig("CellRank/Archetypes/cell fate heatmap.pdf", bbox_inches='tight')

cr.pl.cluster_fates(adata, mode="bar", cluster_key='archetypes', ncols = 3, figsize=(15,25))
plt.show()
plt.savefig("CellRank/Archetypes/cell fate barplot.pdf", bbox_inches='tight')

# Find lineage drivers
cr.tl.lineage_drivers(adata)

# Plot top 5 per lineage
cr.pl.lineage_drivers(adata, lineage="Foamy", n_genes=9, ncols=3)
plt.show()
plt.savefig("CellRank/Archetypes/res lineage drivers.pdf", bbox_inches='tight')

## Gene expression trends
# From Clas mons
root_idx = np.where(adata.obs["initial_states"] == "Monocytes_2")[0][0]
adata.uns["iroot"] = root_idx
sc.tl.dpt(adata)

scv.pl.scatter(
    adata,
    color=["archetypes", root_idx, "latent_time", "dpt_pseudotime"],
    fontsize=16,
    cmap="viridis",
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    title=["clusters", "root cell", "latent time", "dpt pseudotime"],save="CellRank/Archetypes/clas pseudotime.pdf"
)

model = cr.ul.models.GAM(adata)
cr.pl.gene_trends(
    adata,
    model=model,
    data_key="X",
    genes=["PLIN2", "TIMP1", "ABCA1", "FOLR2", "CD9", "PLTP", "FCGR3A", "CDKN1C", "IFITM1"],
    ncols=3,
    time_key="latent_time",
    legend_loc="none",
    same_plot=True,
    hide_cells=True,
    figsize=(15, 15),
    n_test_points=200,
)
plt.show()
plt.savefig("CellRank/Archetypes/clas gene expression trends some genes latent time.pdf", bbox_inches='tight')

cr.pl.gene_trends(
    adata,
    model=model,
    data_key="X",
    genes=["PLIN2", "TIMP1", "ABCA1", "FOLR2", "CD9", "PLTP", "FCGR3A", "CDKN1C", "IFITM1"],
    ncols=3,
    time_key="dpt_pseudotime",
    legend_loc="none",
    same_plot=True,
    hide_cells=True,
    figsize=(15, 15),
    n_test_points=200,
)
plt.show()
plt.savefig("CellRank/Archetypes/clas gene expression trends some genes dpt.pdf", bbox_inches='tight')

cr.pl.heatmap(
    adata,
    model,
    genes=adata.varm['terminal_lineage_drivers']["CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages_corr"].sort_values(ascending=False).index[:100],
    show_absorption_probabilities=True,
    lineages="CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages",
    n_jobs=1,
    backend="loky",
)
plt.show()
plt.savefig("CellRank/Archetypes/clas foam gene expression trends heatmap.pdf", bbox_inches='tight')

cr.pl.heatmap(
    adata,
    model,
    genes=adata.varm['terminal_lineage_drivers']["CD14+CD16+CD64-SLAN+ Non-Classical Monocytes_corr"].sort_values(ascending=False).index[:100],
    show_absorption_probabilities=True,
    lineages="CD14+CD16+CD64-SLAN+ Non-Classical Monocytes",
    n_jobs=1,
    backend="loky",
)
plt.show()
plt.savefig("CellRank/Archetypes/clas non clas gene expression trends heatmap.pdf", bbox_inches='tight')

# From res
root_idx = np.where(adata.obs["initial_states"] == "Resident")[0][0]
adata.uns["iroot"] = root_idx
sc.tl.dpt(adata)

scv.pl.scatter(
    adata,
    color=["archetypes", root_idx, "latent_time", "dpt_pseudotime"],
    fontsize=16,
    cmap="viridis",
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    title=["clusters", "root cell", "latent time", "dpt pseudotime"],save="CellRank/Archetypes/res pseudotime.pdf"
)

model = cr.ul.models.GAM(adata)
cr.pl.gene_trends(
    adata,
    model=model,
    data_key="X",
    genes=["PLIN2", "TIMP1", "ABCA1", "FOLR2", "CD9", "PLTP", "FCGR3A", "CDKN1C", "IFITM1"],
    ncols=3,
    time_key="dpt_pseudotime",
    legend_loc="none",
    same_plot=True,
    hide_cells=True,
    figsize=(15, 15),
    n_test_points=200,
)
plt.show()
plt.savefig("CellRank/Archetypes/res gene expression trends some genes dpt.pdf", bbox_inches='tight')

cr.pl.heatmap(
    adata,
    model,
    genes=adata.varm['terminal_lineage_drivers']["CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages_corr"].sort_values(ascending=False).index[:100],
    show_absorption_probabilities=True,
    lineages="CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages",
    n_jobs=1,
    backend="loky",
)
plt.show()
plt.savefig("CellRank/Archetypes/res foam gene expression trends heatmap.pdf", bbox_inches='tight')

# PLot some random walks
vk.plot_random_walks(
    50,
    start_ixs={"archetypes": "Monocytes"},
    max_iter=100,
    successive_hits=5,
    show_progress_bar=False,
    cmap="viridis",linealpha=0.1,
    seed=420,
    ixs_legend_loc="best"
)
plt.show()
plt.savefig("CellRank/Archetypes/random walks res start.pdf", bbox_inches='tight')

vk.plot_random_walks(
    50,
    start_ixs={"archetypes": "CD14+CD16-CD64+SLAN- Classical Monocytes"},
    max_iter=100,
    successive_hits=5,
    show_progress_bar=False,
    cmap="viridis",linealpha=0.1,
    seed=420,
    ixs_legend_loc="best"
)
plt.show()
plt.savefig("CellRank/Archetypes/random walks clas mon start.pdf", bbox_inches='tight')

vk.plot_random_walks(
    50,
    start_ixs={"archetypes": "CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages"},
    stop_ixs={"archetypes": "CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages"},
    max_iter=100,
    successive_hits=5,
    show_progress_bar=False,
    cmap="viridis",linealpha=0.1,
    seed=420,
    ixs_legend_loc="best"
)
plt.show()
plt.savefig("CellRank/Archetypes/random walks res start foam stop.pdf", bbox_inches='tight')

vk.plot_random_walks(
    50,
    start_ixs={"archetypes": "CD14+CD16-CD64+SLAN- Classical Monocytes"},
    stop_ixs={"archetypes": "CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages"},
     max_iter=200,
    successive_hits=5,
    show_progress_bar=False,
    cmap="viridis",linealpha=0.1,
    seed=420,
    ixs_legend_loc="best"
)
plt.show()
plt.savefig("CellRank/Archetypes/random walks clas mon start foam stop.pdf", bbox_inches='tight')


