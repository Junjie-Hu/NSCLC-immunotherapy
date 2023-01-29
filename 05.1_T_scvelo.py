import anndata
import scanpy
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# read loom file
sample_all = anndata.read_loom("sample_all.loom")
#sample_all = scanpy.read_loom("sample_all.loom")
# read seurat information
sample_obs = pd.read_csv("T8_cellID_obs.csv")
umap = pd.read_csv("T8_cell_embeddings.csv",index_col=0)
cell_clusters = pd.read_csv("T8_clusters.csv")
custer_color = pd.read_csv("T8_color.csv")

# extracted Cell IDs from Seurat
sample_all.obs.index = sample_all.obs.obs_names
sample_use = sample_all[sample_obs["x"]]

# view cell ID order
sample_use.obs.index

# add the UMAP coordinates
sample_use.obsm['X_umap'] = umap.values

# add cluster info 
sample_use.obs['cell_clusters']=cell_clusters.values

## run scvelo
# Basic preprocessing
scv.pp.filter_and_normalize(sample_use)
scv.pp.moments(sample_use)
# Velocity Tools
scv.tl.velocity(sample_use, mode = "deterministic")
#scv.tl.velocity(sample_use, mode = "stochastic")
scv.tl.velocity_graph(sample_use)
# Visualization
plot_color = list(custer_color["color"])
scv.pl.velocity_embedding_stream(sample_use, basis='umap',color="cell_clusters",palette=plot_color,
legend_loc="none",title="",figsize=(6,6),dpi=300,save="T8_scvelo_steady2.png")




