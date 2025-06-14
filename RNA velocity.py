import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

sample_1 = anndata.read_loom("~/loom/D8-CO_count_result.loom")
sample_1.var_names_make_unique()
sample_2 = anndata.read_loom("~/loom/D20-CO_count_result.loom")
sample_2.var_names_make_unique()

sample_obs = pd.read_csv("~/cellID_obs.csv")
umap = pd.read_csv("~/cell_embeddings.csv")
cell_clusters = pd.read_csv("~/clusters_obs.csv")

sample_1 = sample_1[np.isin(sample_1.obs.index,sample_obs["x"])]
sample_2 = sample_2[np.isin(sample_2.obs.index,sample_obs["x"])]

adata = sample_1.concatenate(sample_2)
adata_index = pd.DataFrame(adata.obs.index)
adata_index = adata_index.rename(columns = {0:'Cell ID'})
adata_index = adata_index.rename(columns = {"CellID":'Cell ID'})
adata_index['Cell ID'] = adata_index['Cell ID'].str.replace('-', '_', n=1)
rep=lambda x : x.split("-")[0]
adata_index["Cell ID"]=adata_index["Cell ID"].apply(rep)

umap = umap.rename(columns = {'Unnamed: 0':'Cell ID'})
umap['Cell ID'] = umap['Cell ID'].str.replace('-', '_', n=1)
umap = umap[np.isin(umap["Cell ID"],adata_index["Cell ID"])]
umap=umap.drop_duplicates(subset=["Cell ID"])
umap_ordered = adata_index.merge(umap, on = "Cell ID")
umap_ordered = umap_ordered.iloc[:,1:]
adata.obsm['X_umap'] = umap_ordered.values

cell_clusters = cell_clusters.rename(columns = {'Unnamed: 0':'Cell ID'})
cell_clusters['Cell ID'] = cell_clusters['Cell ID'].str.replace('-', '_', n=1)
cell_clusters = cell_clusters[np.isin(cell_clusters["Cell ID"],adata_index["Cell ID"])]
cell_clusters=cell_clusters.drop_duplicates(subset=["Cell ID"])
cell_clusters_ordered = adata_index.merge(cell_clusters, on = "Cell ID")
cell_clusters_ordered = cell_clusters_ordered.iloc[:,1:]
adata.obs['clusters']=cell_clusters_ordered.values

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode = "stochastic")
scv.tl.velocity_graph(adata)
scv.pl.proportions(adata)
scv.pl.velocity_embedding(adata, basis = 'umap')
scv.pl.velocity_embedding_stream(adata, basis = 'umap')








