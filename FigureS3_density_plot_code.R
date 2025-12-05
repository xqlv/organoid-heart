library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggsci)

merged_EP <- readRDS('~/merged_EP_recluster.rds')
table(merged_EP$seurat_clusters)
DimPlot(merged_EP, 
        group.by = "seurat_clusters",
        reduction = "umap", # tsne, umap, pca
        label = T) + 
		scale_color_manual(values = c('#6892CC','#A07EB7','#C270AA','#19B285','#1DB4B7','#1FA9DE','#E864A0','#27AD3F','#B19A1A','#7AA82B'))

library(Nebulosa)
plot_density(merged_EP, reduction = 'umap', c('TNNT2','TTN','ACTN2','MYH7'),joint = TRUE)
plot_density(merged_EP, reduction = 'umap', c('RGS5','PDGFRB'),joint = TRUE)
plot_density(merged_EP, reduction = 'umap', c('KCNJ8','CSPG4','ACTA2','TAGLN'),joint = TRUE)
plot_density(merged_EP, reduction = 'umap', c('TWIST1','SNAI2','KRT19'),joint = TRUE)
plot_density(merged_EP, reduction = 'umap', c('COL5A1','ITGAV','COL16A1'),joint = TRUE)



