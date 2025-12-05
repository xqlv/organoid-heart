library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggsci)

merged_EP <- readRDS('~/merged_EP.rds')
DimPlot(merged_EP, 
        group.by = "group_seurat_clusters",
        reduction = "umap", # tsne, umap, pca
        label = T) + scale_color_manual(values = c('#7030A0','#9E5ECE','#47A7C9','#318AA9'))

features <- c('TNNT1','SFRP2','SFRP5',
              'LRP2','CALB2','C3',
              'HAND1', 'MAB21L2', 'HOXB6', 'HOXB5', 'BNC2',
              'HEY1','HEY2','CFB','C1R','C1S','ELN','DPT')
avg_exp_merged_EP <- AverageExpression(merged_EP,
                                       assays = 'RNA',
                                       features = unique(features),
                                       return.seurat = FALSE,
                                       group.by = "group_seurat_clusters",
                                       slot = "data") %>% as.data.frame()

library(RColorBrewer)
cols <- c(brewer.pal(11, "RdBu"))
pheatmap::pheatmap(
  avg_exp_merged_EP,
  scale = "row",
  cellwidth = 10,
  cellheight = 10,
  fontsize = 8,
  cluster_rows = F,
  cluster_cols = F,
  angle_col = c("90"),
  # gaps_col= c(4,6,8),
  # gaps_row = c(4,6,8),
  # cutree_cols = c(4,6,8),
  number_color = "white",
  color = colorRampPalette(rev(cols))(100),
  # color = colorRampPalette(rev(cols))(100),
  border_color = "white")

