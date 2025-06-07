library(Seurat)
library(ggplot2)

merged <- readRDS('~/merged.rds')
DimPlot(merged, reduction = "umap",group.by = 'pre_cell_type',label = T,repel = T)
table(merged$pre_cell_type)
merged_CM <- subset(merged, subset = pre_cell_type %in% c('CM'))

features2 <- c('PLK2','RABGAP1L','HEY2','DHRS3','PLN','CKMT2','CGNL1','IRX3','GJA5',
              'CXCL12','IRX2','IRX1')
table(merged_CM$group_seurat_clusters)
merged_CM$group_seurat_clusters <- factor(merged_CM$group_seurat_clusters,
                                          levels = c('D20_CO_6', 'D20_CO_8', 'D20_CO_0', 'D20_CO_4', 'D20_CO_3', 'D20_CO_1', 'D20_CO_2',
                                            'D8_CO_3',  'D8_CO_4',  'D8_CO_5',  'D8_CO_6',  'D8_CO_9'
                                          ))
avg_exp_merged_CM <- AverageExpression(merged_CM,
                                       assays = 'RNA',
                                       features = features2,
                                       return.seurat = FALSE,
                                       group.by = "group_seurat_clusters",
                                       slot = "data") %>% as.data.frame()

library(RColorBrewer)
cols <- c(brewer.pal(11, "RdBu"))
pheatmap::pheatmap(
  avg_exp_merged_CM,
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
