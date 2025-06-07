library(Seurat)
library(ggplot2)
library(tidyverse)

merged <- readRDS('~/merged.rds')

merged$pre_cell_type[merged$seurat_clusters %in% c('3','4','5','6','9','10','13') & merged$group %in% c('D8_CO')] <- 'CM'
merged$pre_cell_type[merged$seurat_clusters %in% c('12') & merged$group %in% c('D8_CO')] <- 'EC'
merged$pre_cell_type[merged$seurat_clusters %in% c('0','7') & merged$group %in% c('D8_CO')] <- 'EP'
merged$pre_cell_type[merged$seurat_clusters %in% c('1','2','8') & merged$group %in% c('D8_CO')] <- 'MESO'
merged$pre_cell_type[merged$seurat_clusters %in% c('11','14') & merged$group %in% c('D8_CO')] <- 'ENDO'

merged$pre_cell_type[merged$seurat_clusters %in% c('0','1','2','3','4','6','8','11') & merged$group %in% c('D20_CO')] <- 'CM'
merged$pre_cell_type[merged$seurat_clusters %in% c('13') & merged$group %in% c('D20_CO')] <- 'EC'
merged$pre_cell_type[merged$seurat_clusters %in% c('7','10') & merged$group %in% c('D20_CO')] <- 'EP'
merged$pre_cell_type[merged$seurat_clusters %in% c('9') & merged$group %in% c('D20_CO')] <- 'MESO'
merged$pre_cell_type[merged$seurat_clusters %in% c('15') & merged$group %in% c('D20_CO')] <- 'ENDO'

merged_order <- merged
merged_order$group <- factor(merged_order$group, levels = c('D8_CO','D20_CO'))
merged_order$pre_cell_type <- factor(merged_order$pre_cell_type, levels = c('CM','MESO','EP','EC','ENDO'))
merged_order$group_seurat_clusters <- factor(merged_order$group_seurat_clusters,
                                             levels = c('D8_CO_3','D8_CO_4','D8_CO_5','D8_CO_6','D8_CO_9','D8_CO_10','D8_CO_13',
                                                        'D20_CO_0','D20_CO_1','D20_CO_2','D20_CO_3','D20_CO_4','D20_CO_6','D20_CO_8','D20_CO_11',
                                                        
                                                        'D8_CO_1','D8_CO_2','D8_CO_8',
                                                        'D20_CO_9',
                                                        
                                                        'D8_CO_0','D8_CO_7',
                                                        'D20_CO_7','D20_CO_10',
                                                        
                                                        'D8_CO_12',
                                                        'D20_CO_13',
                                                        
                                                        'D8_CO_11','D8_CO_14',
                                                        'D20_CO_15'))

features2 <- c('GATA4','TBX20','MYL7','TNNT2','TTN',
               'COL1A1','COL3A1','ACTA2','RGS5','PDGFRB','KCNJ8','PDGFRA',
               'WT1','TBX18','BNC1','BNC2','ALDH1A2','TCF21','DLK1','CXCL12',
               'CDH5','KDR','TIE1','PECAM1','KRT19','CDH1','TTR')

avg_expr <- AverageExpression(merged_order,
                              group.by = 'group_seurat_clusters',
                              features = unique(features2),
                              assays = 'RNA',
                              slot = "data") %>% as.data.frame()
colnames(avg_expr) <- gsub('RNA.','',colnames(avg_expr))

library(RColorBrewer)
cols <- c(brewer.pal(11, "RdBu"))
pheatmap::pheatmap(
  avg_expr,
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



