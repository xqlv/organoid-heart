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


DimPlot(merged_order, reduction = "umap",group.by = 'pre_cell_type',label = T,repel = T)+
  scale_color_manual(values = c('#717070','#1E4C75','#7E2F8E','#0A9A3E','#7C4E23'))
  
prop.table(table(merged_order$pre_cell_type, merged_order$group), margin = 2) * 100
meta <- merged_order@meta.data
table(meta$group)
meta$group <- factor(meta$group, levels = c('D8_CO','D20_CO'))
p1 <- ggplot(data = meta, aes(x = group, fill = pre_cell_type)) +
  geom_bar(position = "fill",
           width = 0.9,
           color = "white") +
  labs(title = "Sample", x = "", y = "Fraction of cells") +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.background = element_blank(),
    axis.title.x = element_text(size = 10),
    axis.line.x = element_line(size = 0.5),
    axis.ticks.x = element_line(size = 0.5, colour = "black"),
    axis.ticks.length = unit(.5, "lines"),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10),
    strip.background = element_rect(color = "white", fill = "white"),
    strip.text.y = element_text(size = 10)
  ) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c('#717070','#1E4C75','#7E2F8E','#0A9A3E','#7C4E23')) 
p1  
  
  