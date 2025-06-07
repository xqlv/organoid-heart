library(tidyverse)
library(Seurat)
library(ggsci)

merged <- readRDS('~/merged_dim30_harmony.rds')
table(merged$cluster)

merged$cell_type <- ''
merged$cell_type[merged$cluster %in% c('Co_D8_c3','Co_D20_c7','Co_D20_c12')] <- 'Endo'
merged$cell_type[merged$cluster %in% c('Co_D8_c9','Co_D20_c17','Co_D8_c0','Co_D8_c2','Co_D8_c8','Co_D8_c11','Co_D8_c14','Co_D8_c16','Co_D20_c2','Co_D20_c5')] <- 'Epi.mes.'
merged$cell_type[merged$cluster %in% c('Co_D8_c5','Co_D8_c7','Co_D8_c4','Co_D20_c3','Co_D20_c4','Co_D20_c8','Co_D20_c9','Co_D20_c14','Co_D20_c13','Co_D8_c12',
                                       'Co_D20_c19','Co_D8_c15')] <- 'Meso.'
merged$cell_type[merged$cluster %in% c('Co_D20_c21')] <- 'Ep.'
merged$cell_type[merged$cluster %in% c('Co_D8_c6','Co_D8_c1','Co_D20_c0','Co_D20_c1','Co_D20_c11','Co_D20_c16')] <- 'CM'
merged$cell_type[merged$cluster %in% c('Co_D8_c10','Co_D20_c6','Co_D20_c18','Co_D20_c15','Co_D20_c10')] <- 'EC'
merged$cell_type[merged$cluster %in% c('Co_D8_c13','Co_D20_c20')] <- 'unknown'

table(merged$cluster, merged$cell_type)


merged$cell_type <- factor(merged$cell_type,
                               levels = rev(c('EC', 'Ep.', 'Meso.','CM',  'Epi.mes.', 'Endo','unknown')))
merged$group <- factor(merged$group,
                           levels = c('Co_D8', 'Co_D20'))
						   
DimPlot(merged, group.by = 'cell_type') + 
  scale_color_manual(values = c('black','#855A31','#F3E18A','#777676','#7A3489','#D62A85','#329947'))						   
						   
prop.table(table(merged$cell_type, merged$orig.ident), margin = 2) * 100
meta <- merged@meta.data
table(meta$orig.ident)
meta$orig.ident <- factor(meta$orig.ident, levels = c('Co_D8', 'Co_D20'))
p1 <- ggplot(data = meta, aes(x = orig.ident, fill = cell_type)) +
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
  scale_fill_manual(values = c('black','#855A31','#D4D832','#777676','#7A3489','#D62A85','#329947')) 
p1						   
						   
						   