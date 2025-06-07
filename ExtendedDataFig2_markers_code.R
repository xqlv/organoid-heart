library(tidyverse)
library(Seurat)
library(ggsci)

merged <- readRDS('~/merged_dim30_harmony.rds')
table(merged$cluster)

sub_merged <- subset(merged,
                     subset = cluster %in% c('Co_D8_c3','Co_D20_c7','Co_D20_c12','Co_D8_c9','Co_D20_c17','Co_D8_c0','Co_D8_c2','Co_D8_c8',
                                             'Co_D8_c11','Co_D8_c14','Co_D8_c16','Co_D20_c2','Co_D20_c5','Co_D8_c5','Co_D8_c7','Co_D8_c4',
                                             'Co_D20_c3','Co_D20_c4','Co_D20_c8','Co_D20_c9','Co_D20_c14','Co_D20_c13','Co_D8_c12',
                                             'Co_D20_c19','Co_D8_c15','Co_D20_c21','Co_D8_c6','Co_D8_c1','Co_D20_c0','Co_D20_c1','Co_D20_c11',
                                             'Co_D20_c16','Co_D8_c10','Co_D20_c6','Co_D20_c18','Co_D20_c15','Co_D20_c10'))

sub_merged$cluster <- factor(sub_merged$cluster,
                             levels = c('Co_D8_c3','Co_D20_c7','Co_D20_c12','Co_D8_c9','Co_D20_c17','Co_D8_c0','Co_D8_c2','Co_D8_c8',
                                        'Co_D8_c11','Co_D8_c14','Co_D8_c16','Co_D20_c2','Co_D20_c5','Co_D8_c5','Co_D8_c7','Co_D8_c4',
                                        'Co_D20_c3','Co_D20_c4','Co_D20_c8','Co_D20_c9','Co_D20_c14','Co_D20_c13','Co_D8_c12',
                                        'Co_D20_c19','Co_D8_c15','Co_D20_c21','Co_D8_c6','Co_D8_c1','Co_D20_c0','Co_D20_c1','Co_D20_c11',
                                        'Co_D20_c16','Co_D8_c10','Co_D20_c6','Co_D20_c18','Co_D20_c15','Co_D20_c10'))

features <-  c('TTR','ALDH1A1','ALB','AFP','HNF4A',
               'KRT19','CDH1','CDH2','EPCAM','PDGFRA','NCAM1',
               'FN1','DCN','COL1A1','COL3A1','COL5A1','LUM','VIM','ACTA2','PDGFRB','RSG5','TAGLN',
               'WT1','TBX18','BNC1','BNC2','SFRP2','ITLN1',
               'GATA4','TBX20','NKX2-5','TBX5','MYL7','TNNT2','TNNI1','NPPA','IRX4','NR2F1','HEY2',
               'CDH5','EPCAM1','KDR','TIE1',
               'TOP2A','MKI67')
			   
avg_expr <- AverageExpression(sub_merged,
                              group.by = 'cluster',
                              features = features,
                              assays = 'RNA',
                              slot = "data") %>% as.data.frame()
colnames(avg_expr) <- gsub('RNA.','',colnames(avg_expr))

library(RColorBrewer)
cols <- c(brewer.pal(11, "RdBu"))
cols <- c(brewer.pal(9, "Purples"))
cols <- c(brewer.pal(11, "PiYG"))
cols <- c(brewer.pal(9, "BuPu"))
cols <- c(brewer.pal(9, "BuGn"))
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
  color = colorRampPalette(cols)(100),
  # color = colorRampPalette(rev(cols))(100),
  border_color = "white")			   
			   
			   

