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

sub_merged$cell_type <- ''
sub_merged$cell_type[sub_merged$cluster %in% c('Co_D8_c3','Co_D20_c7','Co_D20_c12')] <- 'Endo'
sub_merged$cell_type[sub_merged$cluster %in% c('Co_D8_c9','Co_D20_c17','Co_D8_c0','Co_D8_c2','Co_D8_c8','Co_D8_c11','Co_D8_c14','Co_D8_c16','Co_D20_c2','Co_D20_c5')] <- 'Epi.mes.'
sub_merged$cell_type[sub_merged$cluster %in% c('Co_D8_c5','Co_D8_c7','Co_D8_c4','Co_D20_c3','Co_D20_c4','Co_D20_c8','Co_D20_c9','Co_D20_c14','Co_D20_c13','Co_D8_c12',
                                               'Co_D20_c19','Co_D8_c15')] <- 'Meso.'
sub_merged$cell_type[sub_merged$cluster %in% c('Co_D20_c21')] <- 'Ep.'
sub_merged$cell_type[sub_merged$cluster %in% c('Co_D8_c6','Co_D8_c1','Co_D20_c0','Co_D20_c1','Co_D20_c11','Co_D20_c16')] <- 'CM'
sub_merged$cell_type[sub_merged$cluster %in% c('Co_D8_c10','Co_D20_c6','Co_D20_c18','Co_D20_c15','Co_D20_c10')] <- 'EC'

sub_merged$cell_type <- factor(sub_merged$cell_type,
                               levels = rev(c('EC', 'Ep.', 'Meso.','CM',  'Epi.mes.', 'Endo')))
							   
features <-  c('TTR','ALDH1A1','ALB','AFP','HNF4A',
               'KRT19','CDH1','CDH2','EPCAM','PDGFRA','NCAM1',
               'GATA4','TBX20','NKX2-5','TBX5','MYL7','TNNT2','TNNI1','NPPA','IRX4','NR2F1','HEY2',
               'FN1','DCN','COL1A1','COL3A1','COL5A1','LUM','VIM','ACTA2','PDGFRB','RSG5','TAGLN',
               'WT1','TBX18','BNC1','BNC2','SFRP2','ITLN1',
               'CDH5','EPCAM1','KDR','TIE1',
               'TOP2A','MKI67')							   
							   
library(RColorBrewer)
cols <- c(brewer.pal(11, "Blues"))
DotPlot(sub_merged,
        features = rev(features),
        group.by = 'cell_type')+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))+
  scale_color_gradientn(colours= cols)+coord_flip()							   
							   