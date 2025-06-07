library(tidyverse)
library(Seurat)
library(ggsci)
library(RColorBrewer)

merged_regiond <- readRDS('~/merged_dim30_harmony.rds')
table(merged_regiond$cluster)

sub_merged_regiond <- subset(merged_regiond,
                             subset = cluster %in% c('Co_D8_c3','Co_D20_c7','Co_D20_c12','Co_D8_c9','Co_D20_c17','Co_D8_c0','Co_D8_c2','Co_D8_c8',
                                                     'Co_D8_c11','Co_D8_c14','Co_D8_c16','Co_D20_c2','Co_D20_c5','Co_D8_c5','Co_D8_c7','Co_D8_c4',
                                                     'Co_D20_c3','Co_D20_c4','Co_D20_c8','Co_D20_c9','Co_D20_c14','Co_D20_c13','Co_D8_c12',
                                                     'Co_D20_c19','Co_D8_c15','Co_D20_c21','Co_D8_c6','Co_D8_c1','Co_D20_c0','Co_D20_c1','Co_D20_c11',
                                                     'Co_D20_c16','Co_D8_c10','Co_D20_c6','Co_D20_c18','Co_D20_c15','Co_D20_c10'))

sub_merged_regiond$cluster <- factor(sub_merged_regiond$cluster,
                                     levels = c('Co_D8_c3','Co_D20_c7','Co_D20_c12','Co_D8_c9','Co_D20_c17','Co_D8_c0','Co_D8_c2','Co_D8_c8',
                                                'Co_D8_c11','Co_D8_c14','Co_D8_c16','Co_D20_c2','Co_D20_c5','Co_D8_c5','Co_D8_c7','Co_D8_c4',
                                                'Co_D20_c3','Co_D20_c4','Co_D20_c8','Co_D20_c9','Co_D20_c14','Co_D20_c13','Co_D8_c12',
                                                'Co_D20_c19','Co_D8_c15','Co_D20_c21','Co_D8_c6','Co_D8_c1','Co_D20_c0','Co_D20_c1','Co_D20_c11',
                                                'Co_D20_c16','Co_D8_c10','Co_D20_c6','Co_D20_c18','Co_D20_c15','Co_D20_c10'))

sub_merged_regiond$cell_type <- ''
sub_merged_regiond$cell_type[sub_merged_regiond$cluster %in% c('Co_D8_c3','Co_D20_c7','Co_D20_c12')] <- 'Endo'
sub_merged_regiond$cell_type[sub_merged_regiond$cluster %in% c('Co_D8_c9','Co_D20_c17','Co_D8_c0','Co_D8_c2','Co_D8_c8','Co_D8_c11','Co_D8_c14','Co_D8_c16','Co_D20_c2','Co_D20_c5')] <- 'Epi.mes.'
sub_merged_regiond$cell_type[sub_merged_regiond$cluster %in% c('Co_D8_c5','Co_D8_c7','Co_D8_c4','Co_D20_c3','Co_D20_c4','Co_D20_c8','Co_D20_c9','Co_D20_c14','Co_D20_c13','Co_D8_c12',
                                                               'Co_D20_c19','Co_D8_c15')] <- 'Meso.'
sub_merged_regiond$cell_type[sub_merged_regiond$cluster %in% c('Co_D20_c21')] <- 'Ep.'
sub_merged_regiond$cell_type[sub_merged_regiond$cluster %in% c('Co_D8_c6','Co_D8_c1','Co_D20_c0','Co_D20_c1','Co_D20_c11','Co_D20_c16')] <- 'CM'
sub_merged_regiond$cell_type[sub_merged_regiond$cluster %in% c('Co_D8_c10','Co_D20_c6','Co_D20_c18','Co_D20_c15','Co_D20_c10')] <- 'EC'

sub_merged_regiond$cell_type <- factor(sub_merged_regiond$cell_type,
                                       levels = rev(c('EC', 'Ep.', 'Meso.','CM',  'Epi.mes.', 'Endo')))
CM_regiond <- subset(sub_merged_regiond, subset = cell_type %in% c('CM'))
table(CM_regiond$group)
CM_regiond_D20 <- subset(CM_regiond, subset = group %in% c('Co_D20'))

D20_CO <- readRDS("~/sub_reduction_dims30_res1.rds")
CM_D20_CO <- subset(D20_CO, subset = seurat_clusters %in% c('0','1','2','3','4','6','8'))

D20_VCO <- readRDS("~/sub_reduction_dims30_res1.rds")
CM_D20_VCO <- subset(D20_VCO, subset = seurat_clusters %in% c('0','1','3','4','5','16','18'))

CM_regiond_D20[['group']] <- 'D20_regiond'
CM_D20_CO[['group']] <- 'D20_CO'
CM_D20_VCO[['group']] <- 'D20_VCO'

merged_CM <- merge(CM_regiond_D20, y=c(CM_D20_CO, CM_D20_VCO))
table(merged_CM$group)
Idents(merged_CM) <- merged_CM$group

merged_CM <- NormalizeData(merged_CM, normalization.method = "LogNormalize", scale.factor = 10000)

library(scRNAtoolVis)
library(ggplot2)

annoGene <- c('APOE', 'IGFBP7', 'IGF1R','ERBB4','HEY2', 'TBX20', 'NKX2-5', 'GJA5','MLLT3', 'FZD6', 'MAGI2',
              'BMP7', 'SMAD6', 'BMP2', 'FBN2', 'BAMBI',
              'SLC2A3', 'SLC5A3', 'SLC7A2','COL3A1', 'ADAMTS3', 'RBPJ','RPS6')

AverageHeatmap(object = merged_CM,
               #cluster_columns = FALSE,
               column_names_rot = 0,
               cluster.order = c('D20_regiond','D20_CO',"D20_VCO"),
               markerGene = markers$gene,
               showRowNames = F,
               clusterAnnoName = F,
               # htCol = c('#F2B77C',"white","#077E97"),
               htCol = c('#A0A0A4',"white","#AD07E3"),
               markGenes = annoGene,
               fontsize = 9)



