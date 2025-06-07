library(tidyverse)
library(Seurat)
library(ggsci)
library(RColorBrewer)

merged <- readRDS('/home/lvxq/230107A_23908_kdl_yz/result_new/clustering/D8_D20/merged_dim30_harmony.rds')
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
										
library(Nebulosa)
plot_density(sub_merged, reduction = 'umap', c('WT1'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('TBX18'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('SFRP2'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('ITLN1'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('WT1', 'TBX18', 'SFRP2', 'ITLN1'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))

plot_density(sub_merged, reduction = 'umap', c('TNNT2'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('NKX2-5'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('TBX20'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('IRX4'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('TNNT2','NKX2-5', 'TBX20', 'IRX4'),joint = TRUE) +
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))

plot_density(sub_merged, reduction = 'umap', c('ALB'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('AFP'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('HNF4A'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('ALDH1A1'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('ALB', 'AFP', 'HNF4A', 'ALDH1A1'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))

plot_density(sub_merged, reduction = 'umap', c('FN1'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('PDGFRA'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('KRT19'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('CDH1'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('FN1', 'PDGFRA', 'KRT19', 'CDH1'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))

plot_density(sub_merged, reduction = 'umap', c('CDH5'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('KDR'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('TIE1'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('FLI1'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('CD34'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('CDH5', 'KDR', 'TIE1','FLI1','CD34'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))

plot_density(sub_merged, reduction = 'umap', c('COL1A1'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('DCN'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('LUM'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('VIM'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))
plot_density(sub_merged, reduction = 'umap', c('COL1A1', 'DCN', 'LUM', 'VIM'),joint = TRUE)+
  scale_color_gradientn(colors = colorRampPalette(c(brewer.pal(9, "Blues"))[2:8])(100))

										