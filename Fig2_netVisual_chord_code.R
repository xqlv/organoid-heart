library(tidyverse)
library(Seurat)
library(ggsci)
library(CellChat)
library(patchwork)

cellchat <- readRDS("~/cellchat_Co_D8_D20_2celltypes.rds")

netVisual_bubble(cellchat,
                 remove.isolate = FALSE)

pairLR <- data.frame(interaction_name=rev(c('WNT5B_FZD6','WNT5B_FZD4','VTN_ITGAV_ITGB1','TGFB2_TGFBR1_TGFBR2', 'TGFB2_ACVR1_TGFBR1','SEMA6A_PLXNA4', 
                                            'SEMA3D_NRP1_PLXNA4','SEMA3D_NRP1_PLXNA2', 'SEMA3C_NRP1_PLXNA4','SEMA3C_NRP1_PLXNA2','SEMA3A_NRP1_PLXNA4',
                                            'SEMA3A_NRP1_PLXNA2','PTPRM_PTPRM','PLG_PARD3', 'NRXN3_NLGN1', 'NRG3_ERBB4','NRG3_ERBB2_ERBB4','NRG1_ERBB4',
                                            'NRG1_ERBB3','NRG1_ERBB2_ERBB4','NECTIN3_NECTIN2','NCAM1_NCAM1','MDK_SDC2','MDK_NCL','MDK_LRP1', 
                                            'LAMB1_ITGA1_ITGB1','LAMA4_ITGA1_ITGB1','LAMA2_ITGA1_ITGB1','LAMA1_ITGA1_ITGB1','JAM3_JAM3',
                                            'JAM3_F11R','FN1_ITGAV_ITGB1','F2_PARD3','F11R_JAM3','F11R_F11R','EFNA5_EPHA7','EFNA5_EPHA3', 
                                            'COL6A1_ITGA1_ITGB1','COL4A6_ITGA1_ITGB1','COL4A5_ITGA1_ITGB1','COL4A2_ITGA1_ITGB1','COL4A1_ITGA1_ITGB1', 
                                            'CDH2_CDH2','CDH1_CDH1', 'CADM1_CADM1','BMP5_BMPR1B_BMPR2','BMP5_BMPR1B_ACVR2A',
                                            'BMP5_BMPR1A_BMPR2','BMP5_BMPR1A_ACVR2A','BMP5_ACVR1_BMPR2','BMP5 -_ACVR1_ACVR2A', 'ANGPTL1_ITGA1_ITGB1')))
netVisual_chord_gene(cellchat, 
                     color.use = c('#C9BC9C','#DCDDDD'),
                     thresh = 0.05,
                     pairLR.use = pairLR,
                     lab.cex = 0.5,legend.pos.y = 30)