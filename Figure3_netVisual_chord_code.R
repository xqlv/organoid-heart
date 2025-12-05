library(CellChat)

cellchat.D8_CO <- readRDS("~/cellchat_D8_CO_4celltypes.rds")
cellchat.D20_CO <- readRDS("~/cellchat_D20_CO_4celltypes.rds")

cellchat.D8_CO <- netAnalysis_computeCentrality(cellchat.D8_CO, slot.name = "netP",thresh = 0.05)
cellchat.D20_CO <- netAnalysis_computeCentrality(cellchat.D20_CO, slot.name = "netP",thresh = 0.05)

pairLR <- data.frame(interaction_name=rev(c('VEGFA_VEGFR1R2','PTPRM_PTPRM',
                                            'NCAM1_NCAM1','MDK_NCL',
                                            'FN1_ITGAV_ITGB1',
                                            'COL1A1_ITGA9_ITGB1','CDH2_CDH2','CADM1_CADM1',
                                            'BMP5_BMPR1A_BMPR2')))
netVisual_chord_gene(cellchat.D8_CO, pairLR.use = pairLR,
                     color.use = c('grey','#377EB8','#984EA3','#4DAF4A'),
                     lab.cex = 0.7,legend.pos.y = 30)

pairLR <- data.frame(interaction_name=rev(c('VEGFA_VEGFR1R2','PTPRM_PTPRM','NRG3_ERBB4',
                                            'JAG1_NOTCH1','FN1_ITGAV_ITGB1',
                                            'COL1A1_ITGA1_ITGB1','CDH2_CDH2','CADM1_CADM1',
                                            'BMP5_BMPR1A_BMPR2')))
netVisual_chord_gene(cellchat.D20_CO, pairLR.use = pairLR,
                     color.use = c('grey','#377EB8','#984EA3','#4DAF4A'),
                     thresh = 0.05,
                     lab.cex = 0.7,legend.pos.y = 30)




