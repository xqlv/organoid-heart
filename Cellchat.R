library(tidyverse)
library(Seurat)
library(ggsci)

merged <- readRDS('~/merged.rds')
Idents(sub_merged) <- sub_merged$cell_type

output.dir <-
  paste0("~/cellchat/D8_D20", "/") # must have "/"
dir.create(output.dir, recursive = T)

seurat_obj <- subset(sub_merged,
                     subset = orig.ident %in% c('reginoid_day20'))

names(table(Idents(seurat_obj)))
print(table(Idents(seurat_obj)))

expr <- seurat_obj@assays$RNA@data

data.input <- expr
dim(data.input)
data.input[1:4, 1:4]
meta <- as.data.frame(Idents(seurat_obj))
colnames(meta) <- "labels"

unique(meta$labels) # check the cell labels
cellchat <-
  createCellChat(object = data.input,
                 meta = meta,
                 group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <-
  setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <-
  as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <-
  CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <-
  subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
# future::plan("multiprocess", workers = 10) # do parallel
cellchat <-
  identifyOverExpressedGenes(cellchat) # take a short time
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) # take a short time

#### Compute the communication probability and infer cellular communication network ----
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#### Extract the inferred cellular communication network as a data frame ----
#pass
#### Infer the cell-cell communication at a signaling pathway level ----
cellchat <- computeCommunProbPathway(cellchat)
#### Calculate the aggregated cell-cell communication network ----
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Interaction weights/strength"
)

saveRDS(cellchat, file = paste0(output.dir, "cellchat_reginoid_day20.rds"))


pairLR <- data.frame(interaction_name=rev(c('WNT5A_FZD6','WNT2B_FZD6_LRP6','VTN_ITGAV_ITGB1','TGFB2_ACVR1B_TGFBR2','TGFB1_ACVR1B_TGFBR2','SEMA6A_PLXNA4','SEMA6A_PLXNA2',
                                            'SEMA3D_NRP1_PLXNA4','SEMA3C_NRP1_PLXNA4','SEMA3B_NRP1_PLXNA4','SEMA3A_NRP1_PLXNA4','SEMA3A_NRP1_PLXNA2','PTPRM_PTPRM',
                                            'PTN_NCL','PRSS3_PARD3','PLG_PARD3','NRXN3_NLGN1','NRXN1_NLGN1','NRG4_ERBB4','NRG4_ERBB2_ERBB4','NRG3_ERBB4','NRG1_ERBB4',
                                            'NCAM1_NCAM1','NAMPT_INSR','MPZL1_MPZL1','MDK_NCL','JAG1_NOTCH2','FN1_ITGAV_ITGB1','F2_PARD3','EFNA5_EPHA7','EFNA5_EPHA3',
                                            'EFNA1_EPHA7','EFNA1_EPHA3','DLK1_NOTCH2','CDH2_CDH2','CADM1_CADM1','BMP5_BMPR1A_BMPR2','BMP2_BMPR1A_BMPR2','ANGPTL4_CDH11')))

netVisual_bubble(cellchat_D20, 
                 # sources.use = c(3,5,7,8,9), 
                 targets.use = c('CM'), 
                 # signaling = c('ncWNT','VTN',''),
                 pairLR.use = pairLR,
                 remove.isolate = FALSE)


pairLR <- data.frame(interaction_name=rev(c('WNT5B_FZD3','VTN_ITGAV_ITGB1','VCAM1_ITGA9_ITGB1','TGFB1_ACVR1_TGFBR1','SEMA6A_PLXNA2','SEMA3D_NRP1_PLXNA2','SEMA3A_NRP1_PLXNA2',
                                            'PTPRM_PTPRM','PTN_NCL','NRG3_ERBB4','NRG3_ERBB2_ERBB4','NRG2_ERBB4','NRG2_ERBB2_ERBB4','NRG1_ERBB4','NRG1_ERBB2_ERBB4',
                                            'NECTIN3_NECTIN2','NCAM1_NCAM1','NCAM1_FGFR1','NAMPT_INSR','LAMB1_DAG1','LAMA1-DAG1','JAM3_JAM3','JAG1_NOTCH2','HSPG2_DAG1',
                                            'FN1_ITGAV_ITGB1','F2_PARD3','F11R_JAM3','EFNB2_EPHA4','EFNA5_EPHA7','DSC2_DSG2','DLK1_NOTCH2','CDH4_CDH4','CDH2_CDH2',
                                            'BMP7_BMPR1A_BMPR2','BMP5_BMPR1A_BMPR2','BMP2_BMPR1A_BMPR2','ANGPTL4_CDH11')))
netVisual_bubble(cellchat_D8, 
                 # sources.use = c(3,5,7,8,9), 
                 targets.use = c('CM'), 
                 # signaling = c('ncWNT','VTN',''),
                 pairLR.use = pairLR,
                 remove.isolate = FALSE)