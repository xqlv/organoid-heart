library(CellChat)

cellchat_regiond <- readRDS("~/cellchat_regiond.rds")
cellchat_CCO <- readRDS("~/cellchat_CCO.rds")
cellchat_VCCO <- readRDS("~/cellchat_VCCO.rds")

cellchat_regiond <- netAnalysis_computeCentrality(cellchat_regiond, slot.name = "netP",thresh = 0.05)
cellchat_CCO <- netAnalysis_computeCentrality(cellchat_CCO, slot.name = "netP",thresh = 0.05)
cellchat_VCCO <- netAnalysis_computeCentrality(cellchat_VCCO, slot.name = "netP",thresh = 0.05)

object.list <- list(regiond = cellchat_regiond, CCO = cellchat_CCO, VCCO = cellchat_VCCO)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list))
cellchat

saveRDS(cellchat, file = "~/merged_cellchat.rds")

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, comparison = c(1,2,3),do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, comparison = c(1,2,3),do.stat = TRUE)
gg1 + gg2



