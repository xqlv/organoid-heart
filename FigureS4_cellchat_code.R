library(ggplot2)
library(tidyverse)
library(ggsci)
library(CellChat)
library(patchwork)

cellchat_D8 <- readRDS('~/cellchat/D8_D20_CO/cellchat_D8_CO_4celltypes.rds')
netVisual_bubble(cellchat_D8, 
                 # sources.use = 4, 
                 # targets.use = c('CM'), 
                 remove.isolate = FALSE)

cellchat_D20 <- readRDS('~/cellchat/D8_D20_CO/cellchat_D20_CO_4celltypes.rds')
netVisual_bubble(cellchat_D20, 
                 # sources.use = 4, 
                 # targets.use = c('CM'), 
                 remove.isolate = FALSE)
