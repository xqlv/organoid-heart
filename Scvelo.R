library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyverse)
library(Seurat)
merged <- readRDS("~/merged.rds")

table(merged$group)
table(Idents(merged))

D8_CO<-subset(merged, subset = group %in% c("D8_CO"), invert = F)
D20_CO<-subset(merged, subset = group %in% c("D20_CO"), invert = F)

df_Co_D8 <- data.frame(Cells=Cells(Co_D8),sample=paste("Co_D8_SortedByCoordinate_Z3OKN",sep = ""))
df_Co_D20 <- data.frame(Cells=Cells(Co_D20),sample=paste("Co_D20_SortedByCoordinate_RKV38",sep = ""))
df <- rbind(df_Co_D8,df_Co_D20)
df$cell_id <- df$Cells

df$id<-sapply(df$Cells,function(x)paste(unlist(strsplit(x, "_"))[1],"x",sep = ""))
df$Cells<-paste(df$sample,df$id,sep = ":")
write.csv(df$Cells, file = "~/cellID_obs.csv", row.names = FALSE) 

cell_embeddings<-Embeddings(merged, reduction = "umap")
cell_embeddings <- as.data.frame(cell_embeddings)
cell_embeddings$cell_id <- rownames(cell_embeddings)
cell_embeddings <- merge(cell_embeddings, df)
row.names(cell_embeddings) <- cell_embeddings$Cells
cell_embeddings <- cell_embeddings[,c(2,3)]
write.csv(cell_embeddings, file = "~/cell_embeddings.csv")

clusters_obs<-merged$cluster 
clusters_obs <- as.data.frame(clusters_obs)
clusters_obs$cell_id <- rownames(clusters_obs)
clusters_obs <- merge(clusters_obs, df)
row.names(clusters_obs) <- clusters_obs$Cells
clusters_obs2 <-clusters_obs
clusters_obs <- clusters_obs[,-1]
clusters_obs <- clusters_obs[,-c(2,3,4)]
clusters_obs <- as.data.frame(clusters_obs)
rownames(clusters_obs) <- clusters_obs2$Cells
write.csv(clusters_obs, file = "~/clusters_obs.csv")