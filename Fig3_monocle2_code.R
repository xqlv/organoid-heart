library(monocle) 
library(Seurat) 
library(dplyr) 
library(tidyverse)

merged_EP <- readRDS('~/merged_EP.rds')

metadata <- merged_EP@meta.data
table(merged_EP$group_seurat_clusters)
Idents(merged_EP) <- merged_EP$group_seurat_clusters

data <- as(as.matrix(merged_EP@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = merged_EP@meta.data)
fData <-
  data.frame(gene_short_name = row.names(data),
             row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

cds <- newCellDataSet(
  data,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
  )

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <-
  detectGenes(cds, min_expr = 0.1) 
head(fData(cds))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) 

DefaultAssay(merged_EP) <- "RNA"
merged_EP <- NormalizeData(merged_EP)
merged_EP@assays[["RNA"]]@data

deg <- FindAllMarkers(merged_EP,
                      assay = "RNA",only.pos = T)
top_deg <- deg %>% 
  filter(p_val_adj <= 1e-5) %>% 
  arrange("avg_log2FC") %>% 
  group_by(cluster) %>% 
  top_n(100, wt = avg_log2FC) %>% 
  pull(gene)

write.table(
  deg,
  file = "~/monocle_deg.txt",
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
)

ordergene <- top_deg
cds <- setOrderingFilter(cds, ordergene)
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 2) 
saveRDS(cds, "~/cds.rds")

plot_cell_trajectory(cds,
                     color_by = "Pseudotime",
                     size = 1,
                     show_backbone = TRUE)

plot_cell_trajectory(cds,
                     color_by = "group",
                     size = 1,
                     show_backbone = TRUE)

plot_cell_trajectory(cds,
                     color_by = "State",
                     size = 1,
                     show_backbone = TRUE)

plot_cell_trajectory(cds,
                     color_by = "group_seurat_clusters",
                     size = 1,
                     show_backbone = TRUE)+ scale_color_manual(values = c('#7030A0','#9E5ECE','#47A7C9','#318AA9'))

keygenes <- c('TNNT1','SFRP2','SFRP5',
              'LRP2','CALB2','C3',
              'HAND1', 'MAB21L2', 'HOXB6', 'HOXB5', 'BNC2',
              'HEY1','HEY2','CFB','C1R','C1S','ELN','DPT') 

pdata <- Biobase::pData(cds)
s.cells <-
  subset(pdata#, subset = group %in% c("FTO_NK")
  ) %>% rownames()
cds_subset <- cds[keygenes, s.cells]

plot_genes_in_pseudotime(cds_subset, color_by = "State")
plot_genes_in_pseudotime(cds_subset, color_by = "group",ncol = 4)
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")


Time_genes <- c('MAB21L2','HAND1','HOXB6','BNC2','HOXB5',
                'TNNT1','SFRP5','HEY1','LRP2','CALB2',
                'SFRP2','HEY2','C1S','CFB','C1R','C3')
library(RColorBrewer)
cols <- c(brewer.pal(9, "YlGnBu"))
p <- plot_pseudotime_heatmap(
  cds[unique(Time_genes),],
  # cds,
  num_clusters = 1,
  # add_annotation_col = 'group_seurat_clusters',
  cluster_rows = FALSE,
  hmcols = colorRampPalette(rev(cols))(100),
  show_rownames = T,
  return_heatmap = T
)


library(ggplot2)
df <- pdata
sorted_index <- order(df$Pseudotime)
df_sorted <- df[sorted_index, ]
df_sorted$rank <- seq_along(sorted_index)
df_sorted$rank <- as.numeric(df_sorted$rank)

ggplot(df_sorted, aes(x = rank, y = pre_cell_type)) +
  geom_bar(aes(fill = group_seurat_clusters),stat = "identity") +   
  theme_minimal() +
  scale_fill_manual(values = c('#7030A0','#9E5ECE','#47A7C9','#318AA9'))



table(merged_EP$group_seurat_clusters)
merged_EP$group_seurat_clusters <- factor(merged_EP$group_seurat_clusters, levels = c('D8_CO_0','D8_CO_7','D20_CO_7','D20_CO_10'))
features <- c('MAB21L2','HAND1','HOXB6','BNC2','HOXB5',
              'TNNT1','SFRP5','HEY1','LRP2','CALB2',
              'SFRP2','HEY2','C1S','CFB','C1R','C3')
avg_exp_merged_EP <- AverageExpression(merged_EP,
                                       assays = 'RNA',
                                       features = unique(features),
                                       return.seurat = FALSE,
                                       group.by = "group_seurat_clusters",
                                       slot = "data") %>% as.data.frame()

library(RColorBrewer)
cols <- c(brewer.pal(11, "RdBu"))
pheatmap::pheatmap(
  avg_exp_merged_EP,
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
  color = colorRampPalette(rev(cols))(100),
  border_color = "white")

table(merged_EP$group)
merged_EP$group <- factor(merged_EP$group, levels = c('D8_CO','D20_CO'))
avg_exp_merged_EP <- AverageExpression(merged_EP,
                                       assays = 'RNA',
                                       features = unique(features),
                                       return.seurat = FALSE,
                                       group.by = "group",
                                       slot = "data") %>% as.data.frame()

library(RColorBrewer)
cols <- c(brewer.pal(11, "RdBu"))
pheatmap::pheatmap(
  avg_exp_merged_EP,
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
  color = colorRampPalette(rev(cols))(100),
  # color = colorRampPalette(rev(cols))(100),
  border_color = "white")

