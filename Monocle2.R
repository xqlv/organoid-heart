library(Seurat)
library(ggplot2)
library(tidyverse)
library(monocle) 
library(dplyr)
library(RColorBrewer)

merged <- readRDS('~/merged.rds')
merged_CM <- subset(merged, subset = pre_cell_type %in% c('CM'))
merged_CM <- RunUMAP(merged_CM, reduction = "harmony", dims = 1:30)
DimPlot(merged_CM, 
        group.by = "group_seurat_clusters",
        reduction = "umap", # tsne, umap, pca
        label = T) + scale_color_manual(values = c('#00BF77','#00BF93','#00C1AC','#00BFC3','#00BCD8','#00B5EB','#00AEF9','#05A3FF',
                                                   '#9286BE','#A07EB7','#AD77B1','#BC71AC','#CA6CA8','#DB67A3','#EA639D'))
DimPlot(merged_CM, 
        group.by = "group",
        reduction = "umap", # tsne, umap, pca
        label = T) + scale_color_manual(values = c('#1A65A5','#717070'))

metadata <- merged_CM@meta.data
table(merged_CM$group_seurat_clusters)
Idents(merged_CM) <- merged_CM$group_seurat_clusters

data <- as(as.matrix(merged_CM@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = merged_CM@meta.data)
fData <-
  data.frame(gene_short_name = row.names(data),
             row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

cds <- newCellDataSet(
  data,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <-
  detectGenes(cds, min_expr = 0.1) 
head(fData(cds))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) 

DefaultAssay(merged_CM) <- "RNA"
merged_CM <- NormalizeData(merged_CM)
merged_CM@assays[["RNA"]]@data

deg <- FindAllMarkers(merged_CM,
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
# cds <- orderCells(cds, root_state = 1)
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
                     show_backbone = TRUE)+ scale_color_manual(values = c('#00BF77','#00BF93','#00C1AC','#00BFC3','#00BCD8','#00B5EB','#00AEF9','#05A3FF',
                                                                          '#9286BE','#A07EB7','#AD77B1','#BC71AC','#CA6CA8','#DB67A3','#EA639D'))

Time_genes <- c('PLK2','CKMT2','CGNL1','PLN','DHRS3','RABGAP1L','CXCL12','IRX1','HEY2','GJA5','IRX2','IRX3')
library(RColorBrewer)
cols <- c(brewer.pal(9, "YlGnBu"))
p <- plot_pseudotime_heatmap(
  cds[Time_genes,],
  # cds,
  num_clusters = 1,
  # add_annotation_col = 'group_seurat_clusters',
  cluster_rows = FALSE,
  hmcols = colorRampPalette(rev(cols))(100),
  show_rownames = T,
  return_heatmap = T
)