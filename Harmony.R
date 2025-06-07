library(tidyverse)
library(Seurat)
library(ggsci)

CCO_day20 <- readRDS("~/CCO_day20/reduction_dims30_res1.rds")
CCO_day8 <- readRDS("~/CCO_day8/reduction_dims30_res1.rds")

CCO_day20$cluster <- paste('CCO_day20_c', CCO_day20$seurat_clusters, sep = '')
CCO_day8$cluster <- paste('CCO_day8_c', CCO_day8$seurat_clusters, sep = '')
table(CCO_day20$cluster)
table(CCO_day8$cluster)

merged <- merge(CCO_day8, y = CCO_day20)
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged, 
                    # features = all.genes, 
                    vars.to.regress = c("percent.mt","nFeature_RNA", "nCount_RNA",
                                        "percent.rb","G2M.Score","S.Score"))

merged <- RunPCA(merged, features = VariableFeatures(object = merged))
Idents(merged) <- merged$cluster
DefaultAssay(merged)
merged <- merged |>
  harmony::RunHarmony(
    assay.use = "RNA",
    reduction.use = "pca",
    group.by.vars = "group",
    plot_convergence = TRUE
  )
DimPlot(merged,reduction = "harmony",split.by = "group")
VlnPlot(merged,features = "harmony_1",pt.size = 0)

ElbowPlot(merged, ndims = 50)
# t-SNE
merged <- RunTSNE(merged, reduction = "harmony", dims = 1:30)
DimPlot(merged, reduction = "tsne",group.by = "group")+scale_color_npg()
DimPlot(merged, reduction = "tsne",group.by = "group",split.by = "group")+scale_color_npg()
# UMAP
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30)
DimPlot(merged, reduction = "umap",group.by = "group")+scale_color_npg()
DimPlot(merged, reduction = "umap",group.by = "group",split.by = "group")+scale_color_npg()

DimPlot(merged, reduction = "umap",label = T,repel = T)

saveRDS(merged, '~/merged_dim30_harmony.rds')