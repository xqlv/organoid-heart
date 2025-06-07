library(Seurat)
library(SeuratData)
library(tidyverse)
library(SeuratDisk)
library(ggsci)

adult_heart <- LoadH5Seurat("~/2020-nature-Litvinukova.h5seurat",
                            assay ="RNA",slots='counts')
table(adult_heart$cell_type)
DimPlot(adult_heart,reduction = 'umap',raster=FALSE,group.by = 'cell_type')

adult_heart <- subset(adult_heart,
                      subset = cell_type %in% c("Adipocytes",'Atrial_Cardiomyocyte',
                                                'Endothelial','Fibroblast',
                                                'Lymphoid','Mesothelial',
                                                'Myeloid','Neuronal',
                                                'Pericytes','Smooth_muscle_cells',
                                                'Ventricular_Cardiomyocyte'))
DimPlot(adult_heart,reduction = 'umap',raster=FALSE,group.by = 'cell_type')
adult_heart1 <- adult_heart

set.seed(717)
Idents(adult_heart) <- adult_heart$cell_type
adult_heart = adult_heart[, sample(1:ncol(adult_heart),round(ncol(adult_heart)/4))]

adult_heart <- NormalizeData(adult_heart)
adult_heart <- FindVariableFeatures(adult_heart)
adult_heart <- ScaleData(adult_heart)
adult_heart <- RunPCA(adult_heart, npcs = 50, verbose = FALSE)
anchors <- FindTransferAnchors(reference = adult_heart, query = sub_merged_heart, dims = 1:50,
                               reference.reduction = "pca")
adult_heart <- RunUMAP(adult_heart, dims = 1:50, reduction = "pca", return.model = TRUE)
query <- MapQuery(anchorset = anchors, reference = adult_heart, query = sub_merged_heart,
                  refdata = list(celltype = "cell_type"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(adult_heart, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 3,raster=FALSE,
              repel = TRUE) + ggtitle("Reference annotations")
p2 <- DimPlot(query, reduction = "ref.umap", group.by = "cell_type", label = TRUE,raster=FALSE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transferred labels")
p1 + p2

merged_mapping_heart <- merge(adult_heart,query)
merged_mapping_heart[["umap"]] <- merge(adult_heart[["umap"]], query[["ref.umap"]])

DimPlot(merged_mapping_heart, group.by = 'cell_type', shuffle = TRUE,raster=FALSE)
merged_mapping_heart$label <- merged_mapping_heart$cell_type
merged_mapping_heart$label[merged_mapping_heart$cell_type %in% c("Adipocytes",'Atrial_Cardiomyocyte',
                                                                 'Endothelial','Fibroblast',
                                                                 'Lymphoid','Mesothelial',
                                                                 'Myeloid','Neuronal',
                                                                 'Pericytes','Smooth_muscle_cells',
                                                                 'Ventricular_Cardiomyocyte')] <- "reference"
table(merged_mapping_heart$label)
library(RColorBrewer)
cols <- brewer.pal(12, "Paired")
c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
p3 <- DimPlot(merged_mapping_heart, reduction = "umap", group.by = "label", label = TRUE,raster=FALSE,
              cols = c("#A6CEE3", '#329947','#D62A85', '#7A3489', 'grey'),
              label.size = 3, repel = TRUE) 
p1 + p3