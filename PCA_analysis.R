library(Seurat)
library(SeuratData)
library(tidyverse)
library(SeuratDisk)
library(ggsci)

adult_heart <- readRDS('~/adult_heart.rds')
epicardiods <- readRDS('~/epicardiods.rds')
Fetal_heart <- readRDS('~/Fetal_heart.rds')
Co_D8_heart <- readRDS('~/Co_D8_heart.rds')
Co_D20_heart <- readRDS('~/Co_D20_heart.rds')

genes <- Reduce(intersect, list(rownames(adult_heart), rownames(epicardiods), rownames(Fetal_heart),
                                rownames(Co_D8_heart), rownames(Co_D20_heart)))

merged_data <- cbind(adult_heart[genes,],epicardiods[genes,], Fetal_heart[genes,c(10,12,13,1,5,6)],
                     Co_D8_heart[genes,],Co_D20_heart[genes,]) %>% as.data.frame()
rownames(merged_data) <- genes
colnames(merged_data) <- c('adult_heart',
                           'epicardiods',
                           "Fetal_heart_HE5W","Fetal_heart_HE7W" ,"Fetal_heart_HE9W",
                           "Fetal_heart_HE10W", "Fetal_heart_HE20W", "Fetal_heart_HE22W",
                           'D8_heart',
                           'D20_heart')

merged_data <- t(merged_data)
merged_data <- merged_data[ , which(apply(merged_data, 2, var) != 0)]
pca_result <- prcomp(merged_data, scale. = TRUE) 

summary(pca_result)

pca_result$x

pca_result$rotation

biplot(pca_result)

library(ggplot2)
library(ggrepel)

pca_scores <- data.frame(pca_result$x)

ggplot(pca_scores, aes(PC1, PC2)) +
  geom_point() +
  ggtitle("PCA: First and Second Principal Components")+
  geom_text_repel(aes(label = rownames(pca_scores)))  

pc2 <- pca_result$x[, 2]

plot(pc2, main = "Second Principal Component", xlab = "Sample", ylab = "PC2", pch = 19)

library(ggplot2)
pca_df <- data.frame(PC2 = pc2, Sample = rep(1,length(pc2)))
ggplot(pca_df, aes(x = Sample, y = PC2)) +
  geom_point() +
  labs(title = "Second Principal Component", x = "Sample", y = "PC2")+
  geom_text_repel(aes(label = rownames(pca_scores)))  




