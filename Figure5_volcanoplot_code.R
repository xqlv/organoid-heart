library(tidyverse)
library(ggsci)
library(RColorBrewer)

markers <- readRDS('~/deg_CM_D20VCO_D20CO/markers.rds')

source("~/version3/custom_plot_function.R")
p <- markers |>
  catvolcano(
    x = avg_log2FC,
    y = p_val_adj,
    label = gene,
    p_value = 0.05,
    log2FC = 0.05,
    text = c("FTL","PLN","CKMT2","TTN")
  )+
  # scale_color_manual(values=c("#00b0eb","#ffd401","#e20612"))+
  geom_point(alpha = 0.8)
p