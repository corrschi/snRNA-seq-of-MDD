#####-----------1.library-------------#######
library(harmony)
library(Seurat)
library(Matrix)
library(stringr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(parallel)
library(stringi)
library(plyr)
library(ggthemes)
library(cowplot)
library(data.table)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(viridis)
library(reshape2)
library(scales)
library(rlang)
library(dendextend)
library(DoubletFinder)
library(future)
library(ggsci)
library(ggpubr)
library(Matrix)
library(irlba)

#####-----------2.Multiplan-------------#######
options(future.globals.maxSize = 10*1000 * 1024^2)
plan("multicore", workers = 8)
options(future.globals.maxSize = 10000000000000)

#####-----------3.color-------------#######
color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10))[-8]
color_used1 <- c(pal_npg()(10),pal_igv()(9),pal_rickandmorty("schwifty")(12),pal_futurama()(12), pal_aaas()(10))[-5] 
