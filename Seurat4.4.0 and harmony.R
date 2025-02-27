#####-----------1.library and multiplan-------------#######
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
library(Matrix)
library(future)
library(ggsci)
library(ggpubr)
library(Matrix)
library(irlba)

options(future.globals.maxSize = 10*1000 * 1024^2)
plan("multicore", workers = 8)
options(future.globals.maxSize = 10000000000000)

#####-----------2.color-------------#######
color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10))[-8]
color_used1 <- c(pal_npg()(10),pal_igv()(9),pal_rickandmorty("schwifty")(12),pal_futurama()(12), pal_aaas()(10))[-5] 

####-----------3.data load----------########
setwd("./GSE144136/re_analysis")

Counts <- Matrix::readMM("./download_data/GSE144136_GeneBarcodeMatrix_Annotated.mtx")
my_summary <- summary(Counts)

gene_name <- read.csv("./GSE144136/GSE144136_GeneNames.csv")
cell_name <- read.csv("./GSE144136/GSE144136_CellNames.csv")

Counts <- Counts %>% as.sparse
colnames(Counts) <- cell_name$x
rownames(Counts) <- gene_name$x

load("./mdd_count_and_meta.data.RData")

####-----------4.Seurat and Harmony----------########
#create SeuratObject
mdd <- CreateSeuratObject(Counts, project = "MDD", min.features = 200)
mdd@meta.data = mdd_meta.data 

#remove Mix cluster
mdd <- subset(mdd,subset=Anno.paper %in% c("Mix_1","Mix_2","Mix_3","Mix_4","Mix_5"),invert=T)

#harmony and seurat process
dim.use <- 30
res.use <- 0.5
seed.use <- 88
object.name <- "MDD.rmMix" #修改

MERGE <- mdd

MERGE <- NormalizeData(object = MERGE, normalization.method = "LogNormalize", scale.factor = 1e4)
MERGE <- FindVariableFeatures(object = MERGE, selection.method = 'vst', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf), nfeatures = 2000)
MERGE <- ScaleData(object = MERGE, features = rownames(x = MERGE), verbose = T)
MERGE <- RunPCA(object = MERGE, features = VariableFeatures(object = MERGE), verbose = FALSE)

MERGE <- RunHarmony(MERGE,c("sample"), plot_convergence = TRUE,verbose = FALSE) 

MERGE <- FindNeighbors(MERGE, reduction = "harmony", verbose = FALSE, dims = 1:dim.use)
MERGE <- FindClusters(MERGE, resolution = res.use, verbose = FALSE, random.seed = seed.use) 

MERGE <- RunUMAP(MERGE, dims = 1:dim.use, reduction = "harmony",umap.method = "uwot",
                 n.neighbors = 20L, min.dist = 0.5)
mdd.rmMix <- MERGE

#Plots
DimPlot(MERGE,pt.size = 0.05,label=T,label.size=4,cols = color_used,raster=F)

#FindMarkers
MERGE.markers <- FindMarkers_parallel(MERGE, mc.cores = 5)
MERGE.markers %>% TOP_N(50, pct.1 = 0.2) -> top50
MERGE.markers <- MERGE.markers %>% TOP_N(5000)
write.table(top50,file = paste0(object.name,"_res", res.use ,"_dim", dim.use ,"_culster_top50_DEGs.csv"),sep = ",", row.names = T, quote = F)
write.table(MERGE.markers,file = paste0(object.name,"_res",res.use ,"_dim", dim.use ,"_culster_all_DEGs.csv"),sep = ",",row.names = T,quote = F)


mdd.filter1 <- subset(mdd.rmMix,idents=15,invert=T) #c15 is mix
#repeat harmony,seurat and findmarkes process for mdd.filter1

save(mdd.filter1,file = "mdd.RData")
save(MERGE.markers,file='mdd.markers.RData')
