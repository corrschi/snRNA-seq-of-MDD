
library(BiocManager)
library(GSEABase)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

setwd("./01_mdd/gsea")

#2.1 load deg data[Suicide vs Control]
deg <- read.csv("./Source data/differential_expression_results_all_clusters.csv",
                header=T)

deg.excit <- deg%>%dplyr::filter(cluster=='Excit_neuron') 

#2.2 gsea 

geneList= deg.excit$avg_log2FC 
names(geneList)= toupper(rownames(deg.excit))
geneList=sort(geneList,decreasing = T)
head(geneList)
gmtfile ='./Data_Reference/gsea/h.all.v2024.1.Hs.symbols.gmt'
geneset <- read.gmt( gmtfile )
length(unique(geneset$term))
egmt <- GSEA(geneList, TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.05, 
             nPerm = 1000,   
             verbose=FALSE)
head(egmt)
egmt@result 
gsea_results_df <- egmt@result 
rownames(gsea_results_df)
write.csv(gsea_results_df,file = 'Fig8B.gsea_results_df.csv')

#2.3 plots
p <- gseaplot2(egmt,geneSetID = 'HALLMARK_PI3K_AKT_MTOR_SIGNALING',pvalue_table=T)
ggsave("Fig8B.GSEA-PI3K_AKT.pdf", p, width = 10, height = 5)

