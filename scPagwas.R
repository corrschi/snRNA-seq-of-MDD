library(scPagwas)
library(readr)
library(dplyr)
library(plyr)
library(Seurat)


work.dir <- './scPages'
setwd(work.dir)

#######---1 single cell data prepare--------#####
load("/public/home/chidm/Workspace/ouyang/01_mdd/seurat/mdd_rmMix/filter1/mdd.filter1.Seurat.RData")

#subset mdd patients
sc.mdd <- subset(MERGE,subset=group=='Suicide')

#######---2 gwas data prepare--------#####
#Intersect the snp for all gwas file
library(readr)
library(dplyr)
gwas<-read_table("./MDD_gwas_data.txt") 
SNP_prune<- read_table("./EUR_LD0.8.prune") #keeped 2716380 SNPs after pruning
SNP_prune<-SNP_prune[!duplicated(unlist(SNP_prune)),] #（2716292）
colnames(SNP_prune)<-"rsid"

# Left Join using inner_join function 
gwas= gwas %>% inner_join(SNP_prune,by="rsid")
print(nrow(gwas))#（2716456）

gwas <- gwas[,c("chrom", "pos", "rsid", "beta", "se","maf")] %>% as.data.frame
gwas$pos <- gwas$pos %>% as.integer() #must integer，num will call error
write.table(gwas,file="./mdd_prune_gwas_data.txt",row.names=F,quote=F)

#######---3 scPagwas--------#####
sc.test <- sc.mdd

gwas.test <- read.table("./gwas_data/mdd_prune_gwas_data.txt",sep=' ',header=T)

Pagwas <- scPagwas_main(Pagwas = NULL,
                        gwas_data =gwas.test,
                        Single_data =sc.test,
                        output.prefix="MDD",
                        output.dirs="MDD_outputv2.1",
                        block_annotation = block_annotation,
                        assay="RNA",
                        Pathway_list=Genes_by_pathway_kegg,
                        chrom_ld = chrom_ld,
                        singlecell=T,
                        iters_singlecell = 10,
                        celltype=T) 
save(Pagwas,file="./MDD_outputv2.1/Pagwas_mdd.RData")
