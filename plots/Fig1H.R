
###-------------1. data loading-------------######
work.dir="./Source data"
setwd(work.dir)

load("./mdd.filter1.Seurat.RData")
MERGE <- mdd.filter1
### check data
png(paste0("mdd.filter1_all_by_harmony_dimplot-ANNO.chi-CheckData.png"), width = 9, height = 5, res = 400, units = "in")
DimPlot(MERGE,pt.size = 0.1,label=T,label.size=3,cols = color_used,raster=F) 
dev.off()
####-------2.BarPlots------------####
library(ggsci)
library(ggplot2)

dat <- MERGE
dat$cluster <- factor(dat$Anno.chi1,levels = c("Astro", "Endo","Excit_neuron", "Inhib_neuron", "Micro/Macro", "Oligo",
                                               "OPC"))
Idents(dat) <- "cluster"
# dat$sample <- dat$sample %>% as.character()
dat$group <- factor(dat$group,levels= c("Control","Suicide"))
dat$sample_num <- gsub("sample","",dat$sample,ignore.case=F)

## sample change to only number
aa <- dat@meta.data[,c("sample","group")]
bb <- aa%>%group_by(group)%>%unique%>%as.data.frame
bb$sample_num <- bb$sample%>%gsub("sample","",.)
control <- bb%>%filter(group=='Control')%>%.$sample_num

suicide <- bb%>%filter(group=='Suicide')%>%.[["sample_num"]]

dat$sample_num <- factor(dat$sample_num,
                         levels =  c(bb%>%filter(group=='Control')%>%.$sample_num%>%as.numeric%>%sort,
                                     bb%>%filter(group=='Suicide')%>%.$sample_num%>%as.numeric%>%sort ))



dat_plot <- dat


color26 <- c("#D9DD6B","#ECEFA4","#D54C4C","#8D2828","#FDD2BF","#E98580","#DF5E5E","#492F10","#334257","#476072","#548CA8","#00A19D","#ECD662","#5D8233","#284E78","#3E215D","#835151","#F08FC0","#C6B4CE","#BB8760","#FFDADA","#3C5186","#558776","#E99497","#FFBD9B","#0A1D37")


# pdf("Fig1h.clusters Barplots-sample_and_cluster-Anno.chi1.pdf",width = 12,height =6)
png("Fig1H.clusters Barplots-sample_and_cluster-Anno.chi1.png", width = 12, height = 6, units = "in", res = 400)
ggplot(data=dat_plot@meta.data, aes(x=sample_num, fill=cluster)) +
  geom_bar(stat="count", position="fill") +
  labs(x="", y="Proportion of cells") +
  scale_fill_manual(values=color26) +
  theme_bw() +
  theme(
    panel.border = element_rect(color="black", size=1, fill=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour="black"),
    axis.text.x = element_text(angle=0, size=12, family="Times New Roman"),
    axis.text.y = element_text(angle=0, hjust=1, vjust=0.95, size=12, family="Times New Roman"),
    axis.title.x = element_text(size=14, family="Times New Roman"),  
    axis.title.y = element_text(size=14, family="Times New Roman"),  
    legend.text = element_text(size=14, family="Times New Roman"),  
    legend.title = element_text(size=16, family="Times New Roman"), 
    plot.title = element_text(size=16, family="Times New Roman", face="bold")  
  )

dev.off()



######------3.export Source data---------#####
dat <- MERGE
dat$cluster <- factor(dat$Anno.chi1,levels = c("Excit_neuron", "Inhib_neuron", "Oligo", "Astro", "Micro/Macro", 
                                               "OPC", "Endo"))
Idents(dat) <- "cluster"
# dat$sample <- dat$sample %>% as.character()
dat$group <- factor(dat$group,levels= c("Control","Suicide"))

dat <- dat@meta.data[,c("cluster","group","sample")]
write.csv(dat, file = "Fig1H.Source data.csv", row.names =T)
