
###-------------1. data loading-------------######
work.dir="./Source data"
setwd(work.dir)

load("./mdd.filter1.Seurat.RData")
MERGE <- mdd.filter1
### check data
png(paste0("mdd.filter1_all_by_harmony_dimplot-ANNO.chi-CheckData.png"), width = 9, height = 5, res = 400, units = "in")
DimPlot(MERGE,pt.size = 0.1,label=T,label.size=3,cols = color_used,raster=F) 
dev.off()

MERGE #30062 features across 71565 samples 

####-------------2.compare group cell propotion by wilcox test--------#####
dat <- MERGE
dat$cluster <- dat$Anno.chi1 %>% as.character()
Idents(dat) <- "cluster"
dat$sample <- dat$sample %>% as.character()

library(ggpubr)
library(viridis)
Box_dat_produce <- function(x){
  bb <- table(x$sample,x$cluster)%>%as.matrix
  bb_rowSum <- rowSums(bb)
  join_col <- rep(bb_rowSum,length(table(x$cluster)))
  bb <- as.data.frame(bb)
  bb <- cbind(bb,join_col)
  colnames(bb) <- c("sample","variable","freq","sum")
  percent <- bb$freq/bb$sum
  bb <- cbind(bb,percent)
  group <- rep(x$group%>%as.character%>%unique,length(rownames(bb))) 
  bb <- cbind(bb,group)
}


data1 <- subset(dat,subset=group=="Control")
data2 <- subset(dat,subset=group=="Suicide")


dat <- rbind(Box_dat_produce(data1),
             Box_dat_produce(data2))


dat$group <- factor(dat$group,levels=c("Control","Suicide"))

my_comparisons <- list( c("Control","Suicide"))

dat_plot <- dat

# Plots
library(ggpubr)
library(viridis)
png("Fig1g.Propertion of cell clusters_wilcox.test-facet.png", width = 12, height = 4, units = "in", res = 400)
ggboxplot(dat_plot,
          x = "variable", y = "percent",
          color = "group", palette = "jama",
          add = "point", outlier.colour = NULL) +
  facet_wrap(~ variable, scales = "free_x", ncol = 7) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, family = "Times New Roman"),  
    axis.text.y = element_text(size = 12, family = "Times New Roman"),  
    axis.title.x = element_text(size = 14, family = "Times New Roman"), 
    axis.title.y = element_text(size = 14, family = "Times New Roman"), 
    legend.position = "right", 
    legend.text = element_text(size = 12, family = "Times New Roman"),  
    legend.title = element_text(size = 14, family = "Times New Roman"),  
    strip.text = element_text(size = 12, family = "Times New Roman") 
  ) +
  stat_compare_means(aes(group = group), 
                     label = "p.signif", 
                     method = "wilcox.test",
                     label.y = c(0.8)) +  
  scale_colour_manual(values = c("#264653", "#E69F00")) +
  labs(x = "", y = "percent")

dev.off()

####-------------3.wilcox test-Padj--------#####
p_values.BH <- compare_means(percent ~ group, 
                             data = dat_plot, 
                             group.by = "variable", 
                             method = "wilcox.test",
                             p.adjust.method = "BH")  


print(p_values.BH)
write.csv(p_values.BH,file='Fig1g.p_values.BH.csv')

