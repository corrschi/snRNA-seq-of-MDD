
###-------------1. data loading-------------######
work.dir="./Source data"

setwd(work.dir)
load("./mdd.filter1.Seurat.RData")
subset_cells <- mdd.filter1
DefaultAssay(subset_cells) <- "RNA"
subset_cells@assays$RNA@data <- subset_cells@assays$RNA@counts
subset_cells <- NormalizeData(subset_cells)
subset_cells <- ScaleData(subset_cells, features=VariableFeatures(subset_cells))

MERGE <- subset_cells
### recovery
png(paste0("mdd.filter1_all_by_harmony_dimplot-ANNO.chi-CheckData.png"), width = 9, height = 5, res = 400, units = "in")
DimPlot(MERGE,pt.size = 0.1,label=T,label.size=3,cols = color_used,raster=F) 
dev.off()

MERGE #30062 features across 71565 samples 

####-------------2.细胞比例差异 wilcox test--------#####
dat <- MERGE
dat$cluster <- dat$Anno.chi1 %>% as.character()
Idents(dat) <- "cluster"
dat$sample <- dat$sample %>% as.character()

library(ggpubr)
library(viridis)
Box_dat_produce <- function(x){
  bb <- table(x$sample,x$cluster)%>%as.matrix
  bb_rowSum <- rowSums(bb)
  join_col <- rep(bb_rowSum,length(table(x$cluster))) #rep 分群数量，需修改
  bb <- as.data.frame(bb)
  bb <- cbind(bb,join_col)
  colnames(bb) <- c("sample","variable","freq","sum")
  percent <- bb$freq/bb$sum
  bb <- cbind(bb,percent)
  group <- rep(x$group%>%as.character%>%unique,length(rownames(bb))) #根据分组更改
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
  facet_wrap(~ variable, scales = "free_x", ncol = 7) +  # 按 cluster 分面
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, family = "Times New Roman"),  # x 轴刻度文字
    axis.text.y = element_text(size = 12, family = "Times New Roman"),  # y 轴刻度文字
    axis.title.x = element_text(size = 14, family = "Times New Roman"),  # x 轴标题
    axis.title.y = element_text(size = 14, family = "Times New Roman"),  # y 轴标题
    legend.position = "right",  # 将图例位置设置为右侧
    legend.text = element_text(size = 12, family = "Times New Roman"),  # 图例文字
    legend.title = element_text(size = 14, family = "Times New Roman"),  # 图例标题
    strip.text = element_text(size = 12, family = "Times New Roman")  # 分面标题文字
  ) +
  stat_compare_means(aes(group = group), 
                     label = "p.signif",  # 显示显著性符号
                     method = "wilcox.test",
                     label.y = c(0.8)) +  # 调整显著性符号的 y 位置
  scale_colour_manual(values = c("#264653", "#E69F00")) +
  labs(x = "", y = "percent")

dev.off()

####-------------3.wilcox test-Padj--------#####
p_values.BH <- compare_means(percent ~ group, 
                             data = dat_plot, 
                             group.by = "variable", 
                             method = "wilcox.test",
                             p.adjust.method = "BH")  # 选择 FDR 校正方法，默认是holm


print(p_values.BH)
write.csv(p_values.BH,file='Fig1g.p_values.BH.csv')

###---------------4.输出Source data-------#####
write.csv(dat, file = "Fig1g.Source data.csv", row.names = FALSE)

