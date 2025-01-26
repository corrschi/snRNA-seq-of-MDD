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

MERGE 

######------2. ROE---------#####
roeObject <- MERGE
roeObject$cluster <-roeObject$Anno.chi1 %>% as.character() ##去除分群的factor

bb <- table(roeObject$group,roeObject$cluster)  %>% unclass()
# R_O_E_my <- CrossTable(bb %>% as.matrix(), expected = T)
# Multiple_CD4 <- (R_O_E_CD4$CST$observed/R_O_E_CD4$CST$expected)[, 2]
cc <- lapply(rowSums(bb), function(x){
  
  x * colSums(bb) %>% as.data.frame()
})
ww <- Reduce(cbind,cc)
colnames(ww) <- names(cc)
www <- ww/sum(colSums(bb),rowSums(bb))
ROE <- t(bb)/(www)
save(ROE,file = "Fig1F.ROE_all.RData")

### plot ROE
library(ggplot2)
library(hrbrthemes)
library(ggpubr)
library(viridis)

dat_plot <- log2(ROE+0.01)
dat_plot$cluster <- dat_plot %>% rownames()
dat_plot <- melt(dat_plot)
colnames(dat_plot) <- c("cluster", "group", "Roe")

# ------主题设置-----------
text.size = 14
text.angle = 0
text.hjust = 0.5  # 将横轴文字居中
legend.position = "right"

mytheme <- theme(
  plot.title = element_text(size = text.size + 2, color = "black", hjust = 0.5, family = "Times New Roman"), # 标题文字
  axis.title = element_text(size = text.size, color = "black", family = "Times New Roman"),                  # 坐标轴标题文字
  axis.text = element_text(size = text.size, color = "black", family = "Times New Roman"),                  # 坐标轴刻度文字
  axis.text.x = element_text(angle = text.angle, hjust = text.hjust, family = "Times New Roman"),           # 横轴刻度文字，居中
  axis.text.y = element_text(size = text.size, family = "Times New Roman"),                                # 纵轴刻度文字
  panel.grid = element_blank(),                                                                            # 去网格线
  legend.position = legend.position,                                                                       # 图例位置
  legend.text = element_text(size = text.size, family = "Times New Roman"),                                # 图例文字
  legend.title = element_text(size = text.size, family = "Times New Roman")                                # 图例标题
)

# 绘图代码
png("Fig1F.Roe.png", width = 6, height = 5, units = "in", res = 400)

ggplot(data = dat_plot, aes(x = group, y = cluster, fill = Roe)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(
    aes(label = round(Roe, 2)),
    color = "black",
    size = 5,
    show.legend = TRUE
  ) +
  theme_minimal() +
  scale_fill_gradientn(colors = c("#6d8bc3", "white", "#e3716e")) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  mytheme

dev.off()


######------3.输出Source data---------#####
dat <- roeObject@meta.data[,c("cluster","group")]
write.csv(dat, file = "Fig1F.Source data.csv", row.names =T)