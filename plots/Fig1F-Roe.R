###-------------1. data loading-------------######
work.dir="./Source data"

setwd(work.dir)
load("./mdd.filter1.Seurat.RData")
MERGE <- mdd.filter1
# check data
png(paste0("mdd.filter1_all_by_harmony_dimplot-ANNO.chi-CheckData.png"), width = 9, height = 5, res = 400, units = "in")
DimPlot(MERGE,pt.size = 0.1,label=T,label.size=3,cols = color_used,raster=F) 
dev.off()

######------2. ROE---------#####
roeObject <- MERGE
roeObject$cluster <-roeObject$Anno.chi1 %>% as.character() 

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

# ------theme set-----------
text.size = 14
text.angle = 0
text.hjust = 0.5  
legend.position = "right"

mytheme <- theme(
  plot.title = element_text(size = text.size + 2, color = "black", hjust = 0.5, family = "Times New Roman"),
  axis.title = element_text(size = text.size, color = "black", family = "Times New Roman"),                 
  axis.text = element_text(size = text.size, color = "black", family = "Times New Roman"),                
  axis.text.x = element_text(angle = text.angle, hjust = text.hjust, family = "Times New Roman"),       
  axis.text.y = element_text(size = text.size, family = "Times New Roman"),                             
  panel.grid = element_blank(),                                                                           
  legend.position = legend.position,                                                                     
  legend.text = element_text(size = text.size, family = "Times New Roman"),                             
  legend.title = element_text(size = text.size, family = "Times New Roman")                              
)

# Plots
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
