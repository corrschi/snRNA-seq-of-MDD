# Fig1B
dat_plot <- MERGE
png(paste0("Fig1B.UMAP by group.png"), width = 6, height = 5, res = 400, units = "in")
DimPlot(dat_plot, pt.size = 0.3, cols = color_used, raster = FALSE, group = 'group') + 
  theme(
    text = element_text(family = "Times New Roman"),         # 设置全局字体为 Times New Roman
    axis.title = element_text(size = 12, family = "Times New Roman"),  # 坐标轴标题字体
    axis.text = element_text(size = 10, family = "Times New Roman"),   # 坐标轴刻度字体
    legend.text = element_text(size = 12, family = "Times New Roman"), # 图例文字
    legend.title = element_text(size = 12, family = "Times New Roman") # 图例标题
  )
dev.off()



# dat_plot$sample <- factor(dat_plot$sample,levels = c("sample1","sample2", "sample3","sample4", "sample5", "sample6",
#                                                      "sample7", "sample8",  "sample9","sample10", "sample11",
#                                                      "sample12", "sample13","sample14", "sample15", "sample16", "sample17",
#                                                      "sample18", "sample19",  "sample20", "sample21", "sample22",
#                                                       "sample23", "sample24",  "sample25", "sample26", "sample27",
#                                                      "sample28", "sample29",  "sample30", "sample31", "sample32", 
#                                                      "sample33", "sample34"))
png(paste0("Fig1B.UMAP by sample.png"), width = 9, height = 5, res = 400, units = "in")
DimPlot(dat_plot, pt.size = 0.1, cols = color_used[c(1:34)], raster = FALSE, group = 'sample') + 
  theme(
    text = element_text(family = "Times New Roman"),         # 设置全局字体为 Times New Roman
    axis.title = element_text(size = 12, family = "Times New Roman"),  # 坐标轴标题字体
    axis.text = element_text(size = 10, family = "Times New Roman"),   # 坐标轴刻度字体
    legend.text = element_text(size = 12, family = "Times New Roman"), # 图例文字
    legend.title = element_text(size = 12, family = "Times New Roman") # 图例标题
  )

dev.off()



