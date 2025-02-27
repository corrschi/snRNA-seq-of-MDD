# Fig1B
dat_plot <- MERGE
png(paste0("Fig1B.UMAP by group.png"), width = 6, height = 5, res = 400, units = "in")
DimPlot(dat_plot, pt.size = 0.3, cols = color_used, raster = FALSE, group = 'group') + 
  theme(
    text = element_text(family = "Times New Roman"),         # set Times New Roman
    axis.title = element_text(size = 12, family = "Times New Roman"),  
    axis.text = element_text(size = 10, family = "Times New Roman"),   
    legend.text = element_text(size = 12, family = "Times New Roman"), 
    legend.title = element_text(size = 12, family = "Times New Roman") 
  )
dev.off()


png(paste0("Fig1B.UMAP by sample.png"), width = 9, height = 5, res = 400, units = "in")
DimPlot(dat_plot, pt.size = 0.1, cols = color_used[c(1:34)], raster = FALSE, group = 'sample') + 
  theme(
    text = element_text(family = "Times New Roman"),       
    axis.title = element_text(size = 12, family = "Times New Roman"),  
    axis.text = element_text(size = 10, family = "Times New Roman"),   
    legend.text = element_text(size = 12, family = "Times New Roman"), 
    legend.title = element_text(size = 12, family = "Times New Roman") 
  )

dev.off()



