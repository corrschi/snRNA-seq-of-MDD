# Fig1C
library(scRNAtoolVis)
png(paste0("Fig1C.Featureplot-marker.png"), width = 14, height = 6, res = 400, units = "in")
p1 <- FeatureCornerAxes(
  object = MERGE,
  reduction = 'umap',
  groupFacet = NULL,
  relLength = 0.5,
  relDist = 0.1,
  features = c("CX3CR1", "CLDN5", "GFAP", "PCDH15", "MBP"),
  minExp = 0, maxExp = 5,
  show.legend = FALSE 
) + 
  theme(
    text = element_text(family = "Times New Roman"),                                 
    axis.title = element_text(size = 12, family = "Times New Roman"),                 
    # axis.text = element_text(size = 10, family = "Times New Roman"),                 
    legend.text = element_text(size = 10, family = "Times New Roman"),                
    legend.title = element_text(size = 12, family = "Times New Roman"),              
    plot.title = element_text(size = 18, family = "Times New Roman",hjust = 0.5,face = "bold")    
  )

p2 <- FeatureCornerAxes(
  object = MERGE,
  reduction = 'umap',
  groupFacet = NULL,
  relLength = 0.5,
  relDist = 0.1,
  features = c("SATB2", "SLC17A7", "CAMK2A", "GAD1", "GAD2"),
  minExp = 0, maxExp = 5 
) + 
  theme(
    text = element_text(family = "Times New Roman"),                                   
    axis.title = element_text(size = 12, family = "Times New Roman"),                 
    # axis.text = element_text(size = 10, family = "Times New Roman"),                 
    legend.text = element_text(size = 10, family = "Times New Roman"),               
    legend.title = element_text(size = 12, family = "Times New Roman"),               
    plot.title = element_text(size = 18, family = "Times New Roman", hjust = 0.5,face = "bold")     
  )

cowplot::plot_grid(p1, p2, ncol = 1, align = 'hv')

dev.off()
