## Fig1D
inhib <- subset(MERGE,subset=Anno.chi %in% c("Inhib_c1", "Inhib_c2_VIP", "Inhib_c3_SST", "Inhib_c4_SST", 
                                           "Inhib_c5", "Inhib_c6_PVALB", "Inhib_c7_PVALB"))

Idents(inhib) <- 'Anno.chi'

png(paste0( "Fig1D.Dot-inhib.png"),    width = 5, height = 4, units = "in", res = 400)
DotPlot(inhib, 
        # split.by="response",
        features = rev(c("CCK","CALB2","VIP","SST","PVALB","SLC32A1","DLX1","DLX2", "DLX5", "DLX6"))) +  
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.95, size = 10, family = "Times New Roman"), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.95, size = 10, family = "Times New Roman"),
    legend.text = element_text(size = 10, family = "Times New Roman"),                                     
    legend.title = element_text(size = 12, family = "Times New Roman")                                   
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c('#330066', '#336699', '#66CC66', '#FFCC33')
  )
dev.off()

## Fig1E
excit <- subset(MERGE,subset=Anno.chi%in%c("Ex_c1_L2_L4", "Ex_c2_L2_L4", "Ex_c3_L5", "Ex_c4_L2_L4",
                                           "Ex_c5_L4_L6", "Ex_c6_L4_L6", "Ex_c7_L5_L6", "Ex_c8_L5", "Ex_c9_L2_L4",
                                           "Ex_c10_L5_L6", "Ex_c11_L5_L6",   "Ex_c12_L6"))

Idents(excit) <- 'Anno.chi'

png(paste0("Fig1E.Dot-excit.png"), width = 5, height = 4, units = "in", res = 400)
DotPlot(excit, 
        # split.by="response",
        features = rev(c( "CUX2","RASGRF2","PVRL2", "RORB","SULF2",  "PCP4",  "HTR2C","TOX",
                          "ETV1","RXFP1", "FOXP2","NR4A2", "SYNPR","TLE4","NTNG2"))) +  
  
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.95, size = 10, family = "Times New Roman"), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.95, size = 10, family = "Times New Roman"), 
    legend.text = element_text(size = 10, family = "Times New Roman"),                                     
    legend.title = element_text(size = 12, family = "Times New Roman")                                     
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c('#330066', '#336699', '#66CC66', '#FFCC33')
  )
dev.off()
