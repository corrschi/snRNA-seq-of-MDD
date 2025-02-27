# Plots
library(ggpubr)
library(viridis)

load("dat_plot.RData")
####-------------1. wilcox test--------#####
png("Fig1g.Propertion of cell clusters_wilcox.test-facet.png", width = 12, height = 4, units = "in", res = 400)
ggboxplot(dat_plot,
          x = "variable", y = "percent",
          color = "group", palette = "jama",
          add = "point", outlier.colour = NULL) +
  facet_wrap(~ variable, scales = "free_x", ncol = 7) +  # 按 cluster 分面
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

####-------------2.wilcox test and Padj by BH method--------#####
p_values.BH <- compare_means(percent ~ group, 
                             data = dat_plot, 
                             group.by = "variable", 
                             method = "wilcox.test",
                             p.adjust.method = "BH") 


print(p_values.BH)
write.csv(p_values.BH,file='Fig1G.p_values.BH.csv')
