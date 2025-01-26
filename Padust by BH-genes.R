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

####-------------2.wilcox test and Padj by BH method--------#####
p_values.BH <- compare_means(percent ~ group, 
                             data = dat_plot, 
                             group.by = "variable", 
                             method = "wilcox.test",
                             p.adjust.method = "BH")  # 选择 FDR 校正方法，默认是holm


print(p_values.BH)
write.csv(p_values.BH,file='Fig1g.p_values.BH.csv')