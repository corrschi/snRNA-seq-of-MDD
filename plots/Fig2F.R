## Plot Vlnplot of Fig2F
# 2.1 astro.ptn
library(ggplot2)
test_seuratobj <- subset(MERGE ,idents='Astro')
test_seuratobj$group%>%table 


A <- singlecell_gene_test(test_seuratobj,                     
                          genes.use = c('PTN'),
                          group.by = 'group',
                          comp = c("Control", "Suicide"))

anno_pvalue <- format(A$p_val, scientific = T, digits = 3)
anno_sig <- anno_pvalue

plots_violins <- VlnPlot(test_seuratobj,                          
                         pt.size = 0.0001,                         
                         group.by = "group",                         
                         features = c('PTN'),                          
                         ncol = 3,                          
                         log = FALSE,                        
                         combine = FALSE)

for (i in 1:length(plots_violins)) {  
  data <- plots_violins[[i]]$data  
  colnames(data)[1] <- 'gene'  
  plots_violins[[i]] <- plots_violins[[i]] +     
    theme_classic() +    
    theme(
      axis.text.x = element_text(size = 12, color = "black", family = "Times New Roman"),          
      axis.text.y = element_text(size = 10, color = "black", family = "Times New Roman"),          
      axis.title.y = element_text(size = 12, color = "black", family = "Times New Roman"),          
      axis.title.x = element_blank(),          
      plot.title = element_text(size = 14, color = "black", family = "Times New Roman", hjust = 0), 
      legend.position = 'none'
    ) +    
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +    
    scale_x_discrete(labels = c("Control", "Suicide")) +    
    geom_signif(annotations = anno_sig[i],               
                y_position = max(data$gene) + 0.5,               
                xmin = 1,                
                xmax = 2,                
                tip_length = 0)
}

p <- CombinePlots(plots_violins)
ggsave("Fig2F.astro_PTN_ligands.png", p, width = 3, height = 3)

# 2.2 Excit_neuron.PTPRZ1
test_seuratobj <- subset(MERGE ,idents='Excit_neuron')
test_seuratobj$group%>%table

A <- singlecell_gene_test(test_seuratobj,                     
                          genes.use = c('PTPRZ1','SDC3',"SDC4","NCL","ALK"),
                          group.by = 'group',
                          comp = c("Control", "Suicide"))

anno_pvalue <- format(A$p_val, scientific = T, digits = 3)
anno_sig <- anno_pvalue
plots_violins <- VlnPlot(test_seuratobj,                          
                         pt.size = 0.0001,                         
                         group.by = "group",                         
                         features =  c('PTPRZ1','SDC3',"SDC4","NCL","ALK"),                          
                         ncol = 5,                          
                         log = FALSE,                        
                         combine = FALSE)

for (i in 1:length(plots_violins)) {  
  data <- plots_violins[[i]]$data  
  colnames(data)[1] <- 'gene'  
  plots_violins[[i]] <- plots_violins[[i]] +     
    theme_classic() +    
    theme(
      axis.text.x = element_text(size = 12, color = "black", family = "Times New Roman"),          
      axis.text.y = element_text(size = 10, color = "black", family = "Times New Roman"),          
      axis.title.y = element_text(size = 12, color = "black", family = "Times New Roman"),          
      axis.title.x = element_blank(),          
      plot.title = element_text(size = 14, color = "black", family = "Times New Roman", hjust = 0),  # 设置标题字体
      legend.position = 'none'
    ) +    
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +    
    scale_x_discrete(labels = c("Control", "Suicide")) +    
    geom_signif(annotations = anno_sig[i],               
                y_position = max(data$gene) + 0.5,               
                xmin = 1,                
                xmax = 2,                
                tip_length = 0)
}

p <- CombinePlots(plots_violins,ncol=5)
ggsave("Fig2F.excit_neuron_receptors.png", p, width = 15, height = 3)

