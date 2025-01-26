setwd("D:\\workspace\\Weekly_working\\202409week1")
library(ggpubr) 

load("./barplot_data.RData")

p <- ggbarplot(barplot_data, x="Pathway", y="-LogP", 
               fill = "#FFB6C1", ##FFB6C1粉色  #"orchid"#紫色
               color = "white", 
               orientation = "horiz",   
               palette = "nejm",  
               # palette.pals("purple"),
               legend = "right",  
               sort.val = "asc",   
               sort.by.groups=F)+   
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0))
ggsave("Fig8A.Enriched_pathways.pdf", p, width = 6, height = 3)

