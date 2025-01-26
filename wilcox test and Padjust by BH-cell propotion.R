#######---------1.BH of genes----------########
library(Seurat)
library(dplyr)

clusters <- levels(MERGE)

all_degs <- lapply(clusters, function(cluster) {
  subset_data <- subset(MERGE, idents = cluster)
  Idents(subset_data) <- "group"
  
  degs <- FindMarkers(
    subset_data,
    ident.1 = "Suicide",  # 组1
    ident.2 = "Control",  # 组2
    min.pct = 0.01,
    logfc.threshold = 0.01,
    test.use = "wilcox" 
  )
  
  # add pvalue adjust
  degs$p_adj_BH <- p.adjust(degs$p_val, method = "BH")
  
  degs$cluster <- cluster
  degs$gene <- rownames(degs)
  
  return(degs) 
})

final_degs <- bind_rows(all_degs)
head(final_degs)

final_degs <- final_degs[,c("cluster","gene", "p_val","avg_log2FC", "pct.1", "pct.2", 
                            "p_adj_BH")]
colnames(final_degs) <- c("cluster","gene", "p_val","avg_log2FC", "pct.1-Suicide", "pct.2-Control", 
                          "p_adj_BH")

write.csv(final_degs, "differential_expression_results_all_clusters.csv", row.names = F)