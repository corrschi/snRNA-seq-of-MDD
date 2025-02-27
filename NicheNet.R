library(nichenetr) 
library(Seurat)
library(SeuratObject)
library(tidyverse)
rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("./nichenetr")
####--------1. load data--------######
load("/mdd.filter1.Seurat.RData")
seuratObj <- mdd.filter1

# plot Umap check data
p <- DimPlot(seuratObj,pt.size = 0.1,label=T,label.size=4)
ggsave("1.seuratObj-UMAP-checkdata.pdf", p, width = 7, height = 5)

####---------2.Read in NicheNet’s networks----------#######
organism = "human"
ref.path ='./Data_Reference/nichenetr/human/'

ligand_target_matrix <- readRDS(paste0(ref.path,"ligand_target_matrix_nsga2r_final.rds"))
lr_network <- readRDS(paste0(ref.path,"lr_network_human_21122021.rds"))
weighted_networks <- readRDS(paste0(ref.path,"weighted_networks_nsga2r_final.rds"))

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)
ligand_target_matrix[1:5,1:5]
head(weighted_networks$lr_sig)


####---------3.Perform the NicheNet analysis----------#############
# 1. Define a set of potential ligands for both the sender-agnostic and sender-focused approach
#set receiver
receiver = "Excit_neuron"
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

#set sender
sender_celltypes <- c("Astro")

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# 2. Define the gene set of interest
###tips: we focus on suicide vs control down-regulated ligands, nichnetr can only calculate up-regulated, suicide vs control down-regulated ligands = control vs suicide up-regulated ligands

condition_oi <-  "Control"
condition_reference <- "Suicide" 

seurat_obj_receiver <- subset(seuratObj, idents = receiver)

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "group",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val <= 0.05 & abs(avg_log2FC) >= 0.01) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

seurat_obj_receiver <- subset(seuratObj, idents = receiver)

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "group",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val <= 0.05 & abs(avg_log2FC) >= 0.01)  %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

# 3. Define the background genes
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
length(geneset_oi)


# 4. Perform NicheNet ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities
ligand_activities %>% write.csv('1.ligand_activities.csv')
# 
# plot p_hist_lig_activity
p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) +
  geom_histogram(color="black", fill="darkorange")  +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) +
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
ggsave("3.p_hist_lig_activity.pdf", p_hist_lig_activity, width = 8, height = 5)


select_list <- c( "GAST", "SEMA4B", "JAML", "SEMA4A", "SEMA5B", "SEMA5A", "JAM3", "L1CAM", "SEMA7A", "CLDN11", "NRG4", "PTN",
                  "GAS6", "SEMA4C", "SEMA4D", "JAM2", "SEMA4G", "PSAP", "PTPRM") 
best_upstream_ligands <- ligand_activities %>% filter(test_ligand %in% select_list) %>% 
  top_n(length(select_list), aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

# Plot ligand activity measure (AUPR) of these top-ranked ligands:
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

p <- (make_heatmap_ggplot(vis_ligand_aupr,
                          "Prioritized ligands", "Ligand activity", 
                          legend_title = "AUPR", color = "darkorange") + 
        theme(axis.text.x.top = element_blank()))  
ggsave("4.vis_ligand_aupr.pdf", p, width = 6, height = 12)


# 5. Infer target genes and receptors of top-ranked ligands
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 1000) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links_df%>%filter(ligand=='PTN')

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

###Plot
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.1) #多次测试cutoff

nrow(active_ligand_target_links)
## 0.25,31
head(active_ligand_target_links)
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                         color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")
ggsave("5.vis_ligand_target_cutoff0.1.pdf", p, width = 9, height = 6)

# Receptors of top-ranked ligands
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p <-(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                         y_name = "Ligands", x_name = "Receptors",  
                         color = "mediumvioletred", legend_title = "Prior interaction potential"))
ggsave("6.vis_ligand_receptor_network.pdf", p, width = 9, height = 6)

# 6. Sender-focused approach
ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% filter(test_ligand %in% select_list) %>% 
  top_n(length(select_list), aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

ggsave("7.p_ligand_aupr.pdf", p_ligand_aupr, width = 4, height = 8)

# Target gene plot
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 1000) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.1) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

ggsave("8.p_ligand_target.pdf", p_ligand_target, width = 9, height = 6)

# Receptor plot
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

ggsave("9.p_ligand_receptor.pdf", p_ligand_receptor, width = 9, height = 6)

# Visualizing expression and log-fold change in sender cells
# Dotplot of sender-focused approach
p_dotplot <- DotPlot(subset(seuratObj, celltype %in% sender_celltypes),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")
ggsave("10.p_dotplot.pdf", p_dotplot, width = 6, height = 6)


celltype_order <- levels(Idents(seuratObj)) 

# Use this if cell type labels are the identities of your Seurat object
# if not: indicate the celltype_col properly
DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype, 
  seurat_obj = seuratObj,
  condition_colname = "group",
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  celltype_col = "celltype",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                        "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", mid_color = "white",
                                        mid = median(vis_ligand_lfc), high_color = "red",
                                        legend_title = "LFC")
ggsave("11.p_lfc_DEGs.pdf", p_lfc, width = 6, height = 6)

# 7. Summary visualizations of the NicheNet analysis
figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_dotplot + theme(legend.position = "none",
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 12),
                    axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
    ylab("Expression in Sender"),
  p_lfc + theme(legend.position = "none",
                axis.title.y = element_blank()),
  p_ligand_target + theme(legend.position = "none",
                          axis.title.y = element_blank()),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc)+7, ncol(vis_ligand_lfc)+8, ncol(vis_ligand_target)))

legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
ggsave("12.combined_plot_Combine.pdf", combined_plot, width = 24, height = 7)


