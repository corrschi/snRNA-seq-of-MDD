library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
ptm = Sys.time()
options(stringsAsFactors = FALSE)
rm(list = ls())

###-------1.Pre-process------###
load("./mdd.filter1.Seurat.RData")

scRNA <- subset_cells
Idents(scRNA) <- 'celltype'
table(scRNA$celltype)
table(scRNA$group)

seurat.Control <- subset(scRNA, subset=group=="Control") 
seurat.Suicide <- subset(scRNA, subset=group=="Suicide") 

cellchat.Control <- createCellChat(seurat.Control@assays$RNA@data, meta = seurat.Control@meta.data, group.by = "celltype")
cellchat.Suicide <- createCellChat(seurat.Suicide@assays$RNA@data, meta = seurat.Suicide@meta.data, group.by = "celltype")


###------2.cellchat process of cllchat.Control-----####
cellchat <- cellchat.Control
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr_Control.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netp<- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway_Control.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

cellchat.Control <- cellchat
save(cellchat.Control, file="cellchat.Control.RData")

###-------3 cellchat process of cllchat.Suicide-----####
cellchat <- cellchat.Suicide
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr_Suicide.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway_Suicide.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

cellchat.Suicide <- cellchat
save(cellchat.Suicide, file="cellchat.Suicide.RData")

###-------4 mergeCellChat-----####
load("cellchat.Control.RData")
load("cellchat.Control.RData")
object.list <- list(Control = cellchat.Control, Suicide = cellchat.Suicide)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

###-------5 Plots-----####
pdf("1.Astro-Excit-Inhib interaction number.pdf",width = 8,height =6)
par(mfrow = c(1,2))
s.cell <- c("Astro","Excit_neuron", "Inhib_neuron")
count1 <- object.list[[1]]@net$count[s.cell, s.cell]
count2 <- object.list[[2]]@net$count[s.cell, s.cell]

weight.max <- max(max(count1), max(count2))
p1 <- netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                       title.name = paste0("Number of interactions-", names(object.list)[1]))
p2 <- netVisual_circle(count2, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                       title.name = paste0("Number of interactions-", names(object.list)[2]))
dev.off()


pdf("2.netVisual_heatmap_number_strength.pdf",width = 8,height =5)
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
p <- h1+h2
print(p)
dev.off()



gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("3.overall information flow of each signaling pathway.pdf", p, width = 10, height = 5.5)

pdf("4.CCompare_PTN_net.pdf",width = 10,height =6.5)
pathways.show <- c("PTN") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()



levels(cellchat@idents$joint)
gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:13),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Suicide", angle.x = 45, remove.isolate = T)
ggsave("5.Compare_LR_bubble-(Astr-Excit)-up.pdf", gg1, width = 6, height = 4)

gg2 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:13),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Suicide", angle.x = 45, remove.isolate = T)
ggsave("5.Compare_LR_bubble-(Astr-Excit)-down.pdf", gg2, width = 6, height = 4)





