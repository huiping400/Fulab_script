######cell chat analysis

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
options(stringsAsFactors = FALSE)

###在Seurat对象中提取所需要数据
data.input  <- CD8@assays$RNA@data
# create a dataframe consisting of the cell labels
identity = data.frame(group =CD8$cluster_name, row.names = names(CD8$cluster_name)) 
unique(identity$group)
##创建一个Cell Chat对象

#cellchat calculate splitted by adj_clusters
cellchat<-createCellChat(CD8,group.by="cluster_name")
groupsize<-as.numeric(table(cellchat@idents))
####小鼠中的CellChatDB包含2，021个经验证的分子相互作用，包括60%的自分泌/旁分泌信号相互作用、21%的细胞外基质（ECM）受体相互作用和19%的细胞-细胞接触相互作用。
####人的CellChatDB包含1，939个经验证的分子相互作用，包括61.8%的自分泌/旁分泌信号相互作用、21.7%的细胞外基质（ECM）受体相互作用和16.5%的细胞-细胞接触相互作用。
ccdb<-CellChatDB.mouse # use CellChatDB.mouse if running on mouse data(CellChatDB.human)
showDatabaseCategory(ccdb)

cellchat@DB<-ccdb
cellchat<-subsetData(cellchat)

cellchat<-identifyOverExpressedGenes(cellchat)
cellchat<-identifyOverExpressedInteractions(cellchat)
cellchat<-projectData(cellchat,PPI.mouse)
##CellChat 在概率计算中考虑每个细胞组中细胞比例的影响,设置population.size = TRUE
cellchat<-computeCommunProb(cellchat,raw.use=F,population.size=T)
cellchat<-filterCommunication(cellchat,min.cells=10)
#该数据框架由配体/受体级别的所有推断细胞通信组成
df.net<-subsetCommunication(cellchat)

cellchat<-computeCommunProbPathway(cellchat)
##设置slot.name = "netP"可以在信号通路级别访问推断的通信
df.netp<-subsetCommunication(cellchat,slot.name="netP")
###计算链接数或汇总通信概率来计算整合的细胞通信网络
cellchat<-aggregateNet(cellchat)
#可视化整合的细胞通信网络。
groupSize<-as.numeric(table(cellchat@idents))


par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
cellchat@netP$pathways   
#[1] "MHC-I"       "CCL"         "FN1"         "SPP1"        "MIF"         "APP"         "MHC-II"     
#[8] "GALECTIN"    "THY1"        "CD45"        "THBS"        "CXCL"        "ICAM"        "SEMA4"      
#[15] "CD86"        "TGFb"        "ITGAL-ITGB2" "COMPLEMENT"  "PD-L1"       "CSF"         "LCK"        
#[22] "LAMININ"     "CD80"        "JAM"         "PDL2"        "IFN-II"      "VISFATIN"    "ALCAM"      
#[29] "CD6"         "CD48"        "BST2"        "IL2"         "LAIR1"       "NKG2D"       "TNF"        
#[36] "CLEC"        "CD137"       "CD39"        "IL1"         "SEMA6"       "ANNEXIN"     "SELPLG"     
#[43] "IL16"        "APRIL"       "ICOS"        "NOTCH"       "OSM"         "CD200"       "CD96"       
#[50] "PVR"         "CD226"       "SN"          "PARs"        "BAFF"        "RANKL"       "FASLG"      
#[57] "LT"          "SEMA3"       "NRG"         "TWEAK"       "VCAM"        "CADM"        "ANGPTL"     
#[64] "NECTIN"      "PECAM1"      "PROS"        "HGF"         "CDH"         "ACTIVIN"     "CD22"       
#[71] "PDGF"        "LIGHT"       "IL12"        "TIGIT"       "HSPG"        "VISTA" 
#############以"CXCL"为例
pathways.show <- c("CCL")
#vertex.receiver = seq(1,4) # a numeric vector. 

### Circle plot
netVisual_individual(cellchat, signaling = pathways.show,   vertex.receiver = vertex.receiver,vertex.size.max =5,vertex.label.cex = 0.6)
netVisual_aggregate(cellchat, signaling = c("CCL"),  vertex.receiver = vertex.receiver,vertex.size.max =5,vertex.label.cex = 0.6)
# Chord diagram和弦图
group.cellType <- c(rep("CD8", 6)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show,slot.name = "netP", group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
# Plot the aggregated cell-cell communication network at the signaling pathway level
opar <- par(no.readonly=TRUE)

#计算每个配体受体对整体信号通路的贡献，并可视化由单个配体受体对调节的细胞通信
netAnalysis_contribution(cellchat, signaling = pathways.show)
####我们还可以可视化由单个配体受体对调节的细胞-细胞通信。我们提供一个函数extractEnrichedLR来提取给定信号通路的所有重要相互作用（L-R对）和相关信号基因。
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[2,] # show one ligand-receptor pair
# Hierarchy plot
#vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show)
# Circle plot
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")


pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:6), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
pdf("F:/singlecell/3_ZXH_IL2TIM3/0415/ccl_cxcl.pdf",height=8,width=12)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:6), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
dev.off()
#> Comparing communications on a single object
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(m_cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 5)
plotGeneExpression(cellchat, signaling = "CXCL",angle = 45,vjust = 0.5)
cellchat@meta$Cluster_rename <- factor(cellchat@meta$Cluster_rename, levels = c("CD8-01-T_nai","CD8-02-T_ex","CD8-03-T_eff","CD8-04-P_teff","CD8-05-T_exreg","CD4-01-T_nai","CD4-02-P_treg","CD4-03-Sp_treg","CD4-04-T_rm","CD4-05-T_em","CD4-06-T_cm","CD4-07-T_eff","CD4-08-eff_Treg","DP","DN"))
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = TRUE,angle = 45,vjust = 0.5)

##################计算和可视化网络中心分数
pathways.show <- c("CCL")
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
# Heatmap
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat)
# Signaling role analysis on the cell-cell communication networks of interest
netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
setwd("F:/singlecell/寻因生物/测序结果数据/A20_analysis/")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",color.heatmap = "BuPu",title = "netAnalysis_signalingRole_heatmap incoming",font.size = 5)

# Signaling role analysis on the cell-cell communication networks of interest
netAnalysis_signalingRole_heatmap(m_cellchat, signaling = c("CXCL", "CCL"))
cluster.cols = FALSE

library(NMF)
library(pkgmaker)
library(registry)
library(rngtools)
library(cluster)
#NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 15/16
#To enable shared memory capabilities, try: install.extras('NMF')
#The following objects are masked from 'package:igraph':
library(ggalluvial)
###在这里，我们运行selectK推断模式的数量
selectK(cellchat, pattern = "incoming")
#当传出模式数为 3 时，Cophenetic 和Silhouette值都开始突然下降。
nPatterns = 3

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,height=10,width=10)
# river plot
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,height=8,width=10)
netAnalysis_river(cellchat, pattern = "incoming")

# dot plot
netAnalysis_dot(m_cellchat, pattern = "outgoing")

#传入模式显示目标细胞（即信号接收器中的细胞）如何相互协调，以及它们如何与某些信号通路协调以响应传入的信号。
selectK(cellchat, pattern = "incoming")

#当传入模式的数量为 4 时，Cophenetic 值开始下降。
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,height=8,width=10)

# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")
#根据信号组的功能相似性识别信号组
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#install.packages("uwot")
library(uwot)
cellchat<-netEmbedding(cellchat,umap.method='uwot', type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
#基于结构相似性识别信号组
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
####上一句报错，换成下一句
cellchat<-netEmbedding(cellchat,umap.method='uwot', type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space

netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
##第五部分：保存cellchat对象
saveRDS(cellchat, file = "cellchat.rds")

