#load packages
library(stringr)
library(Seurat)
library(cowplot)
library(future)
library(Matrix)
library(ggplot2)
library(ggstatsplot)
#step1
TIM3_PROIL2<-readRDS("F:/3. singlecell/7_TIM3_PROIL2/TIM3_PROIL2_runTSNE.rds")
#### recluster for TIM3_PROIL2
TIM3_PROIL2 <- FindNeighbors(TIM3_PROIL2, reduction = "pca", dims = 1:30, nn.eps = 0.5)
TIM3_PROIL2 <- FindClusters(TIM3_PROIL2, resolution = 2.5)
TIM3_PROIL2<-RunUMAP(TIM3_PROIL2,reduction="pca",dims=1:30)
#head(TIM3_PROIL2@reductions$umap@cell.embeddings) # 提取UMAP坐标值
pdf("F:/3. singlecell/7_TIM3_PROIL2/0415/TIM3_PROIL2_res2.5_umap.pdf",height=10,width=12)
DimPlot(TIM3_PROIL2, reduction = "umap",group.by = "RNA_snn_res.2.5",label = T)
dev.off()
TIM3_PROIL2<-RunTSNE(TIM3_PROIL2,reduction="pca",dims=1:30)
#head(CD3@reductions$tsne@cell.embeddings)
pdf("F:/3. singlecell/7_TIM3_PROIL2/0415/CD3_res1.5_tsne.pdf",height=10,width=12)
DimPlot(TIM3_PROIL2, reduction = "tsne",group.by = "RNA_snn_res.2.5",label = T)
dev.off()

VlnPlot(TIM3_PROIL2,features="Cd3e",group.by="RNA_snn_res.2.5")
####step2   Extract Cd3e data(res2.5 cluster 6, 15, 24, 31, 37, 40, 45)
Idents(TIM3_PROIL2)<-TIM3_PROIL2$RNA_snn_res.2.5
CD3<-subset(TIM3_PROIL2, RNA_snn_res.2.5=="6"| RNA_snn_res.2.5=="15"|RNA_snn_res.2.5=="24"|RNA_snn_res.2.5=="31"|RNA_snn_res.2.5=="37"|RNA_snn_res.2.5=="40"|RNA_snn_res.2.5=="45")
#### recluster for CD3
CD3 <- FindNeighbors(CD3, reduction = "pca", dims = 1:30, nn.eps = 0.5)
CD3 <- FindClusters(CD3, resolution = 1.5)
CD3<-RunUMAP(CD3,reduction="pca",dims=1:30)
pdf("F:/3. singlecell/7_TIM3_PROIL2/0415/CD3_res1.5_umap.pdf",height=10,width=12)
DimPlot(CD3, reduction = "umap",group.by = "RNA_snn_res.1.5",label = T)
dev.off()
CD3<-RunTSNE(CD3,reduction="pca",dims=1:30)
pdf("F:/3. singlecell/7_TIM3_PROIL2/0415/CD3_res1.5_tsne.pdf",height=10,width=12)
DimPlot(CD3, reduction = "tsne",group.by = "RNA_snn_res.1.5",label = T)
dev.off()
saveRDS(CD3,"F:/3. singlecell/7_TIM3_PROIL2/mousedata_res2.5_CD31.5.rds")

#step4 Cd8b1 data 
VlnPlot(CD3_O,features="Cd8b1",group.by="RNA_snn_res.1.5")
FeaturePlot(object = CD3, reduction = "tsne",features = c("Cd8b1"),label = T)


#step5 Extract Cd8b1 data and recluster
Idents(CD3)<-CD3$RNA_snn_res.1.5
CD8<-subset(CD3, RNA_snn_res.1.5=="1"| RNA_snn_res.1.5=="2"|RNA_snn_res.1.5=="3"|RNA_snn_res.1.5=="5"|RNA_snn_res.1.5=="6"|RNA_snn_res.1.5=="7"|RNA_snn_res.1.5=="12"|RNA_snn_res.1.5=="13"|RNA_snn_res.1.5=="14"|RNA_snn_res.1.5=="15"|RNA_snn_res.1.5=="16"|RNA_snn_res.1.5=="18"|RNA_snn_res.1.5=="19"|RNA_snn_res.1.5=="20")
saveRDS(CD8,"F:/3. singlecell/7_TIM3_PROIL2/0415/Tim3IL2_CD8.rds")
saveRDS(CD8,"F:/3. singlecell/7_TIM3_PROIL2/0415/Tim3IL2_CD8.rds")
