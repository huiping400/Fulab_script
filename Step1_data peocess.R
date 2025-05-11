options(future.globals.maxSize = 100000 * 1024^2)
filelist<-c("PBS1","PBS2","PBS3","TIM7","TIM8","TIM9")
namelist<-c("PBS1","PBS2","PBS3","TIM7","TIM8","TIM9")
####数据集中测到的少于200个基因的细胞（min.features = 100）和少于3个细胞覆盖的基因（min.cells = 2）被过滤掉
plan("multisession", workers = 8)
TIM3_PROIL2<-list()
for(i in 1:length(filelist)){
  fl=filelist[i]
  sname=paste0("F:/3. singlecell/7_TIM3_PROIL2/cellranger_result/",fl,"/raw_feature_bc_matrix/")
  sample.raw<-Read10X(sname)
  sample<-CreateSeuratObject(counts=sample.raw,project="TIM3_PROIL2",min.cells=2)
  sample$COMPARE<-namelist[i]
  sample<-subset(sample,subset=nFeature_RNA>100)
  sample<-NormalizeData(sample,verbose=FALSE)
  sample<-FindVariableFeatures(sample,selection.method="vst",nfeatures=4000)
  TIM3_PROIL2<-c(TIM3_PROIL2,list(sample))
}
names(TIM3_PROIL2)<-namelist

TIM3_PROIL2.anchors<-FindIntegrationAnchors(object.list=TIM3_PROIL2,dims=1:30)

TIM3_PROIL2.data<-IntegrateData(anchorset=TIM3_PROIL2.anchors,dims=1:30)
DefaultAssay(TIM3_PROIL2.data)<-"RNA"
TIM3_PROIL2.data<-FindVariableFeatures(TIM3_PROIL2.data)
#####################
# Find  mitochondrial genes, compute the mitochondrial rate
mito.genes <- grep("^mt-", rownames(TIM3_PROIL2.data), value = TRUE)
percent.mito <- Matrix::colSums(TIM3_PROIL2.data[mito.genes, ])/Matrix::colSums(TIM3_PROIL2.data)
TIM3_PROIL2.data <- AddMetaData(object = TIM3_PROIL2.data, metadata = percent.mito, col.name = "percent.mito")
# Find ribosomal genes, compute the mitochondrial rate
ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = TIM3_PROIL2.data), value = TRUE)
percent.ribo <- Matrix::colSums(TIM3_PROIL2.data[ribo.genes, ])/Matrix::colSums(TIM3_PROIL2.data)
TIM3_PROIL2.data <- AddMetaData(object = TIM3_PROIL2.data, metadata = percent.ribo, col.name = "percent.ribo")
# Calculate housekeeping gene score(UMI sum)
housekeeking_marker <- c("Actb","Gapdh","Malat1")
hw_data <- TIM3_PROIL2.data[housekeeking_marker,]
housekeeping_UMI <- colSums(hw_data)
TIM3_PROIL2.data <- AddMetaData(object = TIM3_PROIL2.data, metadata = housekeeping_UMI, col.name = "housekeeping_score")
# summary
dim(TIM3_PROIL2.data);median(TIM3_PROIL2.data$nFeature_RNA);median(TIM3_PROIL2.data$nCount_RNA);median(TIM3_PROIL2.data@meta.data$percent.mito)
table(TIM3_PROIL2.data$nFeature_RNA < 200);table(TIM3_PROIL2.data@meta.data$percent.mito > 0.2);table(TIM3_PROIL2.data@meta.data$percent.ribo > 0.5);table(TIM3_PROIL2.data@meta.data$housekeeping_score < 1)
TIM3_PROIL2.data <- subset(TIM3_PROIL2.data, nFeature_RNA >= 200|percent.mito <= 0.2|percent.ribo <= 0.5|housekeeping_score >= 1)

TIM3_PROIL2.data <- NormalizeData(object = TIM3_PROIL2.data)
TIM3_PROIL2.data <- FindVariableFeatures(object = TIM3_PROIL2.data)
TIM3_PROIL2.data <- ScaleData(object = TIM3_PROIL2.data)
TIM3_PROIL2.data<-RunPCA(TIM3_PROIL2.data)

TIM3_PROIL2.data <- FindNeighbors(TIM3_PROIL2.data, reduction = "pca", dims = 1:30, nn.eps = 0.5)
library(future)
TIM3_PROIL2.data <- FindClusters(TIM3_PROIL2.data, resolution = 0.6)
TIM3_PROIL2.data<-RunUMAP(TIM3_PROIL2.data,reduction="pca",dims=1:20)
head(TIM3_PROIL2.data@reductions$umap@cell.embeddings) 
p1 <-DimPlot(TIM3_PROIL2.data, reduction = "umap",label = T)
p1
TIM3_PROIL2.data2<-RunTSNE(TIM3_PROIL2.data,reduction="pca",dims=1:20)
head(TIM3_PROIL2.data2@reductions$tsne@cell.embeddings)
p2 <- DimPlot(TIM3_PROIL2.data2, reduction = "tsne",label = T)
p2
p1 + p2
saveRDS(TIM3_PROIL2.data,"F:/3. singlecell/7_TIM3_PROIL2/TIM3_PROIL2_runUMAP.rds")
saveRDS(TIM3_PROIL2.data2,"F:/3. singlecell/7_TIM3_PROIL2/TIM3_PROIL2_runTSNE.rds")
