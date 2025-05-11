###########################
CD3<-readRDS("F:/singlecell/IL2TIM3/CD3.rds")

Idents(CD8)<-CD8$GROUP
group.average<-AverageExpression(CD8,group.by = "GROUP")
head(group.average[["RNA"]][,1:3])
#将计算结果返回seurat用于以下步骤计算【将空格换成_，不然ggplot2会fail】
orig.levels<-levels(CD8)
Idents(CD8)<-gsub(pattern = " ",replacement = "_",x=Idents(CD8))
orig.levels<-gsub(pattern = " ",replacement = "_",x=orig.levels)
levels(CD8)<-orig.levels
####levels 和Idents的cluster名称里的空格都要改
group.average<-AverageExpression(CD8,return.seurat = TRUE)#计算结果存在cluster.average[['RNA']]@data中，是matrix
group.average
CellScatter(group.average,cell1 = "PBS",cell2 = "TIM")
#添加样本作为细胞身份
cluster.average<-AverageExpression(CD8,return.seurat = TRUE,group.by = "cluster_name")

###1.inhibitory receptor
genes<-c("Cd160","Lag3","Cd244","Pdcd1","Tigit","Havcr2")
###2.effector genes
genes<-c("Gzma","Gzmb","Ifng")

###cytokines and cytokine receptor
genes<-c("Il2a","Ifngr1","Il7r","Il2rg","Il12rb2","Il18rap","Il1r11","Il2rb","Il12rb1","Il6st","Il10ra","Il21r","Il21")

###Costimulatory molecules
genes<-c("Tnfrsf9","Icos","Cd28","Cd226","Cd27","Cd7","Tnfrsf4")

###Transcription factor
genes<-c("Nr4a2","Foxo1","Satb1","Tbx21","Ikzf2","Irf4","Klf2","Eomes","Batf","Egr2","Maf","Bach2","Nr4a1","Tcf7","Foxo3","Lef1","Nfatc1")

###Migration and adhesion
genes<-c("Cxcr3","Itga1","Itgb7","Ly6c2","Itga4","Cd44","Ccr2","Cx3cr1","S1pr1","Cxcr5","Itgb1")

library(ggplot2)
DoHeatmap(group.average,features = genes,size = 3,draw.lines = FALSE)+scale_fill_gradientn(colors=c("navy","white","firebrick3"))

CD3_meta <- CD3@meta.data %>% select(orig.ident) %>% rownames_to_column("barcodes")
DoHeatmap(CD3_2,features = genes,size = 3,draw.lines = FALSE)+scale_fill_gradientn(colors=c("navy","white","firebrick3"))

##################################
###all cluster差异表达基因
CD8<-readRDS("F:/singlecell/3_ZXH_IL2TIM3/0415/CD8_rename.rds")
library(dplyr)
Idents(CD8)<-CD8$cluster_name
all_markers<-FindAllMarkers(CD8,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
write.csv(all_markers,"F:/singlecell/3_ZXH_IL2TIM3/0415/IL2_Tim3_markers.csv")

top10 <- all_markers %>% group_by(cluster)%>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,"F:/singlecell/3_ZXH_IL2TIM3/0415/IL2_Tim3_markers_top10.csv")

top10_sub<-c("Lyz2","Fth1","Ctsl","Fn1","Arg1","Pid1","Ccl6","Spp1","Apoe","Mmp12","Top2a","Stmn1","Pclaf","Tuba1b","Mki67","Tubb5","Birc5",
             "Ctss","Cd74","Fth1","Tyrobp","Ifitm3","Cebpb","Apoe","Hmgb2","Pclaf","Dut","Gzmb","Gzmc","Gzma",
             "Gzmf","G0s2","S100a9","Hdc","Cstdc4","S100a8","Acod1","Cxcl3","Cxcl2","Ccl3","Il1b")

pdf("F:/singlecell/3_ZXH_IL2TIM3/0415/top10_DEGs.pdf",height=10,width=12)
CD8$cluster_name<-factor(CD8$cluster_name,levels=c("naive like","stem like","early activated","effector","exhausted","central memory")) 
DoHeatmap(subset(CD8,downsample=100), features = top10_sub,size = 3) + NoLegend()
dev.off()

###################DEGs in different tissues
Idents(CD8)<-CD8$GROUP
TIM_PBS<-subset(CD8,GROUP=="PBS" | GROUP=="TIM")
TIM_PRO<-subset(CD8,GROUP=="PBS" | GROUP=="PRO")

TIM_PBS_marker<-FindMarkers(TIM_PBS,ident.1="TIM",ident.2="PBS")

TIM_PRO_marker<-FindMarkers(TIM_PRO,ident.1="PRO",ident.2="PBS")

genes<-c("Gzmf","Gzma","Gzmb","Gzmc","Gzmd","Jund","Ifitm1","Ctla2a","Dusp2","Taf10","Cmss1","Plac8","Gzme","AY036118","Ly6c2","Klrd1","Trbc2","Ms4a4b","Tmem160")
genes<-c("Gzmf","Gzma","Gzmb","Gzmc","Gzmd","Jund")
central_memory<-c("G0s2","S100a9","Hdc","Cstdc4","S100a8","Acod1","Cxcl3","Cxcl2","Ccl3,Il1b")


VlnPlot(CD8, features = "Cmss1", group.by = "GROUP",pt.size = 0)
#
FeaturePlot(TIM_PBS,features=genes,split.by="GROUP")
