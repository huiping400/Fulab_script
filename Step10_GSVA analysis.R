#########################GSVA
library(GSVA)
library(GSEABase)
library(limma)
library(gskb)
library(pheatmap)
library(Seurat)
library(installr)
library(Matrix)


data<-readRDS("F:/singlecell/3_ZXH_IL2TIM3/0415/CD8_rename.rds")
#make gene expression data
diff_gene_file<-read.csv("F:/singlecell/3_ZXH_IL2TIM3/0415/TIM_PBS_markers.csv",header=T,row.names=1)
diff_gene_file$gene<-rownames(diff_gene_file)

diff_gene<-as.vector(diff_gene_file$gene[!duplicated(diff_gene_file$gene)])
Idents(data)<-data$cluster_name
avg_list<-AverageExpression(data,assays="RNA",features=diff_gene)
avg<-avg_list[["RNA"]]
#avg<-read.csv("../tumor/average_exp_tumor_new.csv",header=T,row.names=1)
row.names(avg)<-toupper(row.names(avg))
avg<-as.matrix(avg)
#make gene set data-GO
gmt<-getGmt("/home/wangrj/wangrj/wangrj/mouse_breasttumor/BD_second/gskb_gmt_files/unfiltered/MousePath_GO_BP_gmt.gmt")
#GSVA analysis-GO
GSVA_result <- gsva(avg, gmt, min.sz=10, max.sz=500, verbose=FALSE, parallel.sz=1)
#make heatmap data gsva-GO
gsva_final<-data.frame()
num<-length(GSVA_result[1,])
for(i in 1:num){
  diff_gsva<-GSVA_result[,i] - (apply(GSVA_result[,1:num], 1, sum)-GSVA_result[,i])/(num-1)
  gsva3<-GSVA_result[order(-diff_gsva), ]
  gsva4<-gsva3[1:3,]
  gsva_final<-rbind(gsva_final,gsva4)
}
gsva_final<-unique(gsva_final)
#pheatmap(gsva_final,cluster_rows=F,cluster_cols=F,fontsize = 8,cellwidth=15,cellheight=10,angle_col=90,border_color=NA,filename="/storage/work/wangrj/MOUSE_SEPHIN1_RESULTS/data_reanalysis/GSVA_all_cluster_GO_BP_pathway.pdf")
pheatmap(gsva_final,cluster_rows=F,cluster_cols=F,fontsize = 8,cellwidth=15,cellheight=10,angle_col=90,border_color=NA,filename="/storage/work/wangrj/MOUSE_SEPHIN1_RESULTS/data_reanalysis/GSVA_all_cluster_GO_BP_pathway_by_sample.pdf")

#make gene set data-GO
gmt<-getGmt("/home/wangrj/wangrj/wangrj/mouse_breasttumor/BD_second/gskb_gmt_files/unfiltered/MousePath_GO_BP.gmt")
#GSVA analysis-GO
GSVA_result <- gsva(avg, gmt, min.sz=10, max.sz=500, verbose=FALSE, parallel.sz=1)
#make heatmap data gsva-GO
gsva_final<-data.frame()
num<-length(GSVA_result[1,])
for(i in 1:num){
  diff_gsva<-GSVA_result[,i] - (apply(GSVA_result[,1:num], 1, sum)-GSVA_result[,i])/(num-1)
  gsva3<-GSVA_result[order(-diff_gsva), ]
  gsva4<-gsva3[1:2,]
  gsva_final<-rbind(gsva_final,gsva4)
}
gsva_final<-unique(gsva_final)
pheatmap(gsva_final,cluster_rows=F,cluster_cols=F,fontsize = 8,cellwidth=15,cellheight=10,border_color=NA,filename="/home/wangrj/wangrj/wangrj/mouse_breasttumor/BD_second/gsva_related_analysis/tumor_KEGG__pathway.pdf")
