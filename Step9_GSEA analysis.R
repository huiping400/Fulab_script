###################GSAE
library(tibble)
library(ggplot2)
library(cowplot)
library(Seurat)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
filelist<-dir("diff_gene_by_group/tumor_diff/")

gobp_pathway<-gmtPathways("F:/script/m5.all.v2023.1.Mm.symbols.gmt")
kegg_pathway<-gmtPathways("F:/script/mh.all.v2023.1.Mm.symbols.gmt")

library("msigdbr")
m_df_Hallmark <- msigdbr(species = "Mus musculus", category = "H")
m_df_kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
m_df_reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
m_df_biocarta <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:BIOCARTA")
m_t2g_Hallmark <- m_df_Hallmark %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()  
m_t2g_kegg <- m_df_kegg %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()  
m_t2g_reactome <- m_df_reactome %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()  
m_t2g_biocarta <- m_df_biocarta %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()  

m_t2g_pathways <- rbind.data.frame(m_t2g_Hallmark, m_t2g_kegg, m_t2g_reactome, m_t2g_biocarta)


#draw picture of each cell type in tumor sample
draw_gsva_barplot<-function(fl){
  #get ranks from foldchange data
  diffgene<-read.csv(paste0("diff_gene_by_group/tumor_diff/",fl),row.name=1)
  diffgene<-subset(diffgene,p_val_adj<0.05)
  diffgene<-subset(diffgene,p_val<0.05)
  diffgene$gene<-rownames(diffgene)
  diffgene$gene<-toupper(diffgene$gene)
  diffgene<-diffgene[!duplicated(diffgene$gene),]
  rownames(diffgene)<-diffgene$gene
  rownames(diffgene)<-toupper(rownames(diffgene))
  diff_frame<-data.frame(rownames(diffgene),diffgene$avg_logFC)
  ranks<-deframe(diff_frame)
}
#gsea analysis
fgseaRes <- fgsea(pathways = gobp_pathway, 
                  stats    = ranks,
                  nperm=1000,
                  minSize  = 3,
                  maxSize  = 500)

fgdata<-fgseaRes[order(fgseaRes$NES),]
upgene<-fgdata[(length(rownames(fgdata))-20):length(rownames(fgdata)),]
downgene<-fgdata[1:20,]
geneall<-rbind(downgene,upgene)
#geneall<-subset(fgdata,pval<=0.05)
geneall$Group<-"NA"
for(i in 1:length(geneall$ES)){
  if(geneall$NES[i]>0 & geneall$pval[i]<=0.05){
    geneall$Group[i]="Up"
  }
  else if(geneall$NES[i]<0 & geneall$pval[i]<=0.05){
    geneall$Group[i]="Down"
  }
  else if(geneall$pval[i]>0.05){
    geneall$Group[i]="Unsig"
  }
}
colnames(geneall)[1]<-"Pathway"

geneall$Group<-factor(geneall$Group,levels=c("Up","Down","Unsig"))
geneall$Pathway<-gsub("GO_BP_MM_","",geneall$Pathway)
geneall$Pathway<-factor(geneall$Pathway,levels=as.character(geneall$Pathway))
#barplot
pdfname<-paste0("F:/singlecell/3_ZXH_IL2TIM3/0415/TIM_PBS.pdf")
titlename="GSEA_analysis_TIM_PBS"
ggplot(geneall,aes(x=Pathway,y=NES,fill=Group))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("#FF0099","#6699FF","lightgrey"))+
  #scale_fill_manual(values = c("lightgrey"))+
  coord_flip()+guides(fill = "none")+theme_classic()+
  labs(title=titlename)+
  geom_text(data = subset(geneall, NES>0 & pval<=0.05),aes(x=Pathway, y= 0, label= paste0(Pathway,"  ")),hjust = "outward")+
  geom_text(data = subset(geneall, NES>0 & pval>0.05),aes(x=Pathway, y= 0, label= paste0(Pathway,"  ")),color="darkgrey",hjust = "outward")+
  geom_text(data = subset(geneall, NES<0 & pval>0.05),aes(x=Pathway, y= 0, label= paste0("  ",Pathway)),color="darkgrey",hjust = "inward")+
  geom_text(data = subset(geneall, NES<0 & pval<=0.05),aes(x=Pathway, y= 0, label= paste0("  ",Pathway)),hjust = "inward")+
  theme(plot.title=element_text(hjust=0.5,size=18),axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title=element_text(size=15))

