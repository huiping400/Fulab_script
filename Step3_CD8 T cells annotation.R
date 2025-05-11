#Cell subpopulation annotation
CD8 <- FindNeighbors(CD8, reduction = "pca", dims = 1:30, nn.eps = 0.5)
CD8 <- FindClusters(CD8, resolution = 0.8)
CD8<-RunUMAP(CD8,reduction="pca",dims=1:20)
DimPlot(CD8, reduction = "umap",group.by = "RNA_snn_res.1.5",label = T)
CD8<-RunTSNE(CD8,reduction="pca",dims=1:20)
DimPlot(CD8, reduction = "tsne",group.by = "RNA_snn_res.1.5",label = T)

CD8_markers<-FindAllMarkers(CD8,min.pct=0.25)

#step6 cell annotation
features<-c("Gzma","Gzmb","Pdcd1","Havcr2","Lag3","Mki67","Ifng")
genes<-c("Sell","Tcf7","Cd44","Pdcd1","Cd69","Il2ra","Il7r","Il15ra")
VlnPlot(CD8,features="Cd101",group.by="RNA_snn_res.1.5")

#1. Naïve like T ( CD44low CD62L + Cd25lo Cd69lo )    ##### c("Cd44","Sell","Il2ra","Cd69","Ccr7")
#2. Stem like T (PD-1+ TCF1+ )                   ##### c("Pdcd1","Tcf7")
#3. early activated                              ##### c("Cd44","Sell","Il2ra","Ifng","Gzmb","Pdcd1")
#4. Effector T cell (PD-1+ TIM3mi, Mki67 hi)     ##### c("Ifng","Gzmb","Tnf","Pdcd1","Havcr2","Mki67")
#5. exhausted cell  (PD-1hi TIM3 hi, Mki67mi）   ##### c("Pdcd1","Havcr2","Lag3","Ctla4","Tox3","Mki67")
#6.  Central memory T cell  (CD44 high CD62L +)  ##### c("Sell","Ccr7","Cd27","Cd28","Ptprc","Il7r")

Idents(CD8)<-CD8$RNA_snn_res.0.8
cluster_name<-c("early activated","naive like","exhausted","stem like","early activated","exhausted","effector","effector","naive like","early activated","central memory")
names(cluster_name) <- levels(CD8)
CD8<-RenameIdents(CD8,cluster_name)
CD8$cluster_name<-Idents(CD8)


CD8$cluster_name<-factor(CD8$cluster_name,levels=c("naive like","stem like","early activated","effector","exhausted","central memory")) 
DimPlot(CD8, reduction = "tsne", group.by = "cluster_name",cols=c("#E64B35B2","#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2"),
        label = TRUE,label.size =4,
        label.box = F,
        pt.size = 0.5)+
  theme_ggstatsplot()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank())+guides(color = guide_legend(override.aes = list(size=4), ncol=1) )
dev.off()

