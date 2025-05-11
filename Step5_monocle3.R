library(Seurat)
library(dplyr)
library(monocle3)


data<-GetAssayData(CD8,assay = 'RNA',slot = 'counts')
cell_metadata<-CD8@meta.data
gene_annotation<-data.frame(gene_short_name=rownames(data))
rownames(gene_annotation)<-rownames(data)
cds<-new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)

#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds<-preprocess_cds(cds,num_dim = 50)
plot_pc_variance_explained(cds)
###umap降维
cds<-reduce_dimension(cds,preprocess_method = "PCA")
plot_cells(cds)
###color_cells_by 参数设置umap图的颜色，可以是colData(cds)中的任何一列
colnames(colData(cds))
##以之前的Seurat分群来添加颜色，和原有的Seurat分群对比
p1 <- plot_cells(cds,reduction_method = "UMAP",color_cells_by = "cluster_name",graph_label_size=2, group_label_size=4,cell_size=0.8)+ggtitle('cds.umap')
p1
###从seurat导入整合过的umap坐标
cds.embed<-cds@int_colData$reducedDims$UMAP
int.embed<-Embeddings(CD8,reduction = "umap")
int.embed<-int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP<-int.embed
p2 <- plot_cells(cds,reduction_method = "UMAP",color_cells_by = "cluster_name",graph_label_size=2, group_label_size=4,cell_size=0.8)+ggtitle('int.umap')
p2
p=p1|p2
ggsave("F:/singlecell/3_ZXH_IL2TIM3/0415/reduction_compare.pdf",plot = p,width = 20,height = 10)
##如果细胞数目特别多（>10000细胞或更多），可以设置一些参数来加快UMAP运行速度。在reduce_dimension()
#函数中设置umap.fast_sgd=TRUE可以使用随机梯度下降方法（fast stochastic gradient descent method)加速运行。还可以使用cores参数设置多线程运算。
#可视化指定基因
ciliated_genes<-c("Sell","Cd44","Tcf7","Pdcd1","Havcr2","Mki67","Lag3","Ifng","Gzmb","Il2ra")
ciliated_genes<-c("Sell","Cd44","Tcf7","Pdcd1")
plot_cells(cds,genes = ciliated_genes,label_cell_groups = FALSE,show_trajectory_graph = TRUE)

#也可以使用tSNE降维
cds<-reduce_dimension(cds,reduction_method = "tSNE")
plot_cells(cds,reduction_method = "tSNE",color_cells_by ="cluster_name")

#随后也可使用Monocle3分cluster,鉴定每个cluster的marker基因并进行细胞注释等等。由于在Seurat的操作中
#已经对数据进行了注释，就不再使用Monocle3进行这些操作。
plot_cells(cds,reduction_method = "UMAP",color_cells_by = "cluster_name")
##4.cluster your cells
#这里的cluster其实是做分区，不同分区的细胞会进行单独的轨迹分析
cds<-cluster_cells(cds)
plot_cells(cds,color_cells_by = "partition")
#5构建细胞轨迹
#5.1 轨迹学习Learn the trajectory graph(使用learn_graph()函数)
#识别轨迹
cds<-learn_graph(cds)
p=plot_cells(cds,color_cells_by = "cluster_name",label_groups_by_cluster = FALSE,
             label_leaves = FALSE,label_branch_points = FALSE,graph_label_size=2, group_label_size=4,cell_size=0.8)

#上面这个图将被用于许多下游分析，比如分支分析和差异表达分析
p1=plot_cells(cds,color_cells_by = "cluster_name",label_groups_by_cluster = FALSE,
              label_leaves = TRUE,label_branch_points = TRUE,graph_label_size=4, group_label_size=4,cell_size=0.8)

#黑色的线显示的是graph的结构。数字带白色圆圈表示不同的结局，也就是叶子。数字带黑色圆圈代表分叉点，从这个点开始，细胞可以有多个结局。
#这些数字可以通过   label_leaves  和 label_branch_point参数设置。

#5.2细胞按拟时序排序
#在学习了graph 之后，我们就可以根据学习的发育轨迹（拟时序）排序细胞。
#为了对细胞进行排序我们首先需要告诉 Monocle哪里是这个过程的起始点。也就是需要指定轨迹的'roots'。
#手动选择root
#解决order_cells(cds)报错"object 'V1' not found"
#rownames(cds@principal_graph_aux[["UMAP"]]$dp_mst)<-NULL
#colnames(cds@int_colData@listData$reducedDims@listData$UMAP)<-NULL
cds<-order_cells(cds)
p2<-plot_cells(cds,color_cells_by = "pseudotime",label_groups_by_cluster = FALSE,
               label_leaves = FALSE,label_branch_points = FALSE,graph_label_size=4, group_label_size=4,cell_size=0.8)


#6.差异表达分析
#There are two approaches for differential analysis in Monocle:
#(1)Regression analysis:using fit_models(),you can evaluate whether each gene depends on variables such as time,treatment, etc.
#(2)Graph-autocorrelation analysis: using graph_test(), you can find genes that vary over a trajectory or between clusters.

#6.1寻找拟时轨迹差异基因  
#graph_test分析最重要的结果是莫兰指数（morans_I）,其值在-1至1之间，0代表此基因没有空间共表达效应，
#1代表此基因在空间距离相近的细胞中表达值高度相似
Track_genes<-graph_test(cds,neighbor_graph = "principal_graph",cores = 6)
Track_genes<-Track_genes[,c(5,2,3,4,1,6)]%>% filter(q_value < 1e-3)
write.csv(Track_genes,"F:/singlecell/3_ZXH_IL2TIM3/0415/Trajectory_genes.csv",row.names = F)

#6.2 挑选top10画图展示
Track_genes_sig<-Track_genes %>% top_n(n=10,morans_I) %>% pull(gene_short_name) %>% as.character()

#基因表达趋势图
p3<-plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by = "cluster_name",min_expr = 0.5,ncol = 2)


#FeaturePlot图
p4<-plot_cells(cds,genes = Track_genes_sig,show_trajectory_graph = FALSE,label_cell_groups = FALSE,label_leaves = FALSE)
p4
p4$facet$params$ncol<-5

#寻找共表达基因模块
Track_genes<-read.csv("F:/singlecell/3_ZXH_IL2TIM3/0415/Trajectory_genes.csv")
genelist<-pull(Track_genes,gene_short_name) %>% as.character()
gene_module<-find_gene_modules(cds[genelist,],resolution = 1e-1,cores = 6)
write.csv(gene_module,"F:/singlecell/3_ZXH_IL2TIM3/0415/Genes_Modeles.csv",row.names = F)
cell_group<-tibble::tibble(cell=row.names(colData(cds)),cell_group=colData(cds)$cluster_name)
agg_mat<-aggregate_gene_expression(cds,gene_module,cell_group)
row.names(agg_mat)<-stringr::str_c("Module",row.names(agg_mat))
p5<-pheatmap::pheatmap(agg_mat,scale = "column",clustering_method = "ward.D2")
p5

#提取拟时分析结果返回seurat对象
pseudotime<-pseudotime(cds,reduction_method = 'UMAP')
pseudotime<-pseudotime[rownames(CD8@meta.data)]
CD8$pseudotime<-pseudotime
p6=FeaturePlot(CD8,reduction = "umap",features = "pseudotime")
p6
saveRDS(CD8,file = "F:/singlecell/3_ZXH_IL2TIM3/0415/sco_pseudotime.rds")

##################以下3D plot 暂不完全显示
get_earliest_principal_node  <- function(cds, time_bin="central memory"){
  cell_ids <- which(colData(cds)[, "cluster_name"] == time_bin)
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

#当然也可以选择分支
cds_sub <- choose_graph_segments(cds)
#还可以做酷炫3D图
cds_3d <- reduce_dimension(cds, preprocess_method = 'PCA',max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="cluster_name")
cds_3d_plot_obj
