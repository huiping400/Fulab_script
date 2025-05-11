#1. cell number in diff sample
table(CD8$GROUP)

library(reshape2)
Count2Ratio<-function(sam){
  df<-(t(table(CD8@meta.data[c('GROUP', 'cluster_name')])))
  df<-as.data.frame(df)
  df2<-dcast(df, cluster_name~GROUP)
  write.table(df2, file='F:/singlecell/3_ZXH_IL2TIM3/0415/CD8_count.xls', sep='\t', quote=F, row.names=F)
  for(col in colnames(df2)){
    if(col=='cluster_name'){
      next
    }
    print(col)
    # all_sum = sum(df[col])
    # print(df[col])
    df2[col]<-prop.table(df2[col])
  }
  write.table(df2, file='F:/singlecell/3_ZXH_IL2TIM3/0415/CD8_ratio_6cluster.xls', sep='\t', quote=F, row.names=F)
}
df2$cluster<-factor(df2$cluster, levels=c("naive like","stem like","early activated","effector","exhausted","central memory")) 
df2$group<-factor(df2$group, levels=c("PBS","TIM")) 

ggplot(df2,aes(group,proportion,fill=cluster))+geom_bar(stat="identity",position="fill")+ggtitle(" ")+theme_bw()+theme(axis.ticks.length=unit(0.5,'cm'))+guides(fill=guide_legend(title=NULL))+scale_fill_manual(values=c("#E64B35B2","#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2"))
