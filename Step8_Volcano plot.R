rm(list=ls())
library(ggplot2)
setwd("F:/singlecell/3_ZXH_IL2TIM3/0415")


dataset <- read.table('TIM_PBS_markers.csv',sep = ",",header = TRUE)

# 设置pvalue和logFC的阈值
cut_off_pvalue = 0.0000001
cut_off_logFC = 1
dataset$change = ifelse(dataset$p_val < cut_off_pvalue & abs(dataset$avg_log2FC) >= cut_off_logFC, 
                        ifelse(dataset$avg_log2FC> cut_off_logFC ,'Up','Down'),
                        'Stable')
dataset$Label = ""
# 对差异表达基因的log2FC值进行从大到小排序取top 20
dataset <- dataset[order(abs(dataset$avg_log2FC),decreasing = T), ]
log2FC.genes <- head(dataset$X, 20)
head(log2FC.genes)
dataset <- dataset[order(abs(dataset$p_val),decreasing = T), ]
fdr.genes <- head(dataset$X, 20)

# 将log2FC.genes 和fdr.genes合并，并加入到Label
top20.genes <- c(as.character(log2FC.genes), as.character(fdr.genes))
dataset$Label[match(top20.genes, dataset$X)] <- top20.genes
#print (deg.data$Label)

# plot
ggplot(
  #设置数据
  dataset, 
  aes(x = avg_log2FC, 
      y = -log10(p_val), 
      colour=change),label = dataset$Label) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )
