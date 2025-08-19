rm(list=ls())
library(dplyr)
library(phyloseq)
library(DESeq2)
library(edgeR)
library(ALDEx2)
#install.packages("edgeR")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("edgeR")
#library(randomForest)
library(PRROC)
library(VennDiagram)
#install.packages("RColorBrewer")
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(curatedMetagenomicData)
library(e1071)
library(glmnet)
library(parallel)
library(limma)
library(BiocParallel)
library(scales)
#lefse(version 1.1.2) wilcoxon(stats(version 4.2.1)) DESeq2(version 1.38.3) phyloseq(version 1.42.0) edgeR(version 3.42.4) ALDEx2(version 1.30.0) glmnet(version 4.1.8) PRROC(version 1.3.1) 
#############lefse#############
lefse_analysis = function(count,condition){
  lefse_count = apply(count, 1:2, function(x) as.numeric(x))
  lefse_input <- lefse_count / colSums(lefse_count)
  lefse_input = as.data.frame(lefse_input)
  condition = t(as.data.frame(condition))
  stat = data.frame(t(condition[2,]))
  colnames(stat) = condition[1,]
  # 获取第一个数据框的列名顺序
  col_order <- colnames(lefse_input)
  # 调整第二个数据框的列顺序，使其与第一个数据框保持一致
  condition_new <- stat[, col_order]
  # 合并两个数据框
  lefse_input <- rbind(condition_new, lefse_input)
  rownames(lefse_input)[1] <- "Group"
  new_row <- colnames(lefse_input)
  # 将新的一行添加到数据框的开头
  lefse_input <- rbind(new_row, lefse_input)
  # 设置行名为 "Sample"
  rownames(lefse_input)[1] <- "Sample"
  return(lefse_input)
}

#############other four methods#############
# 传入的condition都是前面是control后面是case，然后微生物数据根据样本进行调整
real_data_meta_analysis <- function(OTU_mat, condition){
  OTU_mat = apply(OTU_mat, 1:2, function(x) as.numeric(x))
  OTU_mat = as.data.frame(OTU_mat)
  # 获取 int 类型的最大值
  max_int <- .Machine$integer.max
  # 将 dataframe 中超出 int 范围的数值转换为 int 类型的最大值
  OTU_mat <- mutate_all(OTU_mat, function(x) ifelse(x > max_int, max_int, ifelse(x < -max_int, -max_int, x)))
  # 获取condition的列名顺序
  col_order <- condition$Run
  # 调整OTU_mat的列顺序，使其按control-case排列
  OTU_mat_new <- OTU_mat[,col_order]
  # 调整meta_data的行顺序，使其按control-case排列
  # meta_data_new <- meta_data[col_order,]
  
  #############Wilcoxon#############
  
  # 计算相对丰度矩阵
  OTU_mat_norm <- scale(OTU_mat,center = F,scale = colSums(OTU_mat))
  #colSums(OTU_mat_norm)#每个样本的相对丰度加和等于1
  
  condition = as.data.frame(condition)
  group_file <- data.frame(sampleID = condition$Run,
                           group = factor(condition$host_disease_stat,levels = c("control","case")))
  
  wilcoxon_result <- NULL
  ##差异分析
  for (i in 1:dim(OTU_mat_norm)[1]) {
    wilcox_test <- wilcox.test(OTU_mat_norm[i,which( colnames(OTU_mat_norm) %in% group_file[which(group_file$group=="control"),]$sampleID )],
                               OTU_mat_norm[i,which( colnames(OTU_mat_norm) %in% group_file[which(group_file$group=="case"),]$sampleID )],
                               paired= F)
    t <- c(rownames(OTU_mat_norm)[i],wilcox_test$p.value)
    wilcoxon_result <-rbind(wilcoxon_result,t)
  }
  
  colnames(wilcoxon_result) <- c("species","p.value")
  
  wilcoxon_p_value = wilcoxon_result[,2]
  wilcoxon_p_value = as.numeric(wilcoxon_p_value)
  selected_rows_wilcox <- wilcoxon_p_value < 0.05
  wilcoxon_DA <- wilcoxon_result[selected_rows_wilcox, 1]
  selected_rows_wilcox <- which(row.names(OTU_mat) %in% wilcoxon_DA)
  wilcoxon_DA_OTU <- OTU_mat[selected_rows_wilcox, ]

  # setwd("../data/PRJDB4176/not_impute/Wilcoxon")
  # write.csv(result,"before_Wilcoxon_biomarker.csv")
  
  #p<0.05 107
  
  #############DESeq2_phyloseq#############
  #需要调整微生物数据及其metadata的样本顺序一致
  #phyloseq
  #otu_table：行：OTU；列：sample
  #tax file：行：界门纲目科属种；列：OTU
  #metadata:row:sample
  # rownames(taxmat) <- rownames(otumat)
  # colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  #构造otu_table
  #行是OTU，列是样本
  OTU_mat_new_rownames = rownames(OTU_mat_new)
  #构造tax文件
  taxmat <- t(sapply(strsplit(OTU_mat_new_rownames, "/|"), function(x) {
    x <- c(x, rep("", 7 - length(x)))
    return(x)
  }))
  rownames(taxmat) <- rownames(OTU_mat_new)
  colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  # View(taxmat)
  
  #otu_table和tax_table函数构造phyloseq对象并加入phyloseq
  OTU = otu_table(OTU_mat_new, taxa_are_rows = TRUE)
  TAX = tax_table(taxmat)
  physeq = phyloseq(OTU, TAX)
  
  #构建sample_data
  # View(sample_data)
  # new_row <- c(rep("disease", 50), rep("healthy", 50))
  # sample_data <- rbind(sample_data, setNames(as.list(new_row), colnames(sample_data)))
  # rownames(sample_data)[nrow(sample_data)] <- "condition"
  
  # stat = data.frame(condition = sample_data$host_disease_stat)
  # rownames(stat) = rownames(sample_data)
  
  sampledata = sample_data(data.frame(
    condition,
    row.names=sample_names(physeq),
    stringsAsFactors=FALSE
  ))
  
  #加入phyloseq体系
  physeq = phyloseq(OTU, TAX, sampledata)
  # View(physeq)
  #第2步,使用phyloseq_to_deseq2将phyloseq对象转变成deseq2对象
  Deseq2_obj <- phyloseq_to_deseq2(physeq, ~host_disease_stat)
  #第2步,使用DESeq()对物种丰度矩阵进行差异物种分析
  #The "poscounts" estimator deals with a gene with some zeros, by calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts.
  Deseq2_dds <- DESeq(Deseq2_obj,sfType = "poscounts")
  # Deseq2_dds <- DESeq(Deseq2_obj,test="Wald", fitType="parametric")
  
  
  #第3步,使用results(),用于在DESeq分析中提取结果
  #使用results()命令提取结果时,“contrast=c(group','disease','healthy')”依次为样本分组信息表中的分组列的列名称
  #本示例中即为“group_file”中的“group”列、样本所属分组中的“disease”分组以及“healthy”分组。
  #此处“'disease”在前，“healthy”在后，即计算组间物种丰度富集下降时，是disease相较于healthy是否发生了物种丰度上调或下调过程：反之相反。
  #pAdjustMethod = "BH"
  Deseq2_res = results(Deseq2_dds, cooksCutoff = FALSE)
  alpha = 0
  sigtab = Deseq2_res[which(Deseq2_res$pvalue > alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(taxmat)[rownames(sigtab), ], "matrix"))
  Deseq2_p_value = sigtab$pvalue
  Deseq2_logFC = sigtab$log2FoldChange
  selected_rows_Deseq2 <- Deseq2_p_value < 0.05 & abs(Deseq2_logFC)>1
  Deseq2_DA <- row.names(sigtab)[selected_rows_Deseq2]
  selected_rows_Deseq2 <- which(row.names(OTU_mat_new) %in% Deseq2_DA)
  Deseq2_DA_OTU <- OTU_mat_new[selected_rows_Deseq2, ]

  #写出差异species文件
  # write.csv(sigtab, file="../data/PRJDB4176/not_impute/DESeq2_phyloseq/before_deseq2_biomarker.csv")
  
  #292个biomarker
  
  #############edgeR#############
  #一般作火山图进行展示
  #行是微生物，列是样本
  #指定分组，注意要保证表达矩阵中的样本顺序和这里的分组顺序是一一对应的
  #对照组在前，处理组在后
  group = unlist(as.vector(condition$host_disease_stat))
  
  #数据预处理
  #构建 DGEList 对象
  OTU_mat_data = as.matrix(OTU_mat_new)
  dgelist <- DGEList(counts = OTU_mat_data, group = group)
  #过滤 low count 数据，例如 CPM 标准化（推荐）
  keep <- rowSums(cpm(dgelist) > 1 ) >= 2
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]

  #标准化，以 TMM 标准化为例
  dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
  
  #差异表达基因分析
  #首先根据分组信息构建试验设计矩阵，分组信息中一定要是对照组在前，处理组在后
  design <- model.matrix(~group)
  
  #（1）估算基因表达值的离散度
  dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
  
  #（2）模型拟合，edgeR 提供了多种拟合算法
  #负二项广义对数线性模型
  fit <- glmFit(dge, design, robust = TRUE)
  lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
  edgeR_p_value = lrt[["table"]]$PValue
  edgeR_logFC = lrt[["table"]]$logFC
  selected_rows_edgeR <- edgeR_p_value < 0.05 & abs(edgeR_logFC)>1
  edgeR_DA <- row.names(lrt[["table"]])[selected_rows_edgeR]
  selected_rows_edgeR <- which(row.names(OTU_mat_new) %in% edgeR_DA)
  edgeR_DA_OTU <- OTU_mat_new[selected_rows_edgeR, ]
  
  # write.csv(lrt,"not_impute/edgeR/edgeR_result.csv")
  
  #在过滤 low count 数据时，去掉三百多个微生物
  
  #############ALDEx2#############
  # 1.用原始输入数据生成每个分类单元的后验概率分布；然后将该分布进行中心对数变换。
  # 2.将变换后的值，用参数或非参数检验进行单变量统计检验，并返回 p 值和 Benjamini-Hochberg 校正后的 p 值。
  #selex.sub:row:OTU; col:sample;
  #conds:要一个condition为向量 储存case和control
  conds = unlist(as.vector(condition$host_disease_stat))
  x.all <- aldex(OTU_mat_new, conds, mc.samples=128, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)
  sig_by_wi <- which(x.all$wi.ep > 0)
  ALDEx2_result = x.all[sig_by_wi,]
  ALDEx2_p_value = x.all$wi.ep
  selected_rows_ALDEx2 <- ALDEx2_p_value < 0.05
  ALDEx2_DA <- row.names(x.all)[selected_rows_ALDEx2]
  selected_rows_ALDEx2 <- which(row.names(OTU_mat_new) %in% ALDEx2_DA)
  ALDEx2_DA_OTU <- OTU_mat_new[selected_rows_ALDEx2, ]
  result1 <- list("wilcoxon_p_value" = wilcoxon_p_value, "Deseq2_p_value" = Deseq2_p_value, "edgeR_p_value" = edgeR_p_value, "ALDEx2_p_value" = ALDEx2_p_value)
  result2 <- list("wilcoxon_DA_OTU" = wilcoxon_DA_OTU, "Deseq2_DA_OTU" = Deseq2_DA_OTU, "edgeR_DA_OTU" = edgeR_DA_OTU, "ALDEx2_DA_OTU" = ALDEx2_DA_OTU)
  result3 <- list("wilcoxon_DA" = wilcoxon_DA, "Deseq2_DA" = Deseq2_DA, "edgeR_DA" = edgeR_DA, "ALDEx2_DA" = ALDEx2_DA)
  
  return(list(result1 = result1, result2 = result2, result3 = result3))
  # write.csv(before_ALDEx2_biomarker_wi,"../data/PRJDB4176/not_impute/ALDEx2/before_ALDEx2_biomarker_wi.csv")
}  

#############bar_plot#############
#得到的P值直接拿来做柱状
draw_barplot <- function(vec) {
  p <- ggplot(data = data.frame(x = vec), aes(x = cut(x, breaks = seq(0, 1, by = 0.05), labels = seq(0, 0.95, by = 0.05)))) +
    geom_bar(fill = "white", color = "black", size = 0.25) +
    scale_x_discrete(labels = seq(0, 0.95, by = 0.05)) +
    ylab("Frequency") +
    xlab("P Value") +
    theme(panel.background = element_rect(fill = "white"),
          axis.line = element_line(color = "black", size = 0.25),
          axis.ticks = element_line(size = 0.25))
  
  
  return(p)
}

############try_logistic_regression_with_elastic_net_penalization########
classification_results <- function(input_data, condition){
  input_data = apply(input_data, 1:2, function(x) as.numeric(x))
  input_data = as.data.frame(input_data)
  # 获取condition的列名顺序
  col_order <- condition$Run
  # 调整OTU_mat的列顺序，使其按control-case排列
  input_data_new <- input_data[,col_order]
  stat = condition$host_disease_stat
  
  input_data_new <- data.frame(t(input_data_new))
  stat <- as.factor(stat)
  
  # single result for all the taxa species.
  ## Define the K folds
  K <- 5
  set.seed(444)
  CV_split_idx <- sample(1:K, size=length(input_data_new[,1]), replace=T)
  ## logistic regression with elastic net penalization
  ### choose the lambda
  input_data_new <- as.matrix(input_data_new)
  lambda.enet <- cv.glmnet(x=input_data_new, y= as.factor(stat), alpha=0.5, family="binomial")$lambda.min
  enet_taxa_list <- list()
  elastic_net_prauc <- 0
  elastic_net_rocauc <- 0
  for(k in 1:K) {
    test_idx <- which(CV_split_idx == k)
    x_train <- input_data_new[-test_idx,]
    x_test <- input_data_new[test_idx,]
    y_train <- stat[-test_idx]
    y_test <- stat[test_idx]
    ## original sample ratio
    lr_model.orig <- glmnet(x=x_train, y=y_train, alpha=0.5, family="binomial")
    lr_pred.orig <- predict(lr_model.orig, newx=x_test, s=lambda.enet, type="response")
    lr_prauc.orig <- pr.curve(scores.class0=lr_pred.orig, weights.class0=as.numeric(as.factor(y_test))-1)$auc.integral
    elastic_net_prauc <- elastic_net_prauc + lr_prauc.orig
    enet_taxa <- coef.glmnet(lr_model.orig, s = lambda.enet)
    weights_data <- data.frame(taxa = enet_taxa@Dimnames[[1]], weight = as.numeric(enet_taxa))
    enet_taxa_list[[k]] <- weights_data
    lr_rocauc.orig <- roc.curve(scores.class0=lr_pred.orig, weights.class0=as.numeric(as.factor(y_test))-1)$auc
    elastic_net_rocauc <- elastic_net_rocauc + lr_rocauc.orig
  }
  elastic_net_prauc = elastic_net_prauc / K
  elastic_net_rocauc = elastic_net_rocauc / K
  for(i in 2 : K) {
    enet_taxa_list[[1]]$weight = enet_taxa_list[[1]]$weight + enet_taxa_list[[i]]$weight
  }
  enet_taxa_list[[1]]$weight = enet_taxa_list[[1]]$weight / K
  # 返回prauc、rocauc和权重稀疏数据框
  return(list("elastic_net.prauc" = elastic_net_prauc, "elastic_net.rocauc" = elastic_net_rocauc, "weight" = enet_taxa_list[[1]]))
}
#############PRJDB4176#############
PRJDB4176_condition = read.table("../data/PRJDB4176/PRJDB4176_condition.csv",header = T,sep = ",",check.names = F)

PRJDB4176 = read.table("../data/PRJDB4176/PRJDB4176.csv",row.names = 1,header = T ,sep = ",",check.names = F)

#未插补的数据需经过以下标准化
scale <- colSums(PRJDB4176) / 10^6 
PRJDB4176 <- floor(PRJDB4176 / scale)
PRJDB4176_condition
#lefse
PRJDB4176_lefse = lefse_analysis(PRJDB4176,PRJDB4176_condition)
write.table(PRJDB4176_lefse,"../data/PRJDB4176/lefse/PRJDB4176_lefse.txt",sep = "\t",col.names = F,quote = F)
get_PRJDB4176_pvalue = real_data_meta_analysis(PRJDB4176,PRJDB4176_condition)
write.csv(get_PRJDB4176_pvalue$result1$wilcoxon_p_value,"../data/PRJDB4176/p_distribution/wilcoxon_p_value.csv")
wilcoxon_p_value = read.csv("../data/PRJDB4176/p_distribution/wilcoxon_p_value.csv",row.names = "X")
write.csv(get_PRJDB4176_pvalue$result1$Deseq2_p_value,"../data/PRJDB4176/p_distribution/Deseq2_p_value.csv")
Deseq2_p_value = read.csv("../data/PRJDB4176/p_distribution/Deseq2_p_value.csv",row.names = "X")
write.csv(get_PRJDB4176_pvalue$result1$edgeR_p_value,"../data/PRJDB4176/p_distribution/edgeR_p_value.csv")
edgeR_p_value = read.csv("../data/PRJDB4176/p_distribution/edgeR_p_value.csv",row.names = "X")
write.csv(get_PRJDB4176_pvalue$result1$ALDEx2_p_value,"../data/PRJDB4176/p_distribution/ALDEx2_p_value.csv")
ALDEx2_p_value = read.csv("../data/PRJDB4176/p_distribution/ALDEx2_p_value.csv",row.names = "X")
LeFse_p_value = read.csv("../data/PRJDB4176/p_distribution/lefse_p_value.csv",row.names = "X")
LeFse_p_barplot = draw_barplot(LeFse_p_value)
Wilcoxon_p_barplot = draw_barplot(wilcoxon_p_value)
DESeq2_p_barplot = draw_barplot(Deseq2_p_value)
edgeR_p_barplot = draw_barplot(edgeR_p_value)
ALDEx2_p_barplot = draw_barplot(ALDEx2_p_value)
setwd("../data/PRJDB4176/p_distribution")
ggsave("before_LeFse_p_barplot.pdf", plot = LeFse_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("before_Wilcoxon_p_barplot.pdf", plot = Wilcoxon_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("before_DESeq2_p_barplot.pdf", plot = DESeq2_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("before_edgeR_p_barplot.pdf", plot = edgeR_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("before_ALDEx2_p_barplot.pdf", plot = ALDEx2_p_barplot, width = 4, height = 3, units = "in",device = "pdf")

PRJDB4176_mbSparse = read.table("../data/PRJDB4176/PRJDB4176_mbSparse.csv",row.names = 1,header = T ,sep = ",",check.names = F)

PRJDB4176_mbSparse_lefse = lefse_analysis(PRJDB4176_mbSparse,PRJDB4176_condition)
write.table(PRJDB4176_mbSparse_lefse,"../data/PRJDB4176/lefse/PRJDB4176_mbSparse_lefse.txt",sep = "\t",col.names = F,quote = F)

get_PRJDB4176_mbSparse_pvalue = real_data_meta_analysis(PRJDB4176_mbSparse,PRJDB4176_condition)
write.csv(get_PRJDB4176_mbSparse_pvalue$result1$wilcoxon_p_value,"../data/PRJDB4176/p_distribution/wilcoxon_mbSparse_p_value.csv")
wilcoxon_mbSparse_p_value = read.csv("../data/PRJDB4176/p_distribution/wilcoxon_mbSparse_p_value.csv",row.names = "X")
write.csv(get_PRJDB4176_mbSparse_pvalue$result1$Deseq2_p_value,"../data/PRJDB4176/p_distribution/Deseq2_mbSparse_p_value.csv")
Deseq2_mbSparse_p_value = read.csv("../data/PRJDB4176/p_distribution/Deseq2_mbSparse_p_value.csv",row.names = "X")
write.csv(get_PRJDB4176_mbSparse_pvalue$result1$edgeR_p_value,"../data/PRJDB4176/p_distribution/edgeR_mbSparse_p_value.csv")
edgeR_mbSparse_p_value = read.csv("../data/PRJDB4176/p_distribution/edgeR_mbSparse_p_value.csv",row.names = "X")
write.csv(get_PRJDB4176_mbSparse_pvalue$result1$ALDEx2_p_value,"../data/PRJDB4176/p_distribution/ALDEx2_mbSparse_p_value.csv")
ALDEx2_mbSparse_p_value = read.csv("../data/PRJDB4176/p_distribution/ALDEx2_mbSparse_p_value.csv",row.names = "X")
lefse_mbSparse_p_value = read.csv("../data/PRJDB4176/p_distribution/lefse_mbSparse_p_value.csv",row.names = "X")
mbSparse_LeFse_p_barplot = draw_barplot(lefse_mbSparse_p_value)
mbSparse_Wilcoxon_p_barplot = draw_barplot(wilcoxon_mbSparse_p_value)
mbSparse_DESeq2_p_barplot = draw_barplot(Deseq2_mbSparse_p_value)
mbSparse_edgeR_p_barplot = draw_barplot(edgeR_mbSparse_p_value)
mbSparse_ALDEx2_p_barplot = draw_barplot(ALDEx2_mbSparse_p_value)
setwd("../data/PRJDB4176/p_distribution")
ggsave("after_LeFse_p_barplot.pdf", plot = mbSparse_LeFse_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("after_Wilcoxon_p_barplot.pdf", plot = mbSparse_Wilcoxon_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("after_DESeq2_p_barplot.pdf", plot = mbSparse_DESeq2_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("after_edgeR_p_barplot.pdf", plot = mbSparse_edgeR_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("after_ALDEx2_p_barplot.pdf", plot = mbSparse_ALDEx2_p_barplot, width = 4, height = 3, units = "in",device = "pdf")

write.table(get_PRJDB4176_pvalue$result3$edgeR_DA,"../data/PRJDB4176/edgeR_DA.csv")
write.table(get_PRJDB4176_mbSparse_pvalue$result3$edgeR_DA,"../data/PRJDB4176/edgeR_mbSparse_DA.csv")

View

#############PRJNA397219#############
PRJNA397219_condition = read.table("../data/PRJNA397219/PRJNA397219_condtion.csv",header = T,sep = ",",check.names = F)

PRJNA397219 = read.table("../data/PRJNA397219/PRJNA397219.csv",row.names = 1,header = T ,sep = ",",check.names = F)
# 未插补的数据已经过以下标准化
scale <- colSums(PRJNA397219) / 10^6
PRJNA397219 <- floor(PRJNA397219 / scale)

#lefse
PRJNA397219_lefse = lefse_analysis(PRJNA397219,PRJNA397219_condition)
write.table(PRJNA397219_lefse,"../data/PRJNA397219/lefse/PRJNA397219_lefse.txt",sep = "\t",col.names = F,quote = F)

get_PRJNA397219_pvalue = real_data_meta_analysis(PRJNA397219,PRJNA397219_condition)
write.csv(get_PRJNA397219_pvalue$result1$wilcoxon_p_value,"../data/PRJNA397219/p_distribution/wilcoxon_p_value.csv")
wilcoxon_p_value = read.csv("../data/PRJNA397219/p_distribution/wilcoxon_p_value.csv",row.names = "X")
write.csv(get_PRJNA397219_pvalue$result1$Deseq2_p_value,"../data/PRJNA397219/p_distribution/Deseq2_p_value.csv")
Deseq2_p_value = read.csv("../data/PRJNA397219/p_distribution/Deseq2_p_value.csv",row.names = "X")
write.csv(get_PRJNA397219_pvalue$result1$edgeR_p_value,"../data/PRJNA397219/p_distribution/edgeR_p_value.csv")
edgeR_p_value = read.csv("../data/PRJNA397219/p_distribution/edgeR_p_value.csv",row.names = "X")
write.csv(get_PRJNA397219_pvalue$result1$ALDEx2_p_value,"../data/PRJNA397219/p_distribution/ALDEx2_p_value.csv")
ALDEx2_p_value = read.csv("../data/PRJNA397219/p_distribution/ALDEx2_p_value.csv",row.names = "X")
LeFse_p_value = read.csv("../data/PRJNA397219/p_distribution/lefse_p_value.csv",row.names = "X")
LeFse_p_barplot = draw_barplot(LeFse_p_value)
Wilcoxon_p_barplot = draw_barplot(wilcoxon_p_value)
DESeq2_p_barplot = draw_barplot(Deseq2_p_value)
edgeR_p_barplot = draw_barplot(edgeR_p_value)
ALDEx2_p_barplot = draw_barplot(ALDEx2_p_value)
setwd("../data/PRJNA397219/p_distribution")
ggsave("before_LeFse_p_barplot.pdf", plot = LeFse_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("before_Wilcoxon_p_barplot.pdf", plot = Wilcoxon_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("before_DESeq2_p_barplot.pdf", plot = DESeq2_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("before_edgeR_p_barplot.pdf", plot = edgeR_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("before_ALDEx2_p_barplot.pdf", plot = ALDEx2_p_barplot, width = 4, height = 3, units = "in",device = "pdf")

PRJNA397219_mbSparse = read.table("../data/PRJNA397219/PRJNA397219_mbSparse.csv",row.names = 1,header = T ,sep = ",",check.names = F)

#lefse
PRJNA397219_mbSparse_lefse = lefse_analysis(PRJNA397219_mbSparse,PRJNA397219_condition)
write.table(PRJNA397219_mbSparse_lefse,"../data/PRJNA397219/lefse/PRJNA397219_mbSparse_lefse.txt",sep = "\t",col.names = F,quote = F)

get_PRJNA397219_mbSparse_pvalue = real_data_meta_analysis(PRJNA397219_mbSparse,PRJNA397219_condition)
write.csv(get_PRJNA397219_mbSparse_pvalue$result1$wilcoxon_p_value,"../data/PRJNA397219/p_distribution/wilcoxon_mbSparse_p_value.csv")
wilcoxon_mbSparse_p_value = read.csv("../data/PRJNA397219/p_distribution/wilcoxon_mbSparse_p_value.csv",row.names = "X")
write.csv(get_PRJNA397219_mbSparse_pvalue$result1$Deseq2_p_value,"../data/PRJNA397219/p_distribution/Deseq2_mbSparse_p_value.csv")
Deseq2_mbSparse_p_value = read.csv("../data/PRJNA397219/p_distribution/Deseq2_mbSparse_p_value.csv",row.names = "X")
write.csv(get_PRJNA397219_mbSparse_pvalue$result1$edgeR_p_value,"../data/PRJNA397219/p_distribution/edgeR_mbSparse_p_value.csv")
edgeR_mbSparse_p_value = read.csv("../data/PRJNA397219/p_distribution/edgeR_mbSparse_p_value.csv",row.names = "X")
write.csv(get_PRJNA397219_mbSparse_pvalue$result1$ALDEx2_p_value,"../data/PRJNA397219/p_distribution/ALDEx2_mbSparse_p_value.csv")
ALDEx2_mbSparse_p_value = read.csv("../data/PRJNA397219/p_distribution/ALDEx2_mbSparse_p_value.csv",row.names = "X")
lefse_mbSparse_p_value = read.csv("../data/PRJNA397219/p_distribution/lefse_mbSparse_p_value.csv",row.names = "X")
mbSparse_LeFse_p_barplot = draw_barplot(lefse_mbSparse_p_value)
mbSparse_Wilcoxon_p_barplot = draw_barplot(wilcoxon_mbSparse_p_value)
mbSparse_DESeq2_p_barplot = draw_barplot(Deseq2_mbSparse_p_value)
mbSparse_edgeR_p_barplot = draw_barplot(edgeR_mbSparse_p_value)
mbSparse_ALDEx2_p_barplot = draw_barplot(ALDEx2_mbSparse_p_value)
setwd("../data/PRJNA397219/p_distribution")
ggsave("after_LeFse_p_barplot.pdf", plot = mbSparse_LeFse_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("after_Wilcoxon_p_barplot.pdf", plot = mbSparse_Wilcoxon_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("after_DESeq2_p_barplot.pdf", plot = mbSparse_DESeq2_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("after_edgeR_p_barplot.pdf", plot = mbSparse_edgeR_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("after_ALDEx2_p_barplot.pdf", plot = mbSparse_ALDEx2_p_barplot, width = 4, height = 3, units = "in",device = "pdf")

write.table(get_PRJNA397219_pvalue$result3$edgeR_DA,"../data/PRJNA397219/edgeR_DA.csv")
write.table(get_PRJNA397219_mbSparse_pvalue$result3$edgeR_DA,"../data/PRJNA397219/edgeR_mbSparse_DA.csv")

#############PRJEB7774#############
PRJEB7774_condition = read.table("../data/PRJEB7774/PRJEB7774_condition.csv",header = T,sep = ",",check.names = F)

PRJEB7774 = read.table("../data/PRJEB7774/PRJEB7774.csv",row.names = 1,header = T,sep = ",",check.names = F)
PRJEB7774 = apply(PRJEB7774, 1:2, function(x) as.numeric(x))
PRJEB7774 = as.data.frame(PRJEB7774)
#未插补的数据需经过以下标准化
scale <- colSums(PRJEB7774) / 10^6 
PRJEB7774 <- floor(PRJEB7774 / scale)

#lefse
PRJEB7774_lefse = lefse_analysis(PRJEB7774,PRJEB7774_condition)
write.table(PRJEB7774_lefse,"../data/PRJEB7774/lefse/PRJEB7774_lefse.txt",sep = "\t",col.names = F,quote = F)

get_PRJEB7774_pvalue = real_data_meta_analysis(PRJEB7774,PRJEB7774_condition)
write.csv(get_PRJEB7774_pvalue$result1$wilcoxon_p_value,"../data/PRJEB7774/p_distribution/wilcoxon_p_value.csv")
wilcoxon_p_value = read.csv("../data/PRJEB7774/p_distribution/wilcoxon_p_value.csv",row.names = "X")
write.csv(get_PRJEB7774_pvalue$result1$Deseq2_p_value,"../data/PRJEB7774/p_distribution/Deseq2_p_value.csv")
Deseq2_p_value = read.csv("../data/PRJEB7774/p_distribution/Deseq2_p_value.csv",row.names = "X")
write.csv(get_PRJEB7774_pvalue$result1$edgeR_p_value,"../data/PRJEB7774/p_distribution/edgeR_p_value.csv")
edgeR_p_value = read.csv("../data/PRJEB7774/p_distribution/edgeR_p_value.csv",row.names = "X")
write.csv(get_PRJEB7774_pvalue$result1$ALDEx2_p_value,"../data/PRJEB7774/p_distribution/ALDEx2_p_value.csv")
ALDEx2_p_value = read.csv("../data/PRJEB7774/p_distribution/ALDEx2_p_value.csv",row.names = "X")
LeFse_p_value = read.csv("../data/PRJEB7774/p_distribution/lefse_p_value.csv",row.names = "X")
LeFse_p_barplot = draw_barplot(LeFse_p_value)
Wilcoxon_p_barplot = draw_barplot(wilcoxon_p_value)
DESeq2_p_barplot = draw_barplot(Deseq2_p_value)
edgeR_p_barplot = draw_barplot(edgeR_p_value)
ALDEx2_p_barplot = draw_barplot(ALDEx2_p_value)
setwd("../data/PRJEB7774/p_distribution")
ggsave("before_LeFse_p_barplot.pdf", plot = LeFse_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("before_Wilcoxon_p_barplot.pdf", plot = Wilcoxon_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("before_DESeq2_p_barplot.pdf", plot = DESeq2_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("before_edgeR_p_barplot.pdf", plot = edgeR_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("before_ALDEx2_p_barplot.pdf", plot = ALDEx2_p_barplot, width = 4, height = 3, units = "in",device = "pdf")

PRJEB7774_mbSparse = read.table("../data/PRJEB7774/PRJEB7774_mbSparse.csv",row.names = 1,header = T ,sep = ",",check.names = F)

PRJEB7774_mbSparse_lefse = lefse_analysis(PRJEB7774_mbSparse,PRJEB7774_condition)
write.table(PRJEB7774_mbSparse_lefse,"../data/PRJEB7774/lefse/PRJEB7774_mbSparse_lefse.txt",sep = "\t",col.names = F,quote = F)

get_PRJEB7774_mbSparse_pvalue = real_data_meta_analysis(PRJEB7774_mbSparse,PRJEB7774_condition)
write.csv(get_PRJEB7774_mbSparse_pvalue$result1$wilcoxon_p_value,"../data/PRJEB7774/p_distribution/wilcoxon_mbSparse_p_value.csv")
wilcoxon_mbSparse_p_value = read.csv("../data/PRJEB7774/p_distribution/wilcoxon_mbSparse_p_value.csv",row.names = "X")
write.csv(get_PRJEB7774_mbSparse_pvalue$result1$Deseq2_p_value,"../data/PRJEB7774/p_distribution/Deseq2_mbSparse_p_value.csv")
Deseq2_mbSparse_p_value = read.csv("../data/PRJEB7774/p_distribution/Deseq2_mbSparse_p_value.csv",row.names = "X")
write.csv(get_PRJEB7774_mbSparse_pvalue$result1$edgeR_p_value,"../data/PRJEB7774/p_distribution/edgeR_mbSparse_p_value.csv")
edgeR_mbSparse_p_value = read.csv("../data/PRJEB7774/p_distribution/edgeR_mbSparse_p_value.csv",row.names = "X")
write.csv(get_PRJEB7774_mbSparse_pvalue$result1$ALDEx2_p_value,"../data/PRJEB7774/p_distribution/ALDEx2_mbSparse_p_value.csv")
ALDEx2_mbSparse_p_value = read.csv("../data/PRJEB7774/p_distribution/ALDEx2_mbSparse_p_value.csv",row.names = "X")
lefse_mbSparse_p_value = read.csv("../data/PRJEB7774/p_distribution/lefse_mbSparse_p_value.csv",row.names = "X")
mbSparse_LeFse_p_barplot = draw_barplot(lefse_mbSparse_p_value)
mbSparse_Wilcoxon_p_barplot = draw_barplot(wilcoxon_mbSparse_p_value)
mbSparse_DESeq2_p_barplot = draw_barplot(Deseq2_mbSparse_p_value)
mbSparse_edgeR_p_barplot = draw_barplot(edgeR_mbSparse_p_value)
mbSparse_ALDEx2_p_barplot = draw_barplot(ALDEx2_mbSparse_p_value)

setwd("../data/PRJEB7774/p_distribution")
ggsave("after_LeFse_p_barplot.pdf", plot = mbSparse_LeFse_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("after_Wilcoxon_p_barplot.pdf", plot = mbSparse_Wilcoxon_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("after_DESeq2_p_barplot.pdf", plot = mbSparse_DESeq2_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("after_edgeR_p_barplot.pdf", plot = mbSparse_edgeR_p_barplot, width = 4, height = 3, units = "in",device = "pdf")
ggsave("after_ALDEx2_p_barplot.pdf", plot = mbSparse_ALDEx2_p_barplot, width = 4, height = 3, units = "in",device = "pdf")

write.table(get_PRJEB7774_pvalue$result3$edgeR_DA,"../data/PRJEB7774/edgeR_DA.csv")
write.table(get_PRJEB7774_mbSparse_pvalue$result3$edgeR_DA,"../data/PRJEB7774/edgeR_mbSparse_DA.csv")

#############PRAUC#############
classify_PRJDB4176 = classification_results(get_PRJDB4176_pvalue$result2$edgeR_DA_OTU,PRJDB4176_condition)
classify_PRJNA397219 = classification_results(get_PRJNA397219_pvalue$result2$edgeR_DA_OTU,PRJNA397219_condition)
classify_PRJEB7774 = classification_results(get_PRJEB7774_pvalue$result2$edgeR_DA_OTU,PRJEB7774_condition)

classify_PRJDB4176_mbSparse = classification_results(get_PRJDB4176_mbSparse_pvalue$result2$edgeR_DA_OTU,PRJDB4176_condition)
classify_PRJNA397219_mbSparse = classification_results(get_PRJNA397219_mbSparse_pvalue$result2$edgeR_DA_OTU,PRJNA397219_condition)
classify_PRJEB7774_mbSparse = classification_results(get_PRJEB7774_mbSparse_pvalue$result2$edgeR_DA_OTU,PRJEB7774_condition)

#############韦恩图#############
#得到的微生物向量做韦恩图
#后续用AI调吧 如果这个图有需要放在文章里的话
PRJEB7774_edgeR_DA = get_PRJEB7774_pvalue[["result3"]][["edgeR_DA"]]
PRJEB7774_mbSparse_edgeR_DA = get_PRJEB7774_mbSparse_pvalue[["result3"]][["edgeR_DA"]]
P = venn.diagram(x=list(get_PRJEB7774_pvalue[["result3"]][["edgeR_DA"]],get_PRJEB7774_mbSparse_pvalue[["result3"]][["edgeR_DA"]]),
                 scaled = T, # 根据比例显示大小
                 alpha= 0.5, #透明度
                 lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
                 label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
                 cex = 2, # 数字大小
                 fontface = "bold",  # 字体粗细；加粗bold
                 fill=c("#156077","#f46f20"), # 填充色 
                 category.names = c("before", "After") , #标签名
                 cat.dist = 0.02, # 标签距离圆圈的远近
                 cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
                 cat.cex = 1, #标签字体大小
                 cat.fontface = "bold",  # 标签字体加粗
                 cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
                 cat.default.pos = "outer",  # 标签位置, outer内;text 外
                 main = "PRJEB7774",
                 output=TRUE,
                 filename="../data/PRJEB7774/venn_plot/venn_plot.png",# 文件保存
                 imagetype="png",  # 类型（tiff png svg）
                 resolution = 400,  # 分辨率
                 compression = "lzw"# 压缩算法
                 
)  
intersection <- intersect(get_PRJEB7774_pvalue[["result3"]][["edgeR_DA"]], get_PRJEB7774_mbSparse_pvalue[["result3"]][["edgeR_DA"]])
difference_before <- setdiff(get_PRJEB7774_pvalue[["result3"]][["edgeR_DA"]], get_PRJEB7774_mbSparse_pvalue[["result3"]][["edgeR_DA"]])
difference_after <- setdiff(get_PRJEB7774_mbSparse_pvalue[["result3"]][["edgeR_DA"]], get_PRJEB7774_pvalue[["result3"]][["edgeR_DA"]])


#############prauc柱状图############
library(ggplot2)
library(tidyr)

# 原始数据
data =  read.table("../data/prauc_result.csv",header = T,sep = ",",check.names = F)


# 将数据转为长格式
data_long <- pivot_longer(data, cols = c(value1, value2), names_to = "variable")

# 绘制柱状图
a = ggplot(data_long, aes(x = project, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylim(0.5, 1) +
  labs(x = "对象", y = "数值", fill = "") +
  theme_classic()
#没画好 先ai修一下
# 保存为PDF文件
pdf("bar_chart.pdf", width = 8, height = 6)
print(a)
dev.off()
getwd()
