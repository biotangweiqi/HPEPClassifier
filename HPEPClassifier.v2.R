#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
library(scales)

#########################
### command arguments ###
#########################
opt <- commandArgs(TRUE)
file_expr_data         <- opt[1]
file_experiment_design <- opt[2]
file_comparison_list   <- opt[3]
file_union_comparison  <- opt[4]
#
dir_DEG_list <- opt[5]
dir_output  <- opt[6]

# example
file_expr_data         <- "../04-1.exprBackground/keepGene_FPKM.xls"
file_experiment_design <- "../experiment_design.xls"
file_comparison_list   <- "../comparison_list.xls"
file_union_comparison  <- "../union_comparison.xls"

dir_DEG_list <- "02.DEGs_list2"
dir_output   <- "../06.familyDEGs/01.pedigreeExprPattern2"

# create output dir
if(!dir.exists(dir_output)) dir.create(dir_output)

#################
### functions ###
#################
## functions obtained from coolmap
cal_zScore <- function(x){
  M <- rowMeans(x, na.rm=TRUE)
  nsamples <- ncol(x)
  DF <- nsamples - 1L
  IsNA <- is.na(x)
  if(any(IsNA)) {
    mode(IsNA) <- "integer"
    DF <-  DF - rowSums(IsNA)
    DF[DF==0L] <- 1L
  }
  x <- x-M
  V <- rowSums(x^2L, na.rm=TRUE) / DF
  x <- x / sqrt(V+0.01)
  return(x)
}

get_DEG_avg <- function(DEG_expr, group_name, sample_list){
  # calculate average of samples in one group
  sample_name <- colnames(DEG_expr)
  DEG_avg <- data.frame()
  for(group in group_name){
    idx <- which(sample_list$group == group)
    samples <- sample_list$reportName[idx]
    #jdx <- which(sample_name %in% sample_list$ids[idx])
    if(length(samples) == 1) value <- DEG_expr[,samples]
    if(length(samples) > 1 ) value <- rowMeans(DEG_expr[,samples])
    #
    if(length(DEG_avg) >1) DEG_avg <- cbind(DEG_avg, value )
    if(length(DEG_avg) == 0) DEG_avg <- value
  }
  colnames(DEG_avg) <- group_name
  rownames(DEG_avg) <- rownames(DEG_expr)
  DEG_avg <- as.data.frame(DEG_avg)
  return(DEG_avg)
  
}

identify_DEGs_classMode <- function(DEG_avg, DEG_list){
  ## first compare P and M
  P_gt_M <- DEG_avg[which(DEG_avg$P > DEG_avg$M & DEG_unique_IDs %in% DEG_list$PvsM),]      # P > M and significant
  P_ge_M <- DEG_avg[which(DEG_avg$P > DEG_avg$M & !(DEG_unique_IDs %in% DEG_list$PvsM) ),]  # P >= M but not significant
  P_lt_M <- DEG_avg[which(DEG_avg$P < DEG_avg$M & DEG_unique_IDs %in% DEG_list$PvsM),]      # P < M and significant
  P_le_M <- DEG_avg[which(DEG_avg$P < DEG_avg$M & !(DEG_unique_IDs %in% DEG_list$PvsM) ),]  # P <= M but not significant
  
  ## second compare P and H, M and H
  ## P > M
  # pattern 1
  H_gt_P_gt_M <- P_gt_M[which( P_gt_M$H > P_gt_M$P &                      # compare expression: H > P
                               (rownames(P_gt_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_gt_M$H > P_gt_M$M &                      # compare expression: H > M
                               (rownames(P_gt_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH
  
  # pattern 2
  H_ge_P_gt_M <- P_gt_M[which( P_gt_M$H >= P_gt_M$P &                     # compare expression: H >= P
                               !(rownames(P_gt_M) %in% DEG_list$PvsH) &   # not existed in DEG.PH
                               P_gt_M$H > P_gt_M$M &                      # compare expression: H > M
                               (rownames(P_gt_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH

  # pattern 3
  P_ge_H_gt_M <- P_gt_M[which( P_gt_M$P >= P_gt_M$H &                     # compare expression: P >= H
                               !(rownames(P_gt_M) %in% DEG_list$PvsH) &   # not existed in DEG.PH
                               P_gt_M$H > P_gt_M$M &                      # compare expression: H > M
                               (rownames(P_gt_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH
  
  # pattern 4
  P_ge_H_ge_M <- P_gt_M[which( P_gt_M$P >= P_gt_M$H &                     # compare expression: P >= H
                               !(rownames(P_gt_M) %in% DEG_list$PvsH) &   # not existed in DEG.PH
                               P_gt_M$H >= P_gt_M$M &                     # compare expression: H >= M
                               !(rownames(P_gt_M) %in% DEG_list$MvsH) ),] # not existed in DEG.MH
  
  # pattern 5
  P_gt_H_gt_M <- P_gt_M[which( P_gt_M$P > P_gt_M$H &                      # compare expression: P > H
                               (rownames(P_gt_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_gt_M$H > P_gt_M$M &                      # compare expression: H > M
                               (rownames(P_gt_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH
  
  # pattern 6
  P_gt_H_ge_M <- P_gt_M[which( P_gt_M$P > P_gt_M$H &                      # compare expression: P > H
                               (rownames(P_gt_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_gt_M$H >= P_gt_M$M &                     # compare expression: H >= M
                               !(rownames(P_gt_M) %in% DEG_list$MvsH) ),] # not existed in DEG.MH
  
  # pattern 7
  P_gt_M_ge_H <- P_gt_M[which( P_gt_M$P > P_gt_M$H &                      # compare expression: P > H
                               (rownames(P_gt_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_gt_M$M >= P_gt_M$H &                     # compare expression: M >= H
                               !(rownames(P_gt_M) %in% DEG_list$MvsH) ),] # not existed in DEG.MH
  
  # pattern 8
  P_gt_M_gt_H <- P_gt_M[which( P_gt_M$P > P_gt_M$H &                      # compare expression: P > H
                               (rownames(P_gt_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_gt_M$M > P_gt_M$H &                      # compare expression: M > H
                               (rownames(P_gt_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH
  
  ## P >= M
  # pattern 9
  H_gt_P_ge_M <- P_ge_M[which( P_ge_M$H > P_ge_M$P &                      # compare expression: H > P
                               (rownames(P_ge_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_ge_M$H > P_ge_M$M &                      # compare expression: H > M
                               (rownames(P_ge_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH
  
  # pattern 10
  H_ge_P_ge_M <- P_ge_M[which( P_ge_M$H >= P_ge_M$P &                     # compare expression: H >= P
                               !(rownames(P_ge_M) %in% DEG_list$PvsH) &   # not existed in DEG.PH
                               P_ge_M$H > P_ge_M$M &                      # compare expression: H > M
                               (rownames(P_ge_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH
  
  # pattern 11
  P_ge_M_gt_H <- P_ge_M[which( P_ge_M$P > P_ge_M$H &                      # compare expression: P > H
                               (rownames(P_ge_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_ge_M$M > P_ge_M$H &                      # compare expression: M > H
                               (rownames(P_ge_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH

  # pattern 12
  P_ge_M_ge_H <- P_ge_M[which( P_ge_M$P > P_ge_M$H &                      # compare expression: P > H
                               (rownames(P_ge_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_ge_M$M >= P_ge_M$H &                     # compare expression: M >= H
                               !(rownames(P_ge_M) %in% DEG_list$MvsH) ),] # not existed in DEG.MH
  
  ## P <= M
  # pattern 13
  H_lt_P_le_M <- P_le_M[which( P_le_M$H < P_le_M$P &                      # compare expression: H < P
                               (rownames(P_le_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_le_M$H < P_le_M$M &                      # compare expression: H < M
                               (rownames(P_le_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH
  
  # pattern 14
  H_le_P_le_M <- P_le_M[which( P_le_M$H <= P_le_M$P &                     # compare expression: H <= P 
                               !(rownames(P_le_M) %in% DEG_list$PvsH) &   # not existed in DEG.PH
                               P_le_M$H < P_le_M$M &                      # compare expression: H < M
                               (rownames(P_le_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH
  
  # pattern 15
  P_le_M_lt_H <- P_le_M[which( P_le_M$P < P_le_M$H &                      # compare expression: P < H
                               (rownames(P_le_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_le_M$M < P_le_M$H &                      # compare expression: M < H
                               (rownames(P_le_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH
  
  # pattern 16
  P_le_M_le_H <- P_le_M[which( P_le_M$P < P_le_M$H &                      # compare expression: P < H
                               (rownames(P_le_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_le_M$M <= P_le_M$H &                     # compare expression: M <= H
                               !(rownames(P_le_M) %in% DEG_list$MvsH) ),] # not existed in DEG.MH
  
  ## P < M
  # pattern 17
  H_lt_P_lt_M <- P_lt_M[which( P_lt_M$H < P_lt_M$P &                      # compare expression: H < P
                               (rownames(P_lt_M) %in% DEG_list$PvsH) &    # existed in DEG.PH 
                               P_lt_M$H < P_lt_M$M &                      # compare expression: H < M
                               (rownames(P_lt_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH
  
  # pattern 18
  H_le_P_lt_M <- P_lt_M[which( P_lt_M$H <= P_lt_M$P &                     # compare expression: H <= P
                               !(rownames(P_lt_M) %in% DEG_list$PvsH) &   # not existed in DEG.PH  
                               P_lt_M$H < P_lt_M$M &                      # compare expression: H < M
                               (rownames(P_lt_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH

  # pattern 19  
  P_le_H_lt_M <- P_lt_M[which( P_lt_M$P <= P_lt_M$H &                     # compare expression: P <= H
                               !(rownames(P_lt_M) %in% DEG_list$PvsH) &   # not existed in DEG.PH  
                               P_lt_M$H < P_lt_M$M &                      # compare expression: H < M
                               (rownames(P_lt_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH
  
  # pattern 20
  P_le_H_le_M <- P_lt_M[which( P_lt_M$P <= P_lt_M$H &                     # compare expression: P <= H
                               !(rownames(P_lt_M) %in% DEG_list$PvsH) &   # not existed in DEG.PH  
                               P_lt_M$H <= P_lt_M$M &                     # compare expression: H <= M
                               !(rownames(P_lt_M) %in% DEG_list$MvsH) ),] # not existed in DEG.MH
  
  # pattern 21
  P_lt_H_lt_M <- P_lt_M[which( P_lt_M$P < P_lt_M$H &                      # compare expression: P < H
                               (rownames(P_lt_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_lt_M$H < P_lt_M$M &                      # compare expression: H < M
                               (rownames(P_lt_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH
  
  # pattern 22
  P_lt_H_le_M <- P_lt_M[which( P_lt_M$P < P_lt_M$H &                      # compare expression: P < H
                               (rownames(P_lt_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_lt_M$H <= P_lt_M$M &                     # compare expression: H <= M
                               !(rownames(P_lt_M) %in% DEG_list$MvsH) ),] # not existed in DEG.MH
  
  # pattern 23
  P_lt_M_le_H <- P_lt_M[which( P_lt_M$P < P_lt_M$H &                      # compare expression: P < H
                               (rownames(P_lt_M) %in% DEG_list$PvsH) &    # existed in DEG.PH
                               P_lt_M$M <= P_lt_M$H &                     # compare expression: M <= H
                               !(rownames(P_lt_M) %in% DEG_list$MvsH) ),] # not existed in DEG.MH
  
  # pattern 24
  P_lt_M_lt_H <- P_lt_M[which( P_lt_M$P < P_lt_M$H &                      # compare expression: P < H
                               (rownames(P_lt_M) %in% DEG_list$PvsH) &    # existed in DEG.PH 
                               P_lt_M$M < P_lt_M$H &                      # compare expression: M < H
                               (rownames(P_lt_M) %in% DEG_list$MvsH) ),]  # existed in DEG.MH
  
  return(list(
    # P>M, eight class-modes
    "H>P>M"   = H_gt_P_gt_M, # 1, overdominant to P
    "H>=P>M"  = H_ge_P_gt_M, # 2, deviation to P
    "P>=H>M"  = P_ge_H_gt_M, # 3, deviation to P
    "P>=H>=M" = P_ge_H_ge_M, # 4, middle
    "P>H>M"   = P_gt_H_gt_M, # 5, middle
    "P>H>=M"  = P_gt_H_ge_M, # 6, deviation to M
    "P>M>=H"  = P_gt_M_ge_H, # 7, deviation to M
    "P>M>H"   = P_gt_M_gt_H, # 8, overdominant to M
    
    # P>=M, four class-modes
    "H>P>=M"  = H_gt_P_ge_M, # 9,  overdominant to P
    "H>=P>=M" = H_ge_P_ge_M, # 10, deviation to P
    "P>=M>H"  = P_ge_M_gt_H, # 11, overdominant to M
    "P>=M>=H" = P_ge_M_ge_H, # 12, deviation to M
    
    # P<=M, four class-modes
    "H<P<=M"  = H_lt_P_le_M, # 13, overdomiant to P
    "H<=P<=M" = H_le_P_le_M, # 14, deviation to P
    "P<=M<H"  = P_le_M_lt_H, # 15, overdominant to M
    "P<=M<=H" = P_le_M_le_H, # 16, deviation to M
    
    # P<M, eight class-modes
    "H<P<M"   = H_lt_P_lt_M, # 17, overdominant to P
    "H<=P<M"  = H_le_P_lt_M, # 18, deviation to P
    "P<=H<M"  = P_le_H_lt_M, # 19, deviation to P
    "P<=H<=M" = P_le_H_le_M, # 20, middle
    "P<H<M"   = P_lt_H_lt_M, # 21, middle
    "P<H<=M"  = P_lt_H_le_M, # 22, deviation to M
    "P<M<=H"  = P_lt_M_le_H, # 23, deviation to M
    "P<M<H"   = P_lt_M_lt_H  # 24, overdominant to M
  ) )  
}

#####################
### main pipeline ###
#####################
# read info of samples
cat("read info of samples\n")
sample_list <- read.table(file_experiment_design, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
print(sample_list)

#
cat("reading information of comparison list\n")
comparison_list <- read.table(file_comparison_list, header = TRUE, stringsAsFactors = FALSE)
rownames(comparison_list) <- comparison_list$CompareID
print(comparison_list)

#
cat("reading unions of comparisons\n")
union_comparison <- read.table(file_union_comparison, header = TRUE, stringsAsFactors = FALSE)
union_list <- strsplit(union_comparison$ComparisonIDs, split = ";")
names(union_list) <- union_comparison$UnionID
print(union_list)

# read data
cat("read expr data\n")
expr_data <- read.table(file_expr_data, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_data) <- expr_data[,1]
head(expr_data)

# 
cat("select samples\n")
expr_data <- expr_data[,sample_list$ids]
colnames(expr_data) <- sample_list$reportName
head(expr_data)

###############################################
#
###############################################
classMode_countTable  <- data.frame()
union_comparison_names <- c()
for(unionID in union_comparison$UnionID){
  cat(unionID, "\n")
  #print(union_list[[unionID]])
  
  cp_ids <- union_list[[unionID]]
  #print(cp_ids)
  
  cp_list <- comparison_list[cp_ids,]
  #print(cp_list)
  
  DEG_list <- list()
  union_group <- c()
  for(i in 1:nrow(cp_list) ){
    # group 1 and group 2
    group1 <- cp_list[i,]$Group1
    group2 <- cp_list[i,]$Group2
    union_group <- c(union_group, group1, group2)
    #cat("Group 1 is", group1, "\n")
    #cat("Group 2 is", group2, "\n")
    cp_name <- paste(cp_list[i,]$CompareID, ".", group1, "_vs_", group2, sep = "")
    #cat(cp_name, "\n")
    #
    file_list <- paste(dir_DEG_list, "/", cp_name, ".all_DEGs.xls", sep = "")
    one_list <- read.table(file_list, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    #
    role1 <- unique(sample_list$pedigree[which(sample_list$group == group1)] )
    role2 <- unique(sample_list$pedigree[which(sample_list$group == group2)] )
    key <- paste(role1, "vs", role2, sep = "")
    DEG_list[[key]] <- one_list[,1]
  }
  #
  #head(DEG_list)
  if(is.null(DEG_list$PvsM) ) DEG_list$PvsM <- DEG_list$MvsP
  if(is.null(DEG_list$PvsH) ) DEG_list$PvsH <- DEG_list$HvsP
  if(is.null(DEG_list$MvsH) ) DEG_list$MvsH <- DEG_list$HvsM
  
  # a list of all DEGs
  cat("get gene list\n")
  DEG_unique_IDs <- unique(unlist(DEG_list) )
  head(DEG_unique_IDs)
  length(DEG_unique_IDs)

  #
  cat("get samples in this union")
  union_group <- unique(union_group)
  union_samples <- sample_list[which(sample_list$group %in% union_group),]
  print(union_samples)
  
  # get gene expression data
  cat("get gene expression data\n")
  DEG_expr <- expr_data[DEG_unique_IDs, union_samples$reportName]
  head(DEG_expr)
  

  
  # calculate mean of group
  cat("calculate mean of group\n")
  group_name <- unique(union_samples$group)
  DEG_avg <- get_DEG_avg(DEG_expr, group_name, sample_list)
  head(DEG_avg)

  # transform to log10(FPKM+1)
  cat("transform to log10(FPKM+1)\n")
  DEG_logavg <- log10(DEG_avg + 1)
  
  # transform to zscore
  cat("transform to zscore\n")
  #DEG_zscore <- cal_zScore(DEG_avg)
  DEG_zscore <- cal_zScore(DEG_logavg)
  #head(DEG_zscore)
  
  # change name
  cat("change name\n")
  colnames(DEG_avg) <- unique(union_samples$pedigree)
  head(DEG_avg)
  
  # union comparison names
  group_ids <- data.frame(group = unique(union_samples$group) )
  rownames(group_ids) <- unique(union_samples$pedigree)
  #P_tag <- as.character(group_ids[c("P"),])
  #M_tag <- as.character(group_ids[c("M"),])
  #union_comparison_names <- c(union_comparison_names, paste(P_tag, "-x-", M_tag, sep = "") )
  H_tag <- as.character(group_ids[c("H"),]) 
  union_comparison_names <- c(union_comparison_names, H_tag )
  
  # class-mode
  cat("get class modes\n")
  classMode_list <- identify_DEGs_classMode(DEG_avg, DEG_list)

  #
  #
  PvsM <- DEG_unique_IDs %in% DEG_list$PvsM
  PvsH <- DEG_unique_IDs %in% DEG_list$PvsH
  MvsH <- DEG_unique_IDs %in% DEG_list$MvsH
  colnames(DEG_expr) <- paste(union_samples$pedigree, union_samples$reportName, sep = ".")
  union_DEG_expr <- data.frame(id = DEG_unique_IDs, DEG_expr, PvsM, PvsH, MvsH)
  #
  exprPattern <- rep(NA, length(DEG_unique_IDs))
  for(mode_name in names(classMode_list) ){
    exprPattern <- ifelse(DEG_unique_IDs %in% rownames(classMode_list[[mode_name]]), mode_name, exprPattern)
  }
  union_DEG_expr <- cbind(union_DEG_expr, exprPattern)
  #

  file_union_DEG_expr <- paste(dir_output, "/", unionID, ".", H_tag, ".exprPattern.xls", sep = "")
  write.table(union_DEG_expr, file_union_DEG_expr, col.names = T, row.names = F, sep = "\t")  
    
  # create dir of subclass
  cat("create dir of subclass\n")
  dir_each_union_subclass <- paste(dir_output, "/", unionID, ".", make.names(H_tag), ".subclass", sep = "")
  if(!dir.exists(dir_each_union_subclass)) dir.create(dir_each_union_subclass, recursive = T)

  # write out data each class
  cat("write out data each class\n")
  classMode_name  <- names(classMode_list)
  classMode_count <- c()
  for(i in 1:length(classMode_list) ){
    data <- as.data.frame(classMode_list[[i]] )
    cat(i, ".", classMode_name[i], "\t", nrow(data), "\n", sep = "")
    classMode_count <- c(classMode_count, nrow(data))
    #
    if(nrow(data) > 0){
      prefix_each_class <- paste(dir_each_union_subclass, "/subclass_", i, sep = "")
      file_each_class_data <- paste(prefix_each_class, ".DEG_avg.xls", sep = "")
      colnames(data) <- group_name
      write.table(cbind(id = rownames(data), data), file_each_class_data, sep = "\t", quote = F, row.names = F)
      #
      file_each_class_zscore <- paste(prefix_each_class, ".DEG_zscore.xls", sep = "")
      write.table(cbind(id = rownames(data), DEG_zscore[rownames(data),]), file_each_class_zscore, sep = "\t", quote = F, row.names = F)
      #
      file_each_class_expr <- paste(prefix_each_class, ".DEG_expr.xls", sep = "")
      write.table(cbind(id = rownames(data), DEG_expr[rownames(data),]), file_each_class_expr, sep = "\t", quote = F, row.names = F)
      #
      #file_linechart_tiff <- paste(prefix_each_class, ".ZScore.linechart.tiff", sep = "")
      #title_linechart <- paste("subclass_", i, ". ", classMode_name[i], ": ", nrow(data), " genes", sep = "")
      #draw_ZScore_linechart(DEG_zscore[rownames(data),], file_linechart_tiff, title_linechart)
    }
  }
  #
  if(nrow(classMode_countTable) > 0 ) classMode_countTable <- cbind(classMode_countTable, classMode_count)
  if(nrow(classMode_countTable) == 0 ) classMode_countTable <- data.frame(classMode_count)
  rownames(classMode_countTable) <- paste(1:length(classMode_name), ". ", classMode_name, sep = "")
}

#colnames(classMode_countTable) <- union_comparison$UnionID 
colnames(classMode_countTable) <- union_comparison_names
print(classMode_countTable)

file_countTable <- paste(dir_output, "/classMode_counts.summary.xls", sep = "")
write.table(cbind("Expression Pattern" = rownames(classMode_countTable), classMode_countTable),
            file = file_countTable, sep = "\t", row.names = F, quote = F)

classMode_countMelt <- melt(as.matrix(classMode_countTable) )
colnames(classMode_countMelt) <- c("pattern", "pedigree", "number")

file_tiff_countTable <- paste(dir_output, "/classMode_counts.bar.tiff", sep = "")
p <- ggplot(data = classMode_countMelt, aes(x = pattern, y = number, fill = pedigree) )
p <- p + geom_bar(stat = "identity")
p <- p + xlab("Expression pattern") + theme(axis.text.x = element_text(angle = -45, hjust = 0) )
p <- p + ylab("Number of differentially expressed genes")
ggsave(file_tiff_countTable, width = 10, height = 6, device = "tiff", dpi = 900, compression = "lzw")

file_tiff_countLog2 <- paste(dir_output, "/classMode_counts.bar_log2.tiff", sep = "")
p <- ggplot(data = classMode_countMelt, aes(x = pattern, y = number, fill = pedigree) )
p <- p + geom_bar(stat = "identity")
p <- p + xlab("Expression pattern") + theme(axis.text.x = element_text(angle = -45, hjust = 0) )
#p <- p + ylab( expression (paste(log[2], "(Number of differentially expressed genes)"), sep = "") )
p <- p + ylab("Number of differentially expressed genes")
p <- p + scale_y_continuous(trans = log2_trans(), 
                            breaks = trans_breaks("log2", function(x) 2^x),
                            labels = trans_format("log2", math_format(2^.x) ) )
ggsave(file_tiff_countLog2, width = 10, height = 6, device = "tiff", dpi = 900, compression = "lzw")

file_tiff_countHeat <- paste(dir_output, "/classMode_counts.heat2.tiff", sep = "")
p <- ggplot(data = classMode_countMelt, aes(x = pedigree, y = pattern, fill = log2(number) ) )
p <- p + geom_raster()
#p <- p + theme_bw()
p <- p + theme_minimal()
p <- p + theme(axis.text.x = element_text(angle = -45, hjust = 0) )
p <- p + scale_fill_continuous(na.value = "grey", low = "white", high = "firebrick3", 
                               name = expression(paste(log[2], "(Number of DEGs)", sep = "")))
p <- p + xlab("Pedigree")
p <- p + ylab("Expression pattern")
ggsave(file_tiff_countHeat, width = 7, height = 7, device = "tiff", dpi = 900, compression = "lzw")

##########################################
# one example for test
##########################################
#unionID <- "Family07"
