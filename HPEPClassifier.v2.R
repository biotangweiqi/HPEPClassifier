#!/usr/bin/env Rscript

#########################
### command arguments ###
#########################
opt <- commandArgs(TRUE)
# input
file_pedigree_conf <- opt[1]
file_pedigree_data <- opt[2]
# output
file_expr_pattern  <- opt[3]

# example
# input
#file_pedigree_conf <- "Cross01.conf.xls"
#file_pedigree_data <- "00.DEGs_union_data/Cross01.DEG_expr.xls"
# output
#file_expr_pattern  <- "01.pedigreeExprPattern_data/Cross01.II32AxMH86.hpep.xls"

# create output dir
if(!dir.exists(dirname(file_expr_pattern) ) ) dir.create(dirname(file_expr_pattern), recursive = T)

###########################
### global arguments    ###
###########################

#################
### functions ###
#################
get_avg <- function(expr, group_name, sample_list){
  # calculate average of samples in one group
  sample_name <- colnames(expr)
  avg <- data.frame()
  for(group in group_name){
    idx <- which(sample_list$group == group)
    samples <- sample_list$sample[idx]
    #jdx <- which(sample_name %in% sample_list$ids[idx])
    if(length(samples) == 1) value <- expr[,samples]
    if(length(samples) > 1 ) value <- rowMeans(expr[,samples])
    #
    if(length(avg) >1) avg <- cbind(avg, value)
    if(length(avg) == 0) avg <- value
  }
  colnames(avg) <- group_name
  rownames(avg) <- rownames(expr)
  avg <- as.matrix(avg)
  return(avg)
}

#####################
### main pipeline ###
#####################
# read conf
cat("read configure file\n")
conf <- read.table(file_pedigree_conf, header = F, stringsAsFactors = F)
head(conf)
# global option from configure
pedigreeID <- conf[conf[,1] == "PedigreeID", 2]
maternal <- conf[conf[,1] == "Maternal", 2]
paternal <- conf[conf[,1] == "Paternal", 2]
hybrid   <- conf[conf[,1] == "Hybrid", 2]
sig.PvsM <- conf[conf[,1] == "PvsM", 2]
sig.PvsH <- conf[conf[,1] == "PvsH", 2]
sig.MvsH <- conf[conf[,1] == "MvsH", 2]

# read data
cat("read expr and significant data\n")
expr_sig <- read.table(file_pedigree_data, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_sig) <- expr_sig[,1]
head(expr_sig)

# check
cat("check configure and data colname\n")
if(length(pedigreeID) != 1) stop("PedigreeID is wrong\n")
if(any(!maternal %in% colnames(expr_sig)) ) stop("Sample name of maternal is not matched\n")
if(any(!paternal %in% colnames(expr_sig)) ) stop("Sample name of paternal is not matched\n")
if(any(!hybrid %in% colnames(expr_sig)) )   stop("Sample name of hybrid is not matched\n")
if(any(!sig.PvsM %in% colnames(expr_sig)) ) stop("Name of PvsM comparison is not matched\n")
if(any(!sig.PvsH %in% colnames(expr_sig)) ) stop("Name of PvsH comparison is not matched\n")
if(any(!sig.MvsH %in% colnames(expr_sig)) ) stop("Name of MvsH comparison is not matched\n")
cat("check is OK\n")

# construct sample list
cat("construct sample list\n")
sample_list <- data.frame(
  sample = c(maternal, paternal, hybrid), 
  group  = c(
    rep("M", length(maternal) ), 
    rep("P", length(paternal) ), 
    rep("H", length(hybrid) )
  ),
  check.names = F,
  stringsAsFactors = F
)
print(sample_list)

###############################################
# pre-process data
###############################################
# calculate mean of group
cat("calculate mean of group\n")
group_name <- unique(sample_list$group)
mtx.avg <- get_avg(as.matrix(expr_sig[,sample_list$sample] ), group_name, sample_list)
head(mtx.avg)

#
cat("get a matrix of significance information\n")
mtx.sig <- as.matrix(expr_sig[,c(sig.PvsM, sig.PvsH, sig.MvsH)])
colnames(mtx.sig) <- c("PvsM", "PvsH", "MvsH")
head(mtx.sig)

# cbind in one
avg_sig_matrix <- cbind(mtx.avg, mtx.sig)
head(avg_sig_matrix)
rm(mtx.avg, mtx.sig)

###############################################
# expression pattern
###############################################
# class-mode
identify_express_pattern <- function(x){
  classCode <- NA
  # AHP
  if(x["P"]> x["M"] && x["H"]> x["P"] && x["H"]> x["M"] &&  x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.01"
  if(x["M"]> x["P"] && x["H"]> x["P"] && x["H"]> x["M"] &&  x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.02"
  if(x["P"]>=x["M"] && x["H"]> x["P"] && x["H"]> x["M"] && !x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.03"
  if(x["M"]>=x["P"] && x["H"]> x["P"] && x["H"]> x["M"] && !x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.04"
  # HP
  if(x["P"]>=x["M"] && x["H"]>=x["P"] && x["H"]> x["M"] && !x["PvsM"] && !x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.05"
  if(x["M"]>=x["P"] && x["H"]> x["P"] && x["H"]>=x["M"] && !x["PvsM"] &&  x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.06"
  if(x["P"]> x["M"] && x["H"]>=x["P"] && x["H"]> x["M"] &&  x["PvsM"] && !x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.07"
  if(x["P"]> x["M"] && x["P"]>=x["H"] && x["H"]> x["M"] &&  x["PvsM"] && !x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.08"
  if(x["M"]> x["P"] && x["H"]> x["P"] && x["H"]>=x["M"] &&  x["PvsM"] &&  x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.09"
  if(x["M"]> x["P"] && x["H"]> x["P"] && x["M"]>=x["H"] &&  x["PvsM"] &&  x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.10"
  # MP
  if(x["P"]> x["M"] && x["P"]>=x["H"] && x["H"]>=x["M"] &&  x["PvsM"] && !x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.11"
  if(x["P"]> x["M"] && x["P"]> x["H"] && x["H"]> x["M"] &&  x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.12"
  if(x["M"]> x["P"] && x["H"]> x["P"] && x["M"]> x["H"] &&  x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.13"
  if(x["M"]> x["P"] && x["H"]>=x["P"] && x["M"]>=x["H"] &&  x["PvsM"] && !x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.14"
  # LP
  if(x["M"]> x["P"] && x["H"]>=x["P"] && x["M"]> x["H"] &&  x["PvsM"] && !x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.15"
  if(x["M"]> x["P"] && x["P"]>=x["H"] && x["M"]> x["H"] &&  x["PvsM"] && !x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.16"
  if(x["P"]> x["M"] && x["P"]> x["H"] && x["H"]>=x["M"] &&  x["PvsM"] &&  x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.17"
  if(x["P"]> x["M"] && x["P"]> x["H"] && x["M"]>=x["H"] &&  x["PvsM"] &&  x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.18"
  if(x["M"]>=x["P"] && x["P"]>=x["H"] && x["M"]> x["H"] && !x["PvsM"] && !x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.19"
  if(x["P"]>=x["M"] && x["P"]> x["H"] && x["M"]>=x["H"] && !x["PvsM"] &&  x["PvsH"] && !x["MvsH"]) classCode <- "ClassNo.20"
  # BLP
  if(x["M"]>=x["P"] && x["P"]> x["H"] && x["M"]> x["H"] && !x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.21"
  if(x["P"]>=x["M"] && x["P"]> x["H"] && x["M"]> x["H"] && !x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.22"
  if(x["M"]> x["P"] && x["P"]> x["H"] && x["M"]> x["H"] &&  x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.23"
  if(x["P"]> x["M"] && x["P"]> x["H"] && x["M"]> x["H"] &&  x["PvsM"] &&  x["PvsH"] &&  x["MvsH"]) classCode <- "ClassNo.24"
  #
  return(classCode)
}

#
cat("identify expression pattern and add ClassCode\n")
classCodes <- apply(avg_sig_matrix, 1, identify_express_pattern)
expr_sig$ClassCode <- classCodes
head(expr_sig)

#
classCode_list <- data.frame(
  ExprPattern = c(
    # AHP
    "H>P>M",  # No.01
    "H>M>P",  # No.02
    "H>P>=M", # No.03
    "H>M>=P", # No.04
    # HP
    "H>=P>=M",# No.05
    "H>=M>=P",# No.06
    "H>=P>M", # No.07
    "P>=H>M", # No.08
    "H>=M>P", # No.09
    "M>=H>P", # No.10
    # MP
    "P>=H>=M",# No.11
    "P>H>M",  # No.12
    "M>H>P",  # No.13
    "M>=H>=P",# No.14
    # LP
    "M>H>=P", # No.15
    "M>P>=H", # No.16
    "P>H>=M", # No.17
    "P>M>=H", # No.18
    "M>=P>=H",# No.19
    "P>=M>=H",# No.20
    # BLP
    "M>=P>H", # No.21
    "P>=M>H", # No.22
    "M>P>H",  # No.23
    "P>M>H"   # No.24
  ),
  Superclass = c(
    "AHP", "AHP", "AHP", "AHP",
    "HP", "HP", "HP", "HP", "HP", "HP",
    "MP", "MP", "MP", "MP",
    "LP", "LP", "LP", "LP", "LP", "LP",
    "BLP", "BLP", "BLP", "BLP"
  ),
  stringsAsFactors = F,
  check.names = F
)
rownames(classCode_list) <- c(
  "ClassNo.01",
  "ClassNo.02",
  "ClassNo.03",
  "ClassNo.04",
  "ClassNo.05",
  "ClassNo.06",
  "ClassNo.07",
  "ClassNo.08",
  "ClassNo.09",
  "ClassNo.10",
  "ClassNo.11",
  "ClassNo.12",
  "ClassNo.13",
  "ClassNo.14",
  "ClassNo.15",
  "ClassNo.16",
  "ClassNo.17",
  "ClassNo.18",
  "ClassNo.19",
  "ClassNo.20",
  "ClassNo.21",
  "ClassNo.22",
  "ClassNo.23",
  "ClassNo.24"
)
classCode_list

# add exprPattern and superclass information
cat("add exprPattern and superclass information\n")
exprPattern <- apply(as.matrix(classCodes), 1, function(x) ifelse(is.na(x), NA, classCode_list[x,1]) )
superClass  <- apply(as.matrix(classCodes), 1, function(x) ifelse(is.na(x), NA, classCode_list[x,2]) )
expr_sig$ExprPattern <- exprPattern
expr_sig$SuperClass  <- superClass
head(expr_sig)

# order by classCode
cat("order by classCode\n")
expr_sig <- expr_sig[order(expr_sig$ClassCode),]
head(expr_sig)

# write out
cat("write out\n")
write.table(expr_sig, file_expr_pattern, sep = "\t", quote = F, row.names = F)
