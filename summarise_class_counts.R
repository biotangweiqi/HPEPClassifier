#!/usr/bin/env Rscript
library(ggplot2)
library(gridExtra)

#########################
### command arguments ###
#########################
opt <- commandArgs(TRUE)
# input
file_expr_pattern <- opt[1]
# output
output_prefix     <- opt[2]

# example
# input
#file_expr_pattern <- "01.pedigreeExprPattern_data/Cross01.II32AxMH86.hpep.xls"
# output
#output_prefix     <- "02.pedigreeExprPattern_summary/Cross01.II32AxMH86"

# output files
file_class_count  <- paste(output_prefix, ".hpep_class_count.xls", sep = "")
file_super_count  <- paste(output_prefix, ".hpep_superClass_count.xls", sep = "")
file_counts_tiff  <- paste(output_prefix, ".hpep_count_bar.tiff", sep = "")
file_rate_tiff    <- paste(output_prefix, ".hpep_rate_bar.tiff", sep = "")

# create output dir
if(!dir.exists(dirname(output_prefix) ) ) dir.create(dirname(output_prefix), recursive = T)

###########################
### global arguments    ###
###########################
total_width  <- 15
total_height <- 8

###########################
### list of classCode   ###
###########################
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
classCode_list$ExprPattern <- factor(classCode_list$ExprPattern, levels = as.character(classCode_list$ExprPattern) )
classCode_list$Superclass <- factor(classCode_list$Superclass, levels = c("AHP", "HP", "MP", "LP", "BLP"))

#
superClass_list <- data.frame(
  SuperClass = c("AHP", "HP", "MP", "LP", "BLP")
)
superClass_list
superClass_list$SuperClass <- factor(superClass_list$SuperClass, levels = superClass_list$SuperClass)

#####################
### main pipeline ###
#####################
# read data
cat("read data\n")
data_expr_pattern <- read.table(file_expr_pattern, header = T, sep = "\t", stringsAsFactors = F, check.names = F)
head(data_expr_pattern)

#
totalNum <- nrow(data_expr_pattern)

#
class_counts <- c()
for(classCode in rownames(classCode_list) ){
  count <- length( which(data_expr_pattern$ClassCode == classCode) )
  if(length(class_counts) > 0)  class_counts <- c(class_counts, count)
  if(length(class_counts) == 0) class_counts <- c(count)
}
classCode_list$count <- class_counts
classCode_list$rate  <- class_counts / totalNum
#classCode_list

#
superClass_counts <- c()
for(superClass in superClass_list$SuperClass){
  count <- length( which(data_expr_pattern$SuperClass == superClass) )
  if(length(superClass_counts) > 0)  superClass_counts <- c(superClass_counts, count)
  if(length(superClass_counts) == 0) superClass_counts <- c(count)
}
superClass_list$count <- superClass_counts
superClass_list$rate  <- superClass_counts / totalNum
#superClass_list

# draw
p <- ggplot(data = classCode_list, aes(x = ExprPattern, y = count) )
p <- p + geom_bar(stat = "identity")
p <- p + facet_wrap(~ Superclass, ncol = 5, scales = "free_x")
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1) )
p <- p + theme(axis.title = element_text(size = rel(1.8) ) )
p <- p + theme(axis.text = element_text(size = rel(1.36) ) )
p <- p + theme(strip.text = element_text(size = rel(1.5) ) )
p <- p + xlab("Expression Pattern") + ylab("Count")

#
p2 <- ggplot(data = superClass_list, aes(x = SuperClass, y = count) )
p2 <- p2 + geom_bar(stat = "identity")
p2 <- p2 + theme(axis.title.y = element_blank() )
p2 <- p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- p2 + theme(axis.title = element_text(size = rel(1.8) ) )
p2 <- p2 + theme(axis.text = element_text(size = rel(1.36) ) )

# open device
tiff(filename = file_counts_tiff,
     width = total_width,
     height = total_height,
     res = 1200,
     units = "in",
     compression = "lzw"
)

# width and height of each grid
grid_width  <- c(9, 3)

# use gridExtra
#g <- gTree(children = gList(diag))
gridExtra::grid.arrange(
  p,  p2,
  ncol = 2, 
  widths = grid_width 
)

# close device
dev.off()

# draw
p <- ggplot(data = classCode_list, aes(x = ExprPattern, y = rate) )
p <- p + geom_bar(stat = "identity")
p <- p + facet_wrap(~ Superclass, ncol = 5, scales = "free_x")
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1) )
p <- p + theme(axis.title = element_text(size = rel(1.8) ) )
p <- p + theme(axis.text = element_text(size = rel(1.36) ) )
p <- p + theme(strip.text = element_text(size = rel(1.5) ) )
p <- p + xlab("Expression Pattern") + ylab("Rate")

#
p2 <- ggplot(data = superClass_list, aes(x = SuperClass, y = rate) )
p2 <- p2 + geom_bar(stat = "identity")
p2 <- p2 + theme(axis.title.y = element_blank() )
p2 <- p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- p2 + theme(axis.title = element_text(size = rel(1.8) ) )
p2 <- p2 + theme(axis.text = element_text(size = rel(1.36) ) )

# open device
tiff(filename = file_rate_tiff,
     width = total_width,
     height = total_height,
     res = 1200,
     units = "in",
     compression = "lzw"
)

# width and height of each grid
grid_width  <- c(9, 3)

# use gridExtra
#g <- gTree(children = gList(diag))
gridExtra::grid.arrange(
  p,  p2,
  ncol = 2, 
  widths = grid_width 
)

# close device
dev.off()

#
data_summary <- data.frame(
  name  = paste(rownames(classCode_list), classCode_list$ExprPattern, sep = "."),
  count = class_counts,
  check.names = F
)

# add number of missing data
numOfNA <- length(which(is.na(data_expr_pattern$ClassCode) ) )
data_summary <- rbind(data_summary, data.frame(t(c(name = "Ambiguous", count = numOfNA))))

# add total number
data_summary <- rbind(data_summary, data.frame(t(c(name = "Total", count = totalNum))))

# add rate
data_summary$rate <- sprintf("%.4f%%", 100*as.numeric(data_summary$count)/totalNum)
head(data_summary)

# write out
write.table(data_summary, file_class_count, quote = F, sep = "\t", row.names = F)

#
superClass_list$rate <- sprintf("%.4f%%", 100*as.numeric(superClass_list$count)/totalNum)
head(superClass_list)

# write out
write.table(superClass_list, file_super_count, quote = F, sep = "\t", row.names = F, append = TRUE)
