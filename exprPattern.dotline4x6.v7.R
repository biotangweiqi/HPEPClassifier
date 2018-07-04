#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# data
# AHP
data.No.01 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(2, 3, 1, 1.5),
  g = c("a", "a", "a", "c")
)
data.No.02 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(1, 3, 2, 1.5),
  g = c("a", "a", "a", "c")
)
data.No.03 <- data.frame(
  x = factor(c("P", "H", "M", "P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(2, 3, 1.5,1.75, 3, 1.75, 1.75),
  g = c("a", "a", "a", "b", "b", "b", "c")
)
data.No.04 <- data.frame(
  x = factor(c("P", "H", "M", "P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(1.5, 3, 2, 1.75, 3, 1.75, 1.75),
  g = c("a", "a", "a", "b", "b", "b", "c")
)
# HP
data.No.05 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(2, 2.5, 1.5, 1.75),
  g = c("a", "a", "a", "c")
)
data.No.06 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(1.5, 2.5, 2, 1.75),
  g = c("a", "a", "a", "c")
)
data.No.07 <- data.frame(
  x = factor(c("P", "H", "M", "P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(2, 2.5, 1, 2, 2, 1, 1.5),
  g = c("a", "a", "a", "b", "b", "b", "c")
)
data.No.08 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(2.2, 2, 1, 1.6),
  g = c("a", "a", "a", "c")
)
data.No.09 <- data.frame(
  x = factor(c("P", "H", "M", "P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(1, 2.5, 2, 1, 2, 2, 1.5),
  g = c("a", "a", "a", "b", "b", "b", "c")
)
data.No.10 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(1, 2, 2.2, 1.6),
  g = c("a", "a", "a", "c")
)
# MP
data.No.11 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(2.5, 2, 1.5, 2),
  g = c("a", "a", "a", "c")
)
data.No.12 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(3, 2, 1, 2),
  g = c("a", "a", "a", "c")
)
data.No.13 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(1, 2, 3, 2),
  g = c("a", "a", "a", "c")
)
data.No.14 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(1.5, 2, 2.5, 2),
  g = c("a", "a", "a", "c")
)

# LP
data.No.15 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(1.8, 2, 3, 2.4),
  g = c("a", "a", "a", "c")
)
data.No.16 <- data.frame(
  x = factor(c("P", "H", "M", "P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(2, 1.5, 3, 2, 2, 3, 2.5),
  g = c("a", "a", "a", "b", "b", "b", "c")
)
data.No.17 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(3, 2, 1.8, 2.4),
  g = c("a", "a", "a", "c")
)
data.No.18 <- data.frame(
  x = factor(c("P", "H", "M", "P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(3, 1.5, 2, 3, 2, 2, 2.5),
  g = c("a", "a", "a", "b", "b", "b", "c")
)
data.No.19 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(2, 1.5, 2.5, 2.25),
  g = c("a", "a", "a", "c")
)
data.No.20 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(2.5, 1.5, 2, 2.25),
  g = c("a", "a", "a", "c")
)

# BLP
data.No.21 <- data.frame(
  x = factor(c("P", "H", "M", "P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(2, 1, 2.5, 2.25, 1, 2.25, 2.25),
  g = c("a", "a", "a", "b", "b", "b", "c")
)
data.No.22 <- data.frame(
  x = factor(c("P", "H", "M", "P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(2.5, 1, 2, 2.25, 1, 2.25, 2.25),
  g = c("a", "a", "a", "b", "b", "b", "c")
)
data.No.23 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(2, 1, 3, 2.5),
  g = c("a", "a", "a", "c")
)
data.No.24 <- data.frame(
  x = factor(c("P", "H", "M", "H"), levels = c("P", "H", "M")),
  y = c(3, 1, 2, 2.5),
  g = c("a", "a", "a", "c")
)

# draw dotline by ggplot2
draw_dotline <- function(data, title){
  p <- ggplot(data, aes(x , y, group = g) ) 
  p <- p + geom_line(aes(linetype = g), color = "orangered", size = 2)
  #p <- p + geom_point(size = 8, color = "black", shape = 20)
  p <- p + geom_point(aes(color = g, size = g, shape = g, alpha = g) )
  p <- p + ylim(0.9,3.1)
  #p <- p + scale_color_discrete(guide = F)  
  #p <- p + scale_size_discrete(guide = F)
  #p <- p + scale_alpha_discrete(guide = F)
  #p <- p + scale_shape_discrete(guide = F)
  p <- p + scale_linetype_discrete(guide = F)
  p <- p + scale_color_manual(guide = F, values = c("black", "black", "navy"), limits = c("a", "b", "c") )
  p <- p + scale_alpha_manual(guide = F, values = c(1, 1, 0.7), limits = c("a", "b", "c") )
  p <- p + scale_shape_manual(guide = F, values = c(20, 20, 18), limits = c("a", "b", "c") )
  p <- p + scale_size_manual(guide = F, values = c(8, 8, 12), limits = c("a", "b", "c") )
  p <- p + theme_minimal()
  p <- p + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, size = rel(2)) )
  p <- p + theme(axis.text.x = element_text(size = rel(2) ) )
  p <- p + theme(axis.text.y = element_blank() )
  p <- p + theme(axis.title = element_blank() )
  p <- p + theme(panel.grid.major.x = element_blank() )
  p
}

#
p1 <- draw_dotline(data.No.01, "ClassNo.01 H>P>M (AHP)")
p2 <- draw_dotline(data.No.02, "ClassNo.02 H>M>P (AHP)")
p3 <- draw_dotline(data.No.03, "ClassNo.03 H>P>=M (AHP)")
p4 <- draw_dotline(data.No.04, "ClassNo.04 H>M>=P (AHP)")
#
p5 <- draw_dotline(data.No.05, "ClassNo.05 H>=P>=M (HP)")
p6 <- draw_dotline(data.No.06, "ClassNo.06 H>=M>=P (HP)")
p7 <- draw_dotline(data.No.07, "ClassNo.07 H>=P>M (HP)")
p8 <- draw_dotline(data.No.08, "ClassNo.08 P>=H>M (HP)")
p9  <- draw_dotline(data.No.09, "ClassNo.09 H>=M>P (HP)")
p10 <- draw_dotline(data.No.10, "ClassNo.10 M>=H>P (HP)")
#
p11 <- draw_dotline(data.No.11, "ClassNo.11 P>=H>=M (MP)")
p12 <- draw_dotline(data.No.12, "ClassNo.12 P>H>M (MP)")
p13 <- draw_dotline(data.No.13, "ClassNo.13 M>H>P (MP)")
p14 <- draw_dotline(data.No.14, "ClassNo.14 M>=H>=P (MP)")
#
p15 <- draw_dotline(data.No.15, "ClassNo.15 M>H>=P (LP)")
p16 <- draw_dotline(data.No.16, "ClassNo.16 M>P>=H (LP)")
p17 <- draw_dotline(data.No.17, "ClassNo.17 P>H>=M (LP)")
p18 <- draw_dotline(data.No.18, "ClassNo.18 P>M>=H (LP)")
p19 <- draw_dotline(data.No.19, "ClassNo.19 M>=P>=H (LP)")
p20 <- draw_dotline(data.No.20, "ClassNo.20 P>=M>=H (LP)")
#
p21 <- draw_dotline(data.No.21, "ClassNo.21 M>=P>H (BLP)")
p22 <- draw_dotline(data.No.22, "ClassNo.22 P>=M>H (BLP)")
p23 <- draw_dotline(data.No.23, "ClassNo.23 M>P>H (BLP)")
p24 <- draw_dotline(data.No.24, "ClassNo.24 P>M>H (BLP)")

# determine total width and heigth
each_width   <- 4
each_height  <- 3
col_grids    <- 4
row_grids    <- 6
total_width  <- col_grids*each_width
total_height <- row_grids*each_height

# open device
tiff(filename = "dotline4x6_output.tiff",
     width = total_width,
     height = total_height,
     res = 1200,
     units = "in",
     compression = "lzw"
)

# width and height of each grid
grid_width  <- rep(each_width, col_grids)
grid_height <- rep(each_height, row_grids) 

# use gridExtra
#g <- gTree(children = gList(diag))
gridExtra::grid.arrange(
  p1,  p2,  p3,  p4,
  p5,  p6,  p7,  p8,
  p9,  p10, p11, p12,
  p13, p14, p15, p16,
  p17, p18, p19, p20,
  p21, p22, p23, p24,
  ncol = col_grids, nrow = row_grids, 
  widths = grid_width, heights = grid_height 
)

# close device
dev.off()

