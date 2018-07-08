#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# data
# AHP
data.AHP <- data.frame(
  x = factor(c("LP", "H", "HP", "H"), levels = c("LP", "H", "HP")),
  y = c(1, 3, 2, 1.5),
  g = c("a", "a", "a", "c")
)
# HP
data.HP <- data.frame(
  x = factor(c("LP", "H", "HP", "H"), levels = c("LP", "H", "HP")),
  y = c(1, 2, 2, 1.5),
  g = c("a", "a", "a", "c")
)
# MP
data.MP <- data.frame(
  x = factor(c("LP", "H", "HP", "H"), levels = c("LP", "H", "HP")),
  y = c(1, 2, 3, 2),
  g = c("a", "a", "a", "c")
)
# LP
data.LP <- data.frame(
  x = factor(c("LP", "H", "HP", "H"), levels = c("LP", "H", "HP")),
  y = c(1, 1, 2, 1.5),
  g = c("a", "a", "a", "c")
)
# BLP
data.BLP <- data.frame(
  x = factor(c("LP", "H", "HP", "H"), levels = c("LP", "H", "HP")),
  y = c(2, 1, 3, 2.5),
  g = c("a", "a", "a", "c")
)

# draw dotline by ggplot2
draw_dotline <- function(data, title){
  p <- ggplot(data, aes(x , y, group = g) ) 
  p <- p + geom_line(aes(linetype = g), color = "orangered", size = 2)
  #p <- p + geom_point(size = 8, color = "black", shape = 20)
  p <- p + geom_point(aes(color = g, size = g, shape = g, alpha = g) )
  p <- p + ylim(0.9,3.1)
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
p1 <- draw_dotline(data.AHP, "Superclass AHP")
p2 <- draw_dotline(data.HP,  "Superclass HP")
p3 <- draw_dotline(data.MP,  "Superclass MP")
p4 <- draw_dotline(data.LP,  "Superclass LP")
p5 <- draw_dotline(data.BLP, "Superclass BLP")

# determine total width and heigth
each_width   <- 4
each_height  <- 3
col_grids    <- 5
row_grids    <- 1
total_width  <- col_grids*each_width
total_height <- row_grids*each_height

# open device
tiff(filename = "dotline1x5_output.tiff",
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
  p1,  p2,  p3,  p4, p5,
  ncol = col_grids, nrow = row_grids, 
  widths = grid_width, heights = grid_height 
)

# close device
dev.off()

