# pct_arrows.R
# Direction arrows for % change panels: explain what lies above / below (or left / right of) the
# 0 line. Drawn just outside the panel edge, anchored at the 0 line, so they never cover data.
# `zero` is the data value of the 0% line: 0 on a linear % scale, 1 on a log-scaled ratio axis
# (annotation positions are given in data units and transformed by the scale).
# Panels using these need coord_cartesian(clip = "off").
#   pct_arrows_y(): vertical, outside the right edge   (+ pct_arrow_margin_y() for room)
#   pct_arrows_x(): horizontal, below the x-axis title  (+ pct_arrow_margin_x() for room)
library(grid)

.arrow_gp <- function(col) gpar(col = col, fill = col, lwd = 1.9)

pct_arrows_y <- function(up, down, zero = 0, col = "#4D4D4D", size = 11.5) {
  ar <- arrow(length = unit(2.8, "mm"), type = "closed")
  xa <- unit(1, "npc") + unit(3.5, "mm"); xt <- unit(1, "npc") + unit(6.5, "mm")
  up_g <- grobTree(
    segmentsGrob(x0 = xa, x1 = xa, y0 = unit(1.2, "mm"), y1 = unit(16, "mm"), arrow = ar, gp = .arrow_gp(col)),
    textGrob(up, x = xt, y = unit(16, "mm"), hjust = 0, vjust = 1, gp = gpar(fontsize = size, col = col, lineheight = 0.9)))
  dn_g <- grobTree(
    segmentsGrob(x0 = xa, x1 = xa, y0 = unit(1, "npc") - unit(1.2, "mm"), y1 = unit(1, "npc") - unit(16, "mm"),
                 arrow = ar, gp = .arrow_gp(col)),
    textGrob(down, x = xt, y = unit(1, "npc") - unit(16, "mm"), hjust = 0, vjust = 0,
             gp = gpar(fontsize = size, col = col, lineheight = 0.9)))
  list(ggplot2::annotation_custom(up_g, xmin = -Inf, xmax = Inf, ymin = zero, ymax = Inf),
       ggplot2::annotation_custom(dn_g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = zero))
}
pct_arrow_margin_y <- function() ggplot2::theme(plot.margin = ggplot2::margin(8, 125, 8, 8))

pct_arrows_x <- function(right, left, zero = 0, col = "#4D4D4D", size = 11.5) {
  ar <- arrow(length = unit(2.8, "mm"), type = "closed")
  ya <- unit(0, "npc") - unit(18, "mm")
  r_g <- grobTree(
    segmentsGrob(x0 = unit(1.2, "mm"), x1 = unit(20, "mm"), y0 = ya, y1 = ya, arrow = ar, gp = .arrow_gp(col)),
    textGrob(right, x = unit(22, "mm"), y = ya, hjust = 0, vjust = 0.5, gp = gpar(fontsize = size, col = col)))
  l_g <- grobTree(
    segmentsGrob(x0 = unit(1, "npc") - unit(1.2, "mm"), x1 = unit(1, "npc") - unit(20, "mm"), y0 = ya, y1 = ya,
                 arrow = ar, gp = .arrow_gp(col)),
    textGrob(left, x = unit(1, "npc") - unit(22, "mm"), y = ya, hjust = 1, vjust = 0.5, gp = gpar(fontsize = size, col = col)))
  list(ggplot2::annotation_custom(r_g, xmin = zero, xmax = Inf, ymin = -Inf, ymax = Inf),
       ggplot2::annotation_custom(l_g, xmin = -Inf, xmax = zero, ymin = -Inf, ymax = Inf))
}
pct_arrow_margin_x <- function() ggplot2::theme(plot.margin = ggplot2::margin(8, 12, 40, 8))

# Box the "Red" tick label of a discrete light-color x-axis in red with white text. The tick
# label itself is blanked and redrawn just below the panel (needs coord_cartesian(clip = "off")).
# `bottom` = data value of the panel's lower edge: -Inf on a linear axis, 0 on a log axis.
red_tick_box <- function(levels, red = "Red", fill = "#B2182B", size = 10, bottom = -Inf) {
  pos <- match(red, levels)
  lab <- grobTree(
    roundrectGrob(x = 0.5, y = unit(0, "npc") - unit(9, "pt"), width = unit(30, "pt"), height = unit(13, "pt"),
                  r = unit(2, "pt"), gp = gpar(col = NA, fill = fill)),
    textGrob(red, x = 0.5, y = unit(0, "npc") - unit(9, "pt"), gp = gpar(fontsize = size, col = "white")))
  list(ggplot2::scale_x_discrete(labels = function(l) ifelse(l == red, "", l)),
       ggplot2::annotation_custom(lab, xmin = pos - 0.5, xmax = pos + 0.5, ymin = bottom, ymax = Inf))
}
