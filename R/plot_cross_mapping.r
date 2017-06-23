# (C) Copyright 2017 Sur Herrera Paredes
# 
# This file is part of wheelP.
# 
# wheelP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# wheelP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with wheelP.  If not, see <http://www.gnu.org/licenses/>.

#' Plot mapping id
#'
#' @export
plot_cross_mapping <- function(Tax, dist.seq, id = 98){

  bar.dat <- data.frame(Strain = Tax$ID[ order(Tax$Block) ], y = 1)
  bar.dat$Strain <- factor(bar.dat$Strain,
                           levels = as.character(Tax$ID[ order(Tax$Block) ]))

  dist.seq <- dist.seq[ levels(bar.dat$Strain), levels(bar.dat$Strain) ]

  dist.seq <- as.dist(dist.seq)
  dist.seq <- dist2df(dist.seq)
  #head(dist.seq)

  dist.plot <- subset(dist.seq, value >= id)
  dist.plot$row <- factor(dist.plot$row, levels = levels(bar.dat$Strain))
  dist.plot$col <- factor(dist.plot$col, levels = levels(bar.dat$Strain))
  dist.plot <- dist.plot[ !apply(is.na(dist.plot),1,any), ]
  #dist.plot

  p1 <- ggplot(bar.dat,
               aes(x = Strain, y = 0.02)) +
    geom_bar(stat = "identity", fill = "black") +
    geom_text(aes(y = 0.01, label = Strain), col = "white",
              angle = 90) +
    geom_curve(data = dist.plot,
               aes(x = row, xend = col, y = 0, yend = 0),
               angle = 90, curvature = -0.5, size = 2, alpha = 0.2) +
    ylim(-0.1,0.02) +
    theme_blackbox +
    theme(panel.background = element_blank(),
          axis.text.x = element_text(angle = 90))

  return(p1)
}
