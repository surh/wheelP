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

#' Plot results from binary association
#'
#' Generates density plot of Pi content in four conditions with
#' and without bacteria
#'
#' @param Dat A data.frame. Requires columns StartP and EndP defining
#' conditions and used for faceting. Column Pi_content indicates the
#' phosphate concentration. Column bacteria indicates whether bacteria
#' was added or not.
#' @param xbreaks vector of values to be indicated in x-axis.
#'
#' @return A ggplot2 plot object.
#'
#' @author Sur Herrera Paredes
#' 
#' @keywords binP plots
#'
#' @export
plot_binary <- function (Dat,xbreaks = c(1, 2, 5, 10, 20, 50)){
  p1 <- ggplot(Dat, aes(x = Pi_content)) +
    facet_grid(StartP ~ EndP) +
    scale_x_log10(breaks = xbreaks) +
    geom_density(aes(fill = Bacteria),alpha = 0.5) +
    theme(panel.background = element_rect(fill = NA),
          panel.border = element_rect(color = "black",fill = NA, size = 3),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(face = "bold", color = "black", size = 24),
          axis.text.x = element_text(color = "black", size = 20),
          strip.text = element_text(color = "black", size = 32, face = "bold"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20))
  return(p1)
}
