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

#' Plot binary association results
#'
#' @param test result from test all strains
#' @param qval_trhes q-value threshold
#'
#' @author Sur Herrera Paredes
#'
#' @return A ggplot2 object
#'
#' @export
plot_res <- function(test,qval_thres = 0.05){
  Dat <- reshape2::melt(test$RES,value.name = "Change",variable.name = "Condition",id.vars = "Strain")
  Dat$Fold <- round(exp(Dat$Change),2)
  Dat$qval <- reshape2::melt(test$QVAL,value.name = "qval",variable.name = "Condition",id.vars = "Strain")$qval
  Dat$Effect <- "none"
  Dat$Effect <- factor(Dat$Effect,levels = c("Decrease","Increase","none"))
  Dat$Effect[ Dat$qval < qval_thres & Dat$Fold > 1 ] <- "Increase"
  Dat$Effect[ Dat$qval < qval_thres & Dat$Fold < 1 ] <- "Decrease"

  p1 <- ggplot2::ggplot(Dat,aes(x = Condition, y = Strain)) +
    ggplot2::geom_tile(aes(fill = Effect)) +
    ggplot2::geom_text(aes(label = Fold)) +
    ggplot2::scale_fill_manual(values = c("steelblue","red","white")) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = element_text(face="bold",size = 8, angle = 90),
                   axis.text.y = element_text(face="bold",size = 8))
  return(p1)
}
