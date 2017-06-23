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

#' Plot the expression of a set of genes
#'
#' Boxplot with dots of set of normalized RPKM values
#' based on a grouping factor.
#'
#' @param Dat Dataset. Should already be normalized
#' @param gene.set character vector of gene IDs to plot
#' @param groups character vector of column names of groups
#' to plot
#'
#' @author Sur Herrera Paredes
#'
#' @return A dataset
#'
#' @export
plot_set <- function(Dat,gene.set,groups){
  # Dat <- Dat.norm
  # gene.set <- pi.core
  # groups <- c("Bacteria","Phosphate")

  # Calculate mean expression per group
  Dat.sub <- AMOR::remove_taxons(Dat,setdiff(AMOR::taxa(Dat),gene.set))
  Dat.sub$Map$Group <- apply(Dat.sub$Map[,groups,drop = FALSE],
                             1,paste,collapse = ".")
  Dat.sub <- AMOR::pool_samples(Dat.sub,groups = "Group",FUN = mean)
  Map <- data.frame(do.call(rbind,strsplit(colnames(Dat.sub$Tab),"[.]")))
  colnames(Map) <- groups
  Map[, groups[1]] <- factor(Map[,groups[1]],
                             levels = levels(Dat$Map[,groups[1]]))
  Map[, groups[2]] <- factor(Map[,groups[2]],
                             levels = levels(Dat$Map[,groups[2]]))
  row.names(Map) <- colnames(Dat.sub$Tab)
  Dat.sub <- AMOR::create_dataset(Dat.sub$Tab,Map)

  # Prepare data for plot
  dat <- cbind(Dat.sub$Map, t(Dat.sub$Tab))
  dat <- reshape2::melt(dat,id.vars = groups,
                        variable.name = "Gene",
                        value.name = "z.score")
  dat$Bacteria <- factor(dat$Bacteria,levels = levels(Dat$Map$Bacteria))
  # dat$Bacteria <- factor(dat$Bacteria,levels = c("No Bacteria", "G1G2","G2G3","G1G3","G1N1",
  #                                                "G2N1","G3N1","N1N2","N2N3",
  #                                                "N3B1","B1B2","B2B3","B3G1","G2B3","G3B3"))
  # head(dat)

  # Plot
  p1 <- ggplot2::ggplot(dat,ggplot2::aes_string(x = groups[2],
                                                y = "z.score",
                                                color = groups[1])) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(),
                        size = 0.3) +
    ggplot2::geom_boxplot(outlier.colour = NA, fill = NA) +
    AMOR::theme_blackbox

  return(p1)
}
