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

#' Plot heatmap of all effects
#'
#' @export
plot_combined_effects <- function(dat, scale = TRUE){

  # Add codition name var
  dat <- lapply(dat, wheelP:::add_condition)

  # Combine
  tab <- lapply(dat,acast, formula = SynCom ~ Condition, value.var = "Estimate")
  tab <- do.call(cbind,tab)
  colnames(tab) <- paste( rep(names(dat),each = 4),
                          colnames(tab),sep = ".")

  # Plot
  pal <- colorRampPalette(colors = c("#8e0152","#de77ae",
                                     "#f7f7f7","#7fbc41","#276419"))
  #tiff("heatmap_syncom_effects.tif",width = 2000, height = 1500, res = 250)
  heatmap.2(x = scale(tab, center = FALSE, scale = scale),
            main = "Scaled SynCom effects\non plant phenotypes",
            dendrogram = "row", trace = "none", col = pal(50),
            margins = c(12,5), Colv = FALSE)
  #dev.off()

  return(tab)
}


#' add condition
#'
#' Internal
add_condition <- function(Dat){

  Dat$Condition <- NULL
  Dat$Condition[ Dat$StartP == "-Pi,0.5%Suc" & Dat$EndP == "100 uM,0%Suc" ] <- "minusP_100uM"
  Dat$Condition[ Dat$StartP == "+Pi,0.5%Suc" & Dat$EndP == "100 uM,0%Suc" ] <- "plusP_100uM"
  Dat$Condition[ Dat$StartP == "-Pi,0.5%Suc" & Dat$EndP == "30 uM,0%Suc" ] <- "minusP_30uM"
  Dat$Condition[ Dat$StartP == "+Pi,0.5%Suc" & Dat$EndP == "30 uM,0%Suc" ] <- "plusP_30uM"

  return(Dat)
}
