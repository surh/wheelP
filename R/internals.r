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

#' Extract taxonomy info for one taxon
#'
#' Internal
extract_taxonomy <- function(taxon,datasets){
  taxonomies <- lapply(datasets,
                       function(x,taxon = taxon){
                         if(taxon %in% AMOR::taxa(x)){
                           return(as.character(as.matrix(x$Tax[ taxon, ])))
                         }else{
                           return(NULL)
                         }
                       } ,
                       taxon = taxon)

  return(taxonomies)
}

#' Make sure all elements in list are identical
#'
#' Internal
all_identical <- function(x){
  #x <- extract_taxonomy("contaminant1", datasets)

  comparisons <- combn(x = 1:length(x),m = 2)
  comparisons
  comparisons <- apply(comparisons,2,
                       function(y,dat = x) identical(dat[[y[1]]],dat[[y[2]]]))
  comparisons <- all(comparisons)
  comparisons

  return(comparisons)
}

#' homogenize taxa in abundance matrix
#'
#' @param Tab Abundance matrix
#' @param taxa.list character vector of taxa
#'
#' Internal
homogenize_taxa <- function(Tab, taxa.list){
  to_add <- setdiff(taxa.list, row.names(Tab))
  mat <- matrix(0, nrow = length(to_add), ncol = dim(Tab)[2])
  colnames(mat) <- colnames(Tab)
  row.names(mat) <- to_add

  Tab <- rbind(Tab,mat)
  Tab <- Tab[ taxa.list, , drop = FALSE]

  return(Tab)
}
