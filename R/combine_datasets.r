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

#' Combine datasets
#'
#' Combine a list of Datasates.
#'
#' All datasets must have the same variables in the same order
#' in their Map attribute. Redundant samples (i.e. samples with the
#' same ID in more than one dataset), have their counts added.
#' Taxa neeed not be the same in all datasets. The union of taxa
#' will be used as long as all the intersections are identical
#' in the Tax elements. Taxa missing from a specific dataset
#' will be added with count zero for the samples in that dataset.
#'
#' IMPORTANT: There might be some issues when map is not null and there are
#' redundant samples
#'
#' @param datasets either a single Dataset object or a list
#' of dataset objects. See \code{\link{create_dataset}} for more info.
#' 
#' @author Sur Herrera Paredes
#' 
#' @keywords utilities
#' 
#' @seealso \code{\link{create_dataset}}
#'
#' @export
combine_datasets <- function(datasets){

  if(class(datasets) == "Dataset"){
    datasets <- list(datasets)
  }
  if(class(datasets) != "list"){
    stop("ERROR: datasets must be a list of datasets",call. = TRUE)
  }
  if(any(sapply(datasets,class) != "Dataset")){
    stop("ERROR: all elements in dataset must be Dataset objects", call. = TRUE)
  }

  # check metadata
  vars <- lapply(datasets,AMOR::variables)
  comparisons <- wheelP:::all_identical(vars)
  if(comparisons){
    # if all datasets have the same variables
    Map <- do.call(rbind,lapply(datasets,function(x) x$Map))
  }else{
    stop("ERROR: All datasets should have the same variables in the same order in the Map attribute.\n",call. = TRUE)
  }

  # Check taxonomy
  taxa.sets <- lapply(datasets,function(x) x$Tax)
  taxa.count <- table(do.call(c,lapply(taxa.sets, row.names)))
  taxa.list <- names(taxa.count)
  taxa.redundant <- names(taxa.count)[ which(taxa.count > 1) ]

  taxonomies <- lapply(taxa.redundant,wheelP:::extract_taxonomy,
                       datasets = datasets)
  # Remove nulls
  taxonomies <- lapply(taxonomies, function(x){ x[!sapply(x,is.null)]})
  # taxonomies
  taxonomies <- sapply(taxonomies, wheelP:::all_identical)
  if(any(!taxonomies)){
    cat(taxa.redundant[ !taxonomies ], "\n")
    stop("ERROR: Taxonomies are different between datasets\n",call. = TRUE)
  }else{
    Tax <- do.call(rbind,lapply(datasets,function(x) x$Tax))
    Tax <- Tax[ taxa.list, ]
  }

  # Homogenize samples
  Tab <- lapply(lapply(datasets, function(x) x$Tab),
                wheelP:::homogenize_taxa,
                taxa.list = taxa.list)
  # add redundant samples
  Tab <- lapply(Tab,reshape2::melt,
                value.name = "count",
                varnames = c("Taxon","Sample"),
                as.is = TRUE)
  Tab <- do.call(rbind,Tab)
  Tab <- reshape2::acast(Taxon ~ Sample,data = Tab,
                         fun.aggregate = sum,value.var = "count")
  Tab <- Tab[ taxa.list, ]
  if (!is.null(Map)) {
    Map <- Map[ colnames(Tab), ]
  }

  Dat <- create_dataset(Tab = Tab, Map = Map, Tax = Tax)

  return(Dat)
}




