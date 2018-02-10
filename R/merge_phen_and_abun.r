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

#' Merge phenotype data with relative abundances of blocks
#' 
#' Merge two data frames containing plant phenotypic data
#' and bacterial relative abundances
#' 
#' The function takes a list of variables that when combined must
#' uniqueley identify each row of both data frames. Those variable smust exist
#' in both data frames. 
#' 
#' It then constructs an index that it uses to match the columns. Abundance
#' data can be missing for some phenotypic observations, in which case
#' it will be replaced by zero.
#' 
#' @param Phen A data frame with phenotypic information
#' @param abun A data frame with abundance information
#' @param varnames Vector of charachter strings indicating which
#' variables to use to combine datasets. Rows where the exact same
#' combination of values for all variables occur are considered identical.
#' 
#' @author Sur Herrera Paredes
#' 
#' @keywords utils syncom colonization
#'
#' @export
merge_phen_and_abun <- function(Phen, abun, columns,
                                varnames = c("Bacteria", "Experiment", "StartP", "EndP")){
  ids <- apply(abun[,varnames],1,paste, collapse = "._.")
  if(length(ids) != nrow(abun))
    stop("ERROR: Variables do not uniqueley identify rows")
  row.names(abun) <- ids
  abun <- abun[,columns]

  abun <- abun[ apply(Phen[,varnames],1,paste, collapse = "._."), ]
  abun[is.na(abun)] <- 0
  row.names(abun) <- NULL

  Phen <- cbind(Phen, abun)

  return(Phen)
}
