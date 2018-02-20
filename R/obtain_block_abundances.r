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

#' Obtain block abundances
#'
#' Takes a Dataset object, pools samples and taxa according to
#' specified variables, and returns  adata frame that contains
#' the average abundance by taxa group per group of samples
#' 
#' @param Dat a dataset object
#' @param varnames Vector of character strings that name the variables
#' to be used to pool samples. Must correspond to variables in
#' the Map attribute of Dat.
#' @param taxa.group Character string indicating the variable for
#' collapsing taxa.
#' @param sep Separator character for creating composite variables.
#' @param taxa2rm Taxa IDs to remove.
#' 
#' @return A data frame that contains the average abundance by taxa group
#' per group of samples.
#'
#' @author Sur Herrera Paredes
#' 
#' @keywords colonization syncom
#'
#' @export
obtain_block_abundances <- function(Dat, varnames = c("Bacteria","Replicate","Experiment","Pre.Pi", "Pos.Pi"),
                                    taxa.group = "Block", sep = "_", taxa2rm = "contaminant"){
  if(class(Dat) != "Dataset")
    stop("ERROR: Dat must be a Dataset", call. = TRUE)

  if(length(varnames) > length(AMOR::variables(Dat)))
    stop("ERROR: More varnames passed than variables in Dataset", call. = TRUE)

  if(!all(varnames %in% AMOR::variables(Dat)))
    stop("ERROR: Some varnames not in variables in Dataset", call. = TRUE)

  if(length(taxa.group) != 1)
    stop("ERROR: Exactly one taxa variable must be passed", call. = TRUE)

  if(!(taxa.group %in% colnames(Dat$Tax)))
    stop("ERROR: Taxa grouping variable not present in Dataset")

  if("sample_group" %in% AMOR::variables(Dat))
    stop("A variable called 'sample_groups' already exisits in the data.frame", call. = TRUE)

  # Create grouping factor for samples
  sample_groups <- Dat$Map[,varnames]
  sample_groups <- apply(sample_groups,1,paste,collapse = sep)
  Dat$Map[,"sample_group"] <- factor(sample_groups)

  # pool
  abun <- AMOR::pool_samples(Dat,
                             groups = "sample_group",
                             FUN = sum)

  # Combine taxa of the same block
  abun <- AMOR::collapse_by_taxonomy(abun,Group = taxa.group)
  abun <- abun[ -which(row.names(abun) %in% taxa2rm), ]

  # Convert abund to  percent
  abun <- apply(abun,2,function(x) x / sum(x))

  # Transpose
  abun <- as.data.frame(t(abun))
  head(abun)

  # Get metadata
  meta <- metadata_from_rowname(Dat = abun, varnames = varnames, sep = sep)

  # Combine and clean row names
  abun <- cbind(meta, abun)
  row.names(abun) <- NULL

  return(abun)
}

#' Extract metadata from rownames
metadata_from_rowname <- function(Dat, varnames = c("Bacteria","Replicate","Experiment","Pre.Pi", "Pos.Pi"),
                                  sep = "_"){

  meta <- row.names(Dat)
  meta <- limma::strsplit2(x = meta, split = sep)

  if(ncol(meta) != length(varnames)){
    stop("ERROR: Wrong number of columns passed",
         call. = TRUE)
  }

  meta <- as.data.frame(meta)
  colnames(meta) <- varnames

  return(meta)
}
