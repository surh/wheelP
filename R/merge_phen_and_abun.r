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
