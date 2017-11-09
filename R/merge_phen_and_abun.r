#' Merge phenotype data with relative abundances of blocks
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
