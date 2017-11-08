#' Obtain block abundances
#'
#' Takes a Dataset object, pools samples and taxa according to
#' specified variables, and returns  adata frame that contains
#' the average abundance by taxa group per group of samples
#'
#' @author Sur Herrera Paredes
#'
#' @export
obtain_block_abundances <- function(Dat, varnames = c("Bacteria","Replicate","Experiment","Pre.Pi", "Pos.Pi"),
                                    taxa.group = "Block", sep = "_", taxa2rm = "contaminant"){
  if(class(Dat) != "Dataset")
    stop("ERROR: Dat must be a Dataset", call. = TRUE)

  if(length(varnames) > length(AMOR::variables(Dat)))
    stop("ERROR: More varnames passed than variables in Dataset", call. = TRUE)

  if(!all(varnames %in% variables(Dat)))
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
  meta <- strsplit2(x = meta, split = sep)

  if(ncol(meta) != length(varnames)){
    stop("ERROR: Wrong number of columns passed",
         call. = TRUE)
  }

  meta <- as.data.frame(meta)
  colnames(meta) <- varnames

  return(meta)
}
