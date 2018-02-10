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

#' Plot GO from edgeR
#' 
#' Creates a metacoder network plot from gene ontologies. Figures from this function
#' where not included in final manuscript.
#' 
#' The script takes a list of differentially expressed genes and their GO annotations
#' and creates a tree representing the gene ontologies present and color codes it by
#' the overal log(fold-change) im expression/
#' 
#' @param dat A data frame containing the results of differential expression analysis on
#' edgeR. The data frame is the result of performing hypothesis testing, for example
#' with the likelihood ratio test (i.e. function \code{\link{glmLRT}}).
#' @param output_folder Charachter string indicating the directory to store the
#' plots.
#' @param output_format Format of the image to save. It must be allowed
#' by the \code{\link{heat_tree}} function.
#' @param min_fdr False discovery rate threshold to filter genes. Only genes
#' with a smaller fdr and their annotations will be included.
#' @param type This must be an object that specifies the parents of each
#' annotation terms. Default is object GOBPPARENTS fromthe GO.db bioconductor
#' package, but any other simular object can be used.
#' @param prefix Character string. Prefix for file names.
#' @param n.supertaxa Number of levels to extend the tree.
#' @param num.changed Minimum number of genes in an annotation category
#' fir that annotation to be included.
#' 
#' @return The data used by \code{link{heat_tree}} function.
#'
#' @author Code from metacoder paper with small adaptations by Sur Herrera Paredes
#' 
#' @references Foster ZSL, Sharpton TJ, Grünwald NJ (2017) Metacoder:
#' An R package for visualization and manipulation of community taxonomic diversity data.
#' PLoS Comput Biol 13(2): e1005404. https://doi.org/10.1371/journal.pcbi.1005404
#' 
#' @keywords rna plots
#' 
#' @seealso \code{\link{parse_tax_and_plot}} \code{\link{heat_tree}}
#'
#' @export
metacoder_plot_go <- function(dat, output_folder = "./", output_format = "svg",
                              min_fdr = 0.01, type = GO.db::GOBPPARENTS, prefix = "go",
                              n.supertaxa = 9, num.changed = 3){
  # dat <- Res
  # output_folder <- "bacteria_vs_nobac/"
  # output_format <- "svg"
  # min_fdr <-0.000001
  # type <- GO.db::GOBPPARENTS
  # prefix <- "gobp"


  # Prepare otuput
  dir.create(output_folder)
  output_file <- paste(output_folder,"/",prefix,"_pub.",output_format,sep = "")
  go_res_file <- paste(output_folder,"/",prefix,"_go_res.txt",sep = "")

  # Keep only significant
  dat <- subset(dat,FDR < min_fdr)

  # Get GO terms
  # dat$go <- AnnotationDbi::mapIds(org.At.tair.db,
  #                                 keys = dat$Gene,
  #                                 column = "GO",
  #                                 keytype = "TAIR",
  #                                 multiVals = "first")

  gos <- AnnotationDbi::mapIds(org.At.tair.db::org.At.tair.db,
                               keys = dat$Gene,
                               column = "GO",
                               keytype = "TAIR",
                               multiVals = "list")
  times <- sapply(gos,length)
  dat <- data.frame(logFC = rep(dat$logFC,times = times),
                    Gene = rep(dat$Gene, times = times),
                    FDR = rep(dat$FDR, times = times),
                    go = unlist(gos), row.names = NULL,
                    stringsAsFactors = FALSE)
  dat <- dat[!is.na(dat$go), ]

  # Get parents
  go_class <- lapply(dat$go, wheelP::term_class, all_paths = FALSE, type = type)

  ## Calculate parent tree and write file.
  go_res <- dat[rep(1:nrow(dat), sapply(go_class, length)), ]
  go_res$class <- unlist(go_class)
  write.table(x = go_res, file = go_res_file, sep = "\t",
              quote = FALSE,row.names = FALSE)
  col <- which(colnames(go_res) == "class")

  data <- wheelP::parse_tax_and_plot(file = go_res_file, col = col,
                                     output_file = output_file,
                                     n.supertaxa = n.supertaxa,
                                     num.changed = num.changed,
                                     min_fdr = min_fdr)

  return(data)
}

