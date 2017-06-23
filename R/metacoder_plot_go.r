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
#' @author from metacoder paper
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

