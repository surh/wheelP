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

#' Read count tables from mapsplit dirs
#'
#' @export
read_mapsplit_dir <- function(dir){

  # dir <- "tables_mapsplit_trim350_0.985/160201/"

  subdirs <- list.dirs(dir)
  subdirs <- subdirs[-1]
  DAT <- NULL
  for(i in 1:length(subdirs)){
    #i <- 1

    subdir <- subdirs[i]

    filename <- paste(subdir,"/","otu_table.txt",sep = "")
    dat <- AMOR::read.am(filename,format = "qiime", taxonomy = "taxonomy")

    if(i == 1){
      DAT <- dat
    }else{
      DAT <- combine_datasets(list(DAT,dat))
    }
  }

  return(DAT)
}
