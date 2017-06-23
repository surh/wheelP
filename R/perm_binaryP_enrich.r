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

#' Enrichment of phosphate modulators
#'
#' @param test result from test_all_strains
#' @param N number of permutations
#' @param qval.thres q-valuethreshold
#' @param max.val Maximum value for log fold chance
#' @param min.val maximum value for log fold change
#'
#' @return vector of p-values
#'
#' @author Sur Herrera Paredes
#'
#' @export
perm_binaryP_enrich <- function(test, N = 999, qval.thres = 0.1, max.val = Inf, min.val = -Inf){
  tab <- test$QVAL[,2:5] < qval.thres & test$RES[,2:5] < max.val & test$RES[,2:5] > min.val

  obs <- colSums(tab)
  pval <- as.numeric(obs <= obs)
  for(i in 1:999){
    perm <- apply(tab,1,sample)
    perm <- rowSums(perm)
    pval <- pval + (obs <= perm)
  }
  pval <- pval / 1000

  return(pval)
}
