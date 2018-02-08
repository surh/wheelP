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

#' Distance to data.frame
#' 
#' Convets an object of class "dist" into an
#' obejct of class "data.frame" where each row is
#' a pairwise distance
#' 
#' @param inDist An object of class "dist". See \code{\link{dist}} for
#' more info
#' 
#' @author Sur Herrera Paredes
#' 
#' @seealso \code{\link{dist}}
#' 
#' @keywords utilities
#'
#' @export
dist2df <- function(inDist){
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  if(is.null(attr(inDist, "Labels"))){
    B <- sequence(A)
  }else{
    B <- attr(inDist, "Labels")
  }
  attr(inDist, "Diag") <- FALSE
  attr(inDist, "Upper") <- FALSE

  res <- data.frame(row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                    col = rep(B[-length(B)], (length(B)-1):1),
                    value = as.vector(inDist))

  return(res)
}
