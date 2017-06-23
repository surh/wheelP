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

#' Bacterial growth curves in 4 different media
#' 
#' A data.frame containing the 72-hour long of 441 strains
#' grown in 4 conditions. Growth curves have been filtered
#' to remove strains that didn't grow by clustering and removing
#' those strains that looked more similar to sample blanks
#' 
#' @format A data.frame with the following columns
#' \describe{
#'   \item{Row}{Row letter in the 96-well plate wee bacteria was grown}
#'   \item{Column}{Column from the 96-well plate were bacteria was grown}
#'   \item{OD600}{Bacterial abundance measurment in optical density at 600nm}
#'   \item{Rep}{Replicate number}
#'   \item{hrs}{Time-point of measurment in hours}
#'   \item{Strain}{Strain ID}
#'   \item{Condition}{Media used to grow the bacteria. minusP is Johnson
#'   media with 0 phosphate added, plusP is Johnson media with 1mM phosphate
#'   added. For minus2plusP, Arabidopsis seeds were germinated in Johnson
#'   media with zero phosphate, and then transferred for 24 hours into
#'   Johnson media with 1mM phosphate. After 24h the media was filter-sterilized
#'   and bacteria were grown on it. For plus2minusP Arabidopsis seeds were
#'   germinated in Johnson media with 1mM phosphate, and then transferred
#'   for 24 hours into Johnson media with no phosphate. After 24h the
#'   media was filter-sterilized}
#' }
#' 
#' @source Dangl Lab
"All.filtered"
