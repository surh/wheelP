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

#' Test AUC
#' 
#' Tests area under the curve for a set of
#' strains
#' 
#' Assumes same number of timepoints per growth curve
#' 
#' @param Dat A data.frame with the same columns as \link{All.filtered}
#' @param thres q-value threshold for significance
#' 
#' @export
test_auc <- function(Dat, thres = 0.05){
  # Dat <- All.filtered
  RES <- NULL
  for(strain in unique(Dat$Strain)){
    #strain <- unique(Dat$Strain)[1]
    # strain <- "113"
    
    # Get strain data and calculate AUD
    Dat.strain <- subset(Dat, Strain == strain)
    Dat.strain <- aggregate(OD600 ~ rep + condition,
                            FUN = sum, data = Dat.strain)
    
    f1 <- OD600 ~ condition
    
    test <- TukeyHSD(aov(f1, data = Dat.strain))$condition
    Res <- data.frame(Strain = strain,
                      comparison = row.names(test),
                      difference = test[,"diff"], qval = test[,"p adj"])
    row.names(Res) <- NULL
    
    RES <- rbind(RES,Res)
  }
  
  RES$significant <- "no"
  RES$significant[ RES$qval < thres ] <- "yes"
  
  return(RES)
}

