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

#' Test all combinations of levels
#'
#' Internal function.
#' 
#' This functions takes the data from one bacterial treatment makes
#' an independent linear model for all observation on each combination 
#' of levels between two variables.
#'
#' @param Dat A data.frame containing the data of one bacterial treatment.
#' Unsupported behaviour would occur if more than one bacterial treatment is
#' included. The data.frame must contain columns: 'Bacteria', 'StartP' and
#' 'EndP'. The 'Bacteria' column must be either a factor or a character vector
#' containing only values 'No Bacteria' and '+Bacteria'
#' @param f1 formula for test. It will be passed to the lm function
#' @param strain A character string, indicating strain name.
#' @param StartP.levels new names for StartP levels
#' @param EndP.lvels new names for EndP levels.
#'
#' @return A list of data frames. Each data.frame contains a column called 'Strain'
#' which indicates the strain ID. An also contains m more columns
#' (m = length(StartP.levels)*lenght(EndP.levels)). Each element of the list stores the
#' either the coefficiennt for the bacterial treatment (Res), its p-value (Pval), or
#' its standard error (SE).
#' 
#' @seealso \code{\link{lm}}
#' 
#' @keywords utilities syncom binP
#'
#' @author Sur Herrera Paredes
independent_lm_tests <- function(Dat,f1,strain,StartP.levels = c("minusP","plusP"),
                                 EndP.levels = c("100uM","30uM")){
  levels(Dat$StartP) <- StartP.levels
  levels(Dat$EndP) <- EndP.levels
  Dat$Bacteria <- factor(Dat$Bacteria,levels=c("No Bacteria","+Bacteria"))

  Models <- list()
  for(startp in StartP.levels){
    for(endp in EndP.levels){
      name <- paste(startp,endp,sep = "_")
      Models[[name]] <- summary(lm(f1, data = subset(Dat,StartP == startp & EndP == endp )))
    }
  }

  Res <- data.frame(Strain = strain)
  Pval <- data.frame(Strain = strain)
  SE <- data.frame(Strain = strain)
  for(name in names(Models)){
    Res[,name] <- Models[[name]]$coefficients["Bacteria+Bacteria",1]
    Pval[,name] <- Models[[name]]$coefficients["Bacteria+Bacteria",4]
    SE[,name] <- Models[[name]]$coefficients["Bacteria+Bacteria",2]
  }

  #   m2.sum <- summary(lm(f1, data = subset(Dat,StartP == "minusP" & EndP == "100uM" )))
  #   m3.sum <- summary(lm(f1, data = subset(Dat,StartP == "minusP" & EndP == "30uM" )))
  #   m4.sum <- summary(lm(f1, data = subset(Dat,StartP == "plusP" & EndP == "100uM" )))
  #   m5.sum <- summary(lm(f1, data = subset(Dat,StartP == "plusP" & EndP == "30uM" )))
  #   Res <- data.frame(Strain = strain, minusP_100uM = m2.sum$coefficients["Bacteria+Bacteria",1],
  #                     minusP_30uM = m3.sum$coefficients["Bacteria+Bacteria",1],
  #                     plusP_100uM = m4.sum$coefficients["Bacteria+Bacteria",1],
  #                     plusP_30uM = m5.sum$coefficients["Bacteria+Bacteria",1])
  #
  #   Pval <- data.frame(Strain = strain, minusP_100uM = m2.sum$coefficients["Bacteria+Bacteria",4],
  #                      minusP_30uM = m3.sum$coefficients["Bacteria+Bacteria",4],
  #                      plusP_100uM = m4.sum$coefficients["Bacteria+Bacteria",4],
  #                      plusP_30uM = m5.sum$coefficients["Bacteria+Bacteria",4])

  return(list(Res = Res, SE = SE, Pval = Pval))
}
