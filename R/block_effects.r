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

#' Estimate block main effects
#' 
#' Fits a single linear model for all data in a specific condition and
#' estimates the main effect of each abcterial block
#' 
#' This function will create a design matrix based on SynCom names. Therefore,
#' it assumes that the SynCom names are named \emph{B#B#}, where 'B' is replaced
#' by the functiona group type ('P', "I' or 'N') and # is the block number
#' within each functional type. It assumes that exactly nine blocks (3 per functional
#' type) exist.
#' 
#' The function will also include any variable named in keep.vars as a covariate
#' in the linear models.
#' 
#' Only the data that matches where Dat$StartP == cond1 & Dat$EndP == cond2 will
#' be used to fit the model. The rest is discarded.
#' 
#' @param Dat A data.frame. Must contain the following columns: 'Bacteria',
#' 'StartP' and 'EndP'. It must further contain a column matching the var.name
#' and keep.var parameters.
#' @param var.name A string character indicating the name of the variable in
#' Dat that has the lef-hand side of the linear model (i.e. the dependent variable)
#' @param cond1 The value for the 'StartP' variable in Dat. Only data where
#' Dat$StartP == cond1 will be kept.
#' @param cond2 The value for the 'EndP' variable in Dat. Only data where
#' Dat$EndP == cond2 will be kept.
#' @param keep.vars A vector of character objects indicating which variables from
#' Dat should be kept as covariates in the model
#' @param create.design logcial indicats whether the desing must be created for
#' the main effects. If FALSE, a linear model with the variables indicated
#' by keep.vars will be produced.
#' 
#' @return Returns an object of class 'lm', see \code{\link{lm}} for more info
#' 
#' @seealso \code{\link{lm}}
#' 
#' @keywords syncom
#' 
#' @author Sur Herrera Paredes
#'
#' @export
block_effects <- function(Dat,var.name, cond1,cond2,
                          keep.vars = c("Plate","Experiment"),
                          create.design = TRUE){

  singlecoms <- paste(rep(c("P","I","N"),each = 3),rep(1:3,times = 3),sep="")

  if(create.design){
    # Build design matrix
    X <- matrix(0, nrow = nrow(Dat), ncol = 10 )
    colnames(X) <- c("Inoculated", singlecoms)
    X[,"Inoculated"] <- 1*(Dat$Bacteria != "none")
    X[ grep(pattern = "P1",x = Dat$Bacteria), "P1"] <- 1
    X[ grep(pattern = "P2",x = Dat$Bacteria), "P2"] <- 1
    X[ grep(pattern = "P3",x = Dat$Bacteria), "P3"] <- 1
    X[ grep(pattern = "I1",x = Dat$Bacteria), "I1"] <- 1
    X[ grep(pattern = "I2",x = Dat$Bacteria), "I2"] <- 1
    X[ grep(pattern = "I3",x = Dat$Bacteria), "I3"] <- 1
    X[ grep(pattern = "N1",x = Dat$Bacteria), "N1"] <- 1
    X[ grep(pattern = "N2",x = Dat$Bacteria), "N2"] <- 1
    X[ grep(pattern = "N3",x = Dat$Bacteria), "N3"] <- 1

    Dat <- cbind(X,Dat)
  }

  dat <- subset(Dat, StartP == levels(Dat$StartP)[cond1] & EndP == levels(Dat$EndP)[cond2])
  dat <- dat[,colnames(dat) %in% c(var.name,singlecoms,keep.vars)]

  f1 <- paste(var.name, " ~ .", sep = "")
  f1 <- formula(f1)
  m1 <- lm(f1 , data = dat )
  
  return(m1)
}
