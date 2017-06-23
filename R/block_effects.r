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
#' @export
block_effects <- function(Dat,var.name, cond1,cond2,
                          keep.vars = c("Plate","Experiment")){
  # cond1 <- 1
  # cond2 <- 1
  # keep.vars <- c("Plate","Experiment")
  # var.name <- "Elongation"

  singlecoms <- paste(rep(c("P","I","N"),each = 3),rep(1:3,times = 3),sep="")

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

  dat <- subset(Dat, StartP == levels(Dat$StartP)[cond1] & EndP == levels(Dat$EndP)[cond2])
  dat <- dat[,colnames(dat) %in% c(var.name,singlecoms,keep.vars)]

  f1 <- paste(var.name, " ~ .", sep = "")
  f1 <- formula(f1)
  m1 <- lm(f1 , data = dat )
  #summary(m1)

  return(m1)
}
