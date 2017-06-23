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

#' Test all strains for change in phosphate accumulation
#'
#' Performs anova in every strain and returns log-fold
#' change versus no bacteria. Performs anova independently
#' in eveyr combination of StartP and EndP levels
#'
#' @param Strains vector of strain names
#' @param treatment_col Column name with variable to test
#' @param treament_control level to use as reference for anova
#' @param Master data.frame with all data
#' @param plot logical produce plot?
#' @param StartP.levels New neames for StartP levels. It must be in the
#' same order as the levels in Master
#' @param EndP.levels New names for EndP levels. Must be in the same order
#' as levels in Master.
#'
#' @return A list.
#'
#' @author Sur Herrera Paredes
#'
#' @export
test_all_strains <- function(Strains, treatment_col,
                             treatment_control = "No Bacteria", Master, plot = FALSE,
                             StartP.levels = c("minusP","plusP"),
                             EndP.levels = c("100uM","30uM")){

  RES <- NULL
  PVAL <- NULL
  for(strain in Strains){

    # Select data
    experiments <- unique(Master$Experiment[Master[,treatment_col] == strain])
    Dat <- Master[ Master[ , treatment_col ] %in% c(strain, treatment_control) &
                     Master$Experiment %in% experiments , ]

    # Define formula
    if(length(unique(Dat$Experiment)) > 1){
      f1 <- formula(log(Pi_content) ~ Bacteria + Experiment)
    }else{
      f1 <- formula(log(Pi_content) ~ Bacteria)
    }

    # Perform test
    res <- independent_lm_tests(Dat = Dat, f1 = f1, strain = strain,
                                StartP.levels = StartP.levels,
                                EndP.levels = EndP.levels)
    RES <- rbind(RES,res$Res)
    PVAL <- rbind(PVAL,res$Pval)

    if(plot){
      filename <- paste(strain,"_mono.png",sep="")
      p1 <- plot_binary(Dat)
      ggsave(filename,p1,width = 10,height = 10)
    }

  }

  # Multiple testing correction
  QVAL <- data.frame(Strain = PVAL$Strain, apply(PVAL[,-1],2,p.adjust,method='fdr'))
  return(list(RES = RES, PVAL = PVAL, QVAL = QVAL))
}
