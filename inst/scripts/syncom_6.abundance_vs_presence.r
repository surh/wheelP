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

library(AMOR)
library(wheelP)

# setwd("~/rhizogenomics/github/wheelP/")
# devtools::document()
setwd("~/rhizogenomics/experiments/2017/today5/")
outdir <- "~/rhizogenomics/experiments/2017/today5/"

data(Elongation)
data(wheelP.mapsplit)

Dat <- obtain_block_abundances(Dat = wheelP.mapsplit,
                               varnames = c("Bacteria","Experiment","Pre.Pi", "Pos.Pi"),
                               taxa.group = "Block", sep = "_", taxa2rm = "contaminant")
head(Dat)
head(Elongation)

# Remove validation experiments since there is no pohenotypic data
Dat <- droplevels(subset(Dat, Experiment %in% c("1","2")))

# Rename experiments to match real batches
Dat$Bacteria <- factor(Dat$Bacteria, levels = levels(Elongation$Bacteria))
ftable(Experiment ~ Bacteria, data = Elongation)
ftable(Experiment ~ Bacteria, data = Dat)
Dat$EXP <- NULL
Dat$EXP[ Dat$Bacteria %in% c("P1P2","P2P3","P3I1","I1I2","I2I3") & Dat$Experiment == "1" ] <- "A"
Dat$EXP[ Dat$Bacteria %in% c("P1P2","P2P3","P3I1","I1I2","I2I3") & Dat$Experiment == "2" ] <- "B"
Dat$EXP[ Dat$Bacteria %in% c("I3N1","N1N2","N2N3","P1N3") & Dat$Experiment == "1" ] <- "C"
Dat$EXP[ Dat$Bacteria %in% c("I3N1","N1N2","N2N3","P1N3") & Dat$Experiment == "2" ] <- "D"
Dat$EXP[ Dat$Bacteria %in% c("P1I1","P2I1","P1P3") & Dat$Experiment == "1" ] <- "F"
Dat$EXP[ Dat$Bacteria %in% c("P1I1","P2I1","P1P3") & Dat$Experiment == "2" ] <- "H"
Dat$EXP[ Dat$Bacteria %in% c("P2N3","P3N3") & Dat$Experiment == "1" ] <- "E"
Dat$EXP[ Dat$Bacteria %in% c("P2N3","P3N3") & Dat$Experiment == "2" ] <- "G"
ftable(EXP ~ Bacteria, data = Dat)
Dat$Experiment <- Dat$EXP
Dat$EXP <- NULL

# Rename condition names
Dat$StartP <- Dat$Pre.Pi
Dat$EndP <- Dat$Pos.Pi
Dat$Pre.Pi <- NULL
Dat$Pos.Pi <- NULL

# Collapse elongation to match sequencing data
Phen <- Elongation
Phen <- aggregate(Elongation ~ Bacteria + Experiment + StartP + EndP, data = Phen, FUN = mean)
head(Phen)

# Test phenotype on collapsed data
Res.sc <- test_single_community_phenotype(Dat = Phen,dir = outdir,
                                          var.name = "Elongation",
                                          bacteria.col = "Bacteria", ref.level = 'none',
                                          plot = FALSE, f1.extra = "+ Experiment")
summary(qvalue::qvalue(Res.sc$p.value))

# Now test block effects
singlecoms <- paste(rep(c("P","I","N"),each = 3),rep(1:3,times = 3),sep="")
Phen <- merge_phen_and_abun(Phen = Phen, abun = Dat, columns = singlecoms,)
head(Phen)

Res2 <- NULL
for(i in 1:2){
  for(j in 1:2){
    m1 <- block_effects(Dat = Dat, cond1 = i, cond2 = j,
                        var.name = "Elongation",
                        keep.vars = c("Plate","Experiment"))
    m1.sum <- summary(m1)

    res <- data.frame(SynCom = singlecoms, StartP = levels(Dat$StartP)[i],
                      EndP = levels(Dat$EndP)[j],
                      Estimate = m1.sum$coefficients[ singlecoms, 1 ],
                      SE = m1.sum$coefficients[ singlecoms, 2 ],
                      t.value = m1.sum$coefficients[ singlecoms, 3 ],
                      p.value = m1.sum$coefficients[ singlecoms, 4 ])
    Res2 <- rbind(Res2,res)
    # p1 <- coefplot(m1, intercept = FALSE,innerCI = 1, outerCI = 2,
    #                coefficients = singlecoms,
    #                lwdOuter = 0.5, lwdInner = 2.5, pointSize = 4,
    #                color = "black", zeroColor = "red", zeroType = 1)
    # p1 <- p1 + theme(axis.text.y = element_text(color = "black"),
    #                  panel.border = element_rect(fill = NA, color = "black", size = 2),
    #                  panel.background = element_rect(fill = "white")) +
    #   ggtitle(paste(levels(Dat$StartP)[i],"=>",levels(Dat$EndP)[j]))
    # p1
    #
    # filename <- paste(dir,"/coefplot_",i,j,".png",sep = "")
    # ggsave(filename = filename, p1 , width = 3.5, height = 5)
  }
}
