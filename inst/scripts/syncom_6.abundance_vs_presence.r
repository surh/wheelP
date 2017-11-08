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

setwd("~/rhizogenomics/github/wheelP/")
setwd("~/rhizogenomics/experiments/2017/today4/")

data(Elongation)
data(wheelP.mapsplit)

Dat <- obtain_block_abundances(Dat = wheelP.mapsplit,
                               varnames = c("Bacteria","Replicate","Experiment","Pre.Pi", "Pos.Pi"),
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


head(Dat)
head(Elongation)
