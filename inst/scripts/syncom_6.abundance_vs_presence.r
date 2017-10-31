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

setwd("~/rhizogenomics/experiments/2017/today4/")

data(Elongation)
data(wheelP.mapsplit)


# Remove validation experiment for which there is only one phenotype
wheelP.mapsplit <- subset(wheelP.mapsplit,
                          !(Experiment %in% c("Validation1","Validation2")),
                          drop = TRUE,clean = TRUE)





# head(Elongation)
# head(wheelP.mapsplit$Map)
#
# Elongation$Experiment
# wheelP.mapsplit$Map$Experiment
#
# Elongation$Bacteria
# wheelP.mapsplit$Map$Bacteria


# Combine samples of the same group
wheelP.mapsplit$Map$Group <- paste(wheelP.mapsplit$Map$Bacteria,
                                   wheelP.mapsplit$Map$Replicate,
                                   wheelP.mapsplit$Map$Experiment,
                                   wheelP.mapsplit$Map$Pre.Pi,
                                   wheelP.mapsplit$Map$Pos.Pi, sep = "_")
abun <- pool_samples(wheelP.mapsplit,
                     groups = "Group",
                     FUN = sum)

# Combine taxa of the same block
abun <- collapse_by_taxonomy(abun,Group = "Block")
abun <- abun[ -which(row.names(abun) == "contaminant"), ]

# Convert abund to  percent
abun <- apply(abun,2,function(x) x / sum(x))


abun <- as.data.frame(t(abun))
head(abun)

# Reformat Abundances
abun <- melt(abun)




