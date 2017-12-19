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
# setwd("~/rhizogenomics/experiments/2017/today3/")

data(wheelP.mapsplit)
Dat <- wheelP.mapsplit
rm(wheelP.mapsplit)

Dat <- remove_taxons(Dat, as.character(Dat$Tax$ID[ Dat$Tax$Type == "contaminant" ]))
Dat$Map$Usable <- colSums(Dat$Tab)
Dat <- subset(Dat, Usable > 400)
Dat <- clean(Dat)
Dat <- subset(Dat, !(Experiment %in% c("Validation1","Validation2")),
              drop = TRUE, clean = TRUE)

# set.seed(427)
# Dat <- rarefaction(Dat,sample = 400)

# distfun <- function(x) vegan::vegdist(x,method = "cao")
# Dat.pco <- PCO(Dat,dim = 2, distfun = distfun)
#
#
# plotgg(Dat.pco, col = "Bacteria")
# head(Dat$Map)

cap <- vegan::capscale(t(Dat$Tab) ~ Condition(Plate) + Condition(Usable) +
                         Condition(Experiment),
                       data = Dat$Map,distance = "bray")
cap
cap.sum <- summary(cap)
Dat$Map <- cbind(Dat$Map,cap.sum$sites)
# percvar <- round(100 * cap$CCA$eig / cap$CCA$tot.chi,2)
percvar <- round(100 * cap$CA$eig / cap$CA$real.tot.chi,2)

p1 <- ggplot(Dat$Map, aes(x = MDS1, y = MDS2, col = Fraction, fill = Fraction)) +
  geom_point(aes(shape = Fraction), size = 3) +
  scale_color_manual(values = c("grey35","lawngreen")) +
  scale_fill_manual(values = c("grey35","lawngreen")) +
  scale_shape_manual(values = c(21,22)) +
  xlab(paste("MDS1 (",percvar[1],"%)",sep = "")) +
  ylab(paste("MDS2 (",percvar[2],"%)",sep = "")) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        panel.background = element_rect(color = "black", size = 3, fill = NA),
        panel.grid = element_blank())

p1
ggsave("colonization_mds_conditioned.svg",p1, width = 4, height = 4)
dir.create("figuredata/")
figS6A <- p1$data
save(figS6A, file = "figuredata/figS6A.rda")

