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
Dat <- subset(Dat, (Experiment %in% c("Validation1","Validation2")),
              drop = TRUE, clean = TRUE)
Dat <- subset(Dat, Fraction != "T0Agar",
              drop = TRUE, clean = TRUE)

Dat.ord <- collapse_by_taxonomy(Dat,level = 5, sepchar = ";",FUN = sum)

temp <- subset(Dat.ord, Fraction == "Root")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 8) +
  scale_fill_brewer(palette = "Paired") +
  ggtitle(label = "Root") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12, angle = 90),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("validation.phylogram.svg",p1, width = 14, height = 4)

temp <- subset(Dat.ord, Fraction == "Agar")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 8) +
  scale_fill_brewer(palette = "Paired") +
  ggtitle(label = "Agar") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12, angle = 90),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("validation.phylogram.agar.svg",p1, width = 14, height = 4)

temp <- subset(Dat.ord, Fraction == "Inoculum")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 8) +
  scale_fill_brewer(palette = "Paired") +
  ggtitle(label = "Inoculum") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12,angle = 90),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("validation.phylogram.inoculum.svg",p1, width = 14, height = 4)

rm(Dat.ord, p1, temp)


#### Blocks
Blocks <- create_dataset(Tab = collapse_matrix(Dat$Tab,
                                               groups = Dat$Tax$Block,
                                               FUN = sum,
                                               dim = 1),
                         Map = Dat$Map[,c("ID","Bacteria","Pre.Pi",
                                          "Pos.Pi","Fraction")])

temp <- subset(Blocks, Fraction == "Root")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 9) +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 9,name = "PiYG"))) +
  ggtitle(label = "Root") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12, angle = 90),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("validation.phylogram.blocks.svg",p1, width = 14, height = 4)

temp <- subset(Blocks, Fraction == "Agar")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 9) +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 9,name = "PiYG"))) +
  ggtitle(label = "Agar") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12, angle = 90),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("validation.phylogram.blocks.agar.svg",p1, width = 14, height = 4)

temp <- subset(Blocks, Fraction == "Inoculum")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 9) +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 9,name = "PiYG"))) +
  ggtitle(label = "Inoculum") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12, angle = 90),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("validation.phylogram.blocks.inoculum.svg",p1, width = 14, height = 4)
