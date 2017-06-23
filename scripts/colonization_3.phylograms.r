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

Dat.ord <- collapse_by_taxonomy(Dat,level = 5, sepchar = ";",FUN = sum)

temp <- subset(Dat.ord, Pre.Pi == "+Pi,0.5%Suc" &
                 Pos.Pi == "100 uM,0%Suc" &
                 Fraction == "Root")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 8) +
  scale_fill_brewer(palette = "Paired") +
  ggtitle(label = "+Pi,0.5%Suc => 100 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("plusP_100uM.phylogram.svg",p1, width = 14, height = 4)

temp <- subset(Dat.ord, Pre.Pi == "+Pi,0.5%Suc" &
                 Pos.Pi == "100 uM,0%Suc" &
                 Fraction == "Agar")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 8) +
  scale_fill_brewer(palette = "Paired") +
  ggtitle(label = "+Pi,0.5%Suc => 100 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("plusP_100uM.phylogram.agar.svg",p1, width = 14, height = 4)

temp <- subset(Dat.ord, Pre.Pi == "+Pi,0.5%Suc" &
                 Pos.Pi == "30 uM,0%Suc" &
                 Fraction == "Root")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 8) +
  scale_fill_brewer(palette = "Paired") +
  ggtitle(label = "+Pi,0.5%Suc => 30 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("plusP_30uM.phylogram.svg",p1, width = 14, height = 4)

temp <- subset(Dat.ord, Pre.Pi == "+Pi,0.5%Suc" &
                 Pos.Pi == "30 uM,0%Suc" &
                 Fraction == "Agar")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 8) +
  scale_fill_brewer(palette = "Paired") +
  ggtitle(label = "+Pi,0.5%Suc => 30 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("plusP_30uM.phylogram.agar.svg",p1, width = 14, height = 4)

temp <- subset(Dat.ord, Pre.Pi == "-Pi,0.5%Suc" &
                 Pos.Pi == "100 uM,0%Suc" &
                 Fraction == "Root")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 8) +
  scale_fill_brewer(palette = "Paired") +
  ggtitle(label = "-Pi,0.5%Suc => 100 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("minusP_100uM.phylogram.svg",p1, width = 14, height = 4)

temp <- subset(Dat.ord, Pre.Pi == "-Pi,0.5%Suc" &
                 Pos.Pi == "100 uM,0%Suc" &
                 Fraction == "Agar")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 8) +
  scale_fill_brewer(palette = "Paired") +
  ggtitle(label = "-Pi,0.5%Suc => 100 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("minusP_100uM.phylogram.agar.svg",p1, width = 14, height = 4)

temp <- subset(Dat.ord, Pre.Pi == "-Pi,0.5%Suc" &
                 Pos.Pi == "30 uM,0%Suc" &
                 Fraction == "Root")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 8,space = "fixed") +
  scale_fill_brewer(palette = "Paired") +
  ggtitle(label = "-Pi,0.5%Suc => 30 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("minusP_30uM.phylogram.svg",p1, width = 14, height = 4)

temp <- subset(Dat.ord, Pre.Pi == "-Pi,0.5%Suc" &
                 Pos.Pi == "30 uM,0%Suc" &
                 Fraction == "Agar")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 8,space = "fixed") +
  scale_fill_brewer(palette = "Paired") +
  ggtitle(label = "-Pi,0.5%Suc => 30 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("minusP_30uM.phylogram.agar.svg",p1, width = 14, height = 4)

rm(Dat.ord, p1, temp)

## Theoretical
data(Tax.colonization)
Tax <- subset(Tax.colonization, Type != "contaminant")
syncoms <- levels(Dat$Map$Bacteria)

Tab <- matrix(0, nrow = nrow(Tax), ncol = length(syncoms))
row.names(Tab) <- row.names(Tax)
colnames(Tab) <- syncoms
for(syncom in syncoms){
  # syncom <- syncoms[1]

  b1 <- substr(syncom,1,2)
  b2 <- substr(syncom,3,4)

  Tab[ row.names(subset(Tax,Block %in% c(b1,b2))), syncom ] <- 1
}
Theo <- create_dataset(Tab = Tab,
                       Map = data.frame(ID = syncoms,
                                        Bacteria = factor(syncoms,
                                                           levels = levels(Dat$Map$Bacteria)),
                                        row.names = syncoms),
                       Tax = Tax)
Theo.ord <- collapse_by_taxonomy(Theo,level = 5, sepchar = ";",FUN = sum)

p1 <- phylogram(Theo.ord,facet = ~ Bacteria,nrow.legend = 8) +
  scale_fill_brewer(palette = "Paired") +
  ggtitle(label = "-Pi,0.5%Suc => 30 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("theoretical.phylogram.svg",p1, width = 14, height = 4)

rm(Tab,Tax,Tax.colonization,b1,b2,p1,syncom,syncoms,Theo,Theo.ord)

### Blocks
Blocks <- create_dataset(Tab = collapse_matrix(Dat$Tab,
                                               groups = Dat$Tax$Block,
                                               FUN = sum,
                                               dim = 1),
                         Map = Dat$Map[,c("ID","Bacteria","Pre.Pi","Pos.Pi","Fraction")])

temp <- subset(Blocks, Pre.Pi == "+Pi,0.5%Suc" &
                 Pos.Pi == "100 uM,0%Suc" &
                 Fraction == "Root")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 9) +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 9,name = "PiYG"))) +
  ggtitle(label = "+Pi,0.5%Suc => 100 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("plusP_100uM.phylogram.blocks.svg",p1, width = 14, height = 4)

temp <- subset(Blocks, Pre.Pi == "+Pi,0.5%Suc" &
                 Pos.Pi == "100 uM,0%Suc" &
                 Fraction == "Agar")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 9) +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 9,name = "PiYG"))) +
  ggtitle(label = "+Pi,0.5%Suc => 100 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("plusP_100uM.phylogram.blocks.agar.svg",p1, width = 14, height = 4)

temp <- subset(Blocks, Pre.Pi == "+Pi,0.5%Suc" &
                 Pos.Pi == "30 uM,0%Suc" &
                 Fraction == "Root")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 9) +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 9,name = "PiYG"))) +
  ggtitle(label = "+Pi,0.5%Suc => 30 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("plusP_30uM.phylogram.blocks.svg",p1, width = 14, height = 4)

temp <- subset(Blocks, Pre.Pi == "+Pi,0.5%Suc" &
                 Pos.Pi == "30 uM,0%Suc" &
                 Fraction == "Agar")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 9) +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 9,name = "PiYG"))) +
  ggtitle(label = "+Pi,0.5%Suc => 30 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("plusP_30uM.phylogram.blocks.agar.svg",p1, width = 14, height = 4)

temp <- subset(Blocks, Pre.Pi == "-Pi,0.5%Suc" &
                 Pos.Pi == "100 uM,0%Suc" &
                 Fraction == "Root")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 9) +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 9,name = "PiYG"))) +
  ggtitle(label = "-Pi,0.5%Suc => 100 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("minusP_100uM.phylogram.blocks.svg",p1, width = 14, height = 4)

temp <- subset(Blocks, Pre.Pi == "-Pi,0.5%Suc" &
                 Pos.Pi == "100 uM,0%Suc" &
                 Fraction == "Agar")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 9) +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 9,name = "PiYG"))) +
  ggtitle(label = "+Pi,0.5%Suc => 100 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("minusP_100uM.phylogram.blocks.agar.svg",p1, width = 14, height = 4)

temp <- subset(Blocks, Pre.Pi == "-Pi,0.5%Suc" &
                 Pos.Pi == "30 uM,0%Suc" &
                 Fraction == "Root")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 9) +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 9,name = "PiYG"))) +
  ggtitle(label = "-Pi,0.5%Suc => 30 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("minusP_30uM.phylogram.blocks.svg",p1, width = 14, height = 4)

temp <- subset(Blocks, Pre.Pi == "-Pi,0.5%Suc" &
                 Pos.Pi == "30 uM,0%Suc" &
                 Fraction == "Agar")
p1 <- phylogram(temp,facet = ~ Bacteria,nrow.legend = 9) +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 9,name = "PiYG"))) +
  ggtitle(label = "-Pi,0.5%Suc => 30 uM,0%Suc") +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 14, color = "black", angle = 0),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank())
p1
ggsave("minusP_30uM.phylogram.blocks.agar.svg",p1, width = 14, height = 4)

rm(p1,Blocks,temp)
