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

library(wheelP)
library(reshape2)
library(gplots)

# setwd("~/rhizogenomics/experiments/2017/today4/")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")

# Read data
Elongation <- read.table("elongation_single_community_test.txt", header = TRUE, sep = "\t")
Area <- read.table("area_single_community_test.txt", header = TRUE, sep = "\t")
Pi <- read.table("pi_single_community_test.txt", header = TRUE, sep = "\t")
Totalroot <- read.table("totalroot_single_community_test.txt", header = TRUE, sep = "\t")

# Plot
dat <- list(Pi = Pi, Elongation = Elongation,
            Area = Area, Totalroot = Totalroot)
svglite::svglite("syncom_effects_clustered.svg", width = 6, height = 6)
tab <- plot_combined_effects(dat = dat)
rm(Elongation,Area,Pi,Totalroot)

# Read data
Elongation <- read.table("elongation_block_test.txt", header = TRUE, sep = "\t")
Area <- read.table("area_block_test.txt", header = TRUE, sep = "\t")
Pi <- read.table("pi_block_test.txt", header = TRUE, sep = "\t")
Totalroot <- read.table("totalroot_block_test.txt", header = TRUE, sep = "\t")

# Plot
dat <- list(Pi = Pi, Elongation = Elongation,
            Area = Area, Totalroot = Totalroot)
svglite::svglite("block_effects_clustered.svg", width = 6, height = 5)
tab <- plot_combined_effects(dat = dat)
dev.off()


tab <- scale(tab, center = FALSE)
tab <- as.data.frame(tab)
tab$Block <- row.names(tab)
tab <- reshape2::melt(tab,id.vars = "Block")
tab <- cbind(tab, do.call(rbind, strsplit(x = as.character(tab$variable), split = "[.]")))
colnames(tab)[4:5] <- c("Phen","cond")
tab$Block <- factor(tab$Block, levels = c("P1","P2","P3","I1","I2","I3","N1","N2","N3"))

p1 <- ggplot2::ggplot(tab, ggplot2::aes(x =  variable, y = Block)) +
  ggplot2::geom_tile(ggplot2::aes(fill = value)) +
  ggplot2::scale_fill_gradient2(low = c("#8e0152","#de77ae"),
                                min = "#f7f7f7",
                                high = c("#7fbc41","#276419")) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                 axis.text = ggplot2::element_text(color = "black"))
p1
ggplot2::ggsave("block_effects_sorted.svg", p1, width = 5, height = 4)

