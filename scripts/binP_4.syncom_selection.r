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
library(ggplot2)
library(reshape2)

# setwd("~/rhizogenomics/experiments/2017/today8/")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")
# source("/home/sur/rhizogenomics/src/trunk/phosphate_code/monoP/functions.r")

# Get data
data(binP)
Dat <- binP
rm(binP)

# Get data and reorder strain factor by mean, and groups by effect
Dat <- subset(Dat, Group != "none",drop = TRUE)
Dat$Strain <- factor(Dat$Strain, levels = unique(Dat$Strain))
Dat$Group <- factor(Dat$Group, levels = c("Positive","Indifferent","Negative"))

# Remove un-needed variables
Dat$minusP_100uM.qval <- NULL
Dat$minusP_30uM.qval <- NULL
Dat$plusP_100uM.qval <- NULL
Dat$plusP_30uM.qval <- NULL
Dat$Mean <- NULL

# Reformat
head(Dat)
dat <- melt(data = Dat, varnames = c("minusP_100uM","minusP_30uM","plusP_100uM","plusP_30uM"),
            id.vars = c("Strain","Group"), value.name = "Estimate", variable.name = "Condition")
head(dat)

p1 <- ggplot(dat,aes(x = Strain, y = Condition)) +
  facet_grid(~ Group, scales = "free") +
  geom_tile(aes(fill = Estimate)) +
  # scale_fill_gradient2(low = "#d01c8b",mid = "white",
  #                      high = "#4dac26",midpoint = 0,
  #                      na.value = "#404040",
  #                      guide = guide_colorbar(title = "log(Fold Change)")) +
  scale_fill_gradient2(low = c("#8e0152","#de77ae"),
                       mid = "#f7f7f7", high = c("#7fbc41","#276419"),
                       midpoint = 0,
                       guide = guide_colorbar(title = "log(Fold Change)")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
p1
ggsave("heatmap_syncom_strains_monoP.svg",width = 12, height = 4)
