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
setwd("~/rhizogenomics/experiments/2017/today7/")
#devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")

# Read data
Master <- read.table("~/rhizogenomics/data/phosphate/exudate_isolate_screening/features/2015-03-31_monomaster.txt",
                     sep="\t",header = TRUE,na.strings = "NA",
                     stringsAsFactors = FALSE)
Master <- Master[!is.na(Master$Pi_content),]

# Fix variables
Master$Replicate[ Master$Replicate == "Rep I" ] <- "RepI"
Master$Replicate[ Master$Replicate == "Rep II" ] <- "RepII"
Master$Replicate[ Master$Replicate == "REpII" ] <- "RepII"
Master$Replicate[ Master$Replicate == "Rep III" ] <- "RepIII"
Master$Replicate[ Master$Replicate == "Rep I, III" ] <- "RepI,III"
Master$Replicate[ Master$Replicate == "Rep I, II" ] <- "RepI,II"
Master$Replicate[ Master$Replicate == "Rep I,III" ] <- "RepI,III"
table(Master$Replicate)
unique(Master$Replicate)

Master$StartP[ Master$StartP == " -Pi,0.5%Suc" ] <- "-Pi,0.5%Suc"
Master$StartP[ Master$StartP == " -P, 0.5%Suc" ] <- "-Pi,0.5%Suc"
Master$StartP[ Master$StartP == " +Pi,0.5%Suc" ] <- "+Pi,0.5%Suc"
Master$StartP[ Master$StartP == " +P, 0.5%Suc" ] <- "+Pi,0.5%Suc"
Master$StartP <- factor(Master$StartP)
unique(Master$StartP)
levels(Master$StartP)

Master$EndP[ Master$EndP == "30 uM,0%Suc"  ] <- "30uM,0%Suc"
Master$EndP[ Master$EndP == "30uM, 0% Suc"  ] <- "30uM,0%Suc"
Master$EndP[ Master$EndP == "100 uM,0%Suc"  ] <- "100uM,0%Suc"
Master$EndP[ Master$EndP == "100 uM, 0%Suc"  ] <- "100uM,0%Suc"
Master$EndP[ Master$EndP == "100uM, 0% Suc"  ] <- "100uM,0%Suc"
Master$EndP[ Master$EndP == "100 uM, 0% Suc"  ] <- "100uM,0%Suc"
Master$EndP <- factor(Master$EndP)
unique(Master$EndP)
levels(Master$EndP)

# Fix names
Master$Treatment[ Master$Treatment == "Cl11" ] <- "CL11"
Master$Treatment[ Master$Treatment == "cl89" ] <- "CL89"
Master$Treatment[ Master$Treatment == "Cl21" ] <- "CL21"
Master$Treatment[ Master$Treatment == "cl28" ] <- "CL28"
Master$Treatment[ Master$Treatment == "cl151" ] <- "CL151"
Master$Treatment[ Master$Treatment == "cl75" ] <- "CL75"
Master$Treatment[ Master$Treatment == "Cl59" ] <- "CL59"
Master$Treatment[ Master$Treatment == "M79_Natalie" ] <- "M79"
Master$Treatment[ Master$Treatment == "PfWSC" ] <- "pf"
Master$Treatment[ Master$Treatment == "R219" ] <- "MR219"
Master$Treatment[ Master$Treatment == "cl130" ] <- "CL130"
Master$Treatment[ Master$Treatment == "DH5alpha" ] <- "Ecoli"
#Master$Treatment[ Master$Treatment == "mR219" ] <- "M219" # Not the same
Master$Treatment[ Master$Treatment == "no bacteria" ] <- "No Bacteria"
Master$Treatment[ Master$Treatment == "no bacteria " ] <- "No Bacteria"
Master$Treatment[ Master$Treatment == "No bacteria" ] <- "No Bacteria"
Master$Treatment[ Master$Treatment == "No bacteria " ] <- "No Bacteria"
Master$Treatment[ Master$Treatment == " no bacteria" ] <- "No Bacteria"
Master$Treatment[ Master$Treatment == " no bacteria " ] <- "No Bacteria"
Master$Treatment <- sub("^m","M",Master$Treatment)

unique(Master$Treatment)[ grep("219",unique(Master$Treatment)) ]
unique(Master$Treatment)
sort(unique(Master$Treatment))
sort(table(Master$Treatment))

# Remove data
Master <- Master[ Master$Treatment != "CL59 + M1",]
#Master <- Master[ Master$Treatment != "M79_Natalie",]
Master$Treatment <- sub(pattern = "^M",replacement = "",x = Master$Treatment)
# sort(unique(Master$Treatment))
Master$Treatment <- factor(Master$Treatment)

# Add main testing variable
Master$Bacteria <- "+Bacteria"
Master$Bacteria[ Master$Treatment == "No Bacteria" ] <- "No Bacteria"
Master$Bacteria <- factor(Master$Bacteria,levels=c("No Bacteria","+Bacteria"))

# Save data in package
binP.all <- Master
rm(Master)
devtools::use_data(binP.all, pkg = "~/rhizogenomics/src/trunk/phosphate_code/wheelP/",
                   overwrite = TRUE)

###########
data(binP.all)

# Plot all strains together
p1 <- plot_binary(binP.all)
p1 <- p1 + geom_vline(xintercept = 10)
p1
ggsave("twoway_allstrains.svg",p1 + geom_vline(xintercept = 10),width = 10, height = 8)
dat <- p1$data

dat$Condition <- NULL
dat$Condition[ dat$StartP == "-Pi,0.5%Suc" & dat$EndP == "30uM,0%Suc"] <- '-P 30uM'
dat$Condition[ dat$StartP == "-Pi,0.5%Suc" & dat$EndP == "100uM,0%Suc"] <- '-P 100uM'
dat$Condition[ dat$StartP == "+Pi,0.5%Suc" & dat$EndP == "30uM,0%Suc"] <- '+P 30uM'
dat$Condition[ dat$StartP == "+Pi,0.5%Suc" & dat$EndP == "100uM,0%Suc"] <- '+P 100uM'
dat$Condition <- factor(dat$Condition, levels = c('-P 100uM','-P 30uM','+P 100uM','+P 30uM'))
dat$Bacteria <- factor(dat$Bacteria, levels = c('+Bacteria','No Bacteria'))
head(dat)
p1 <- ggplot(dat,aes(x = Condition, y = Pi_content, color = Bacteria)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(), size = 0.7, alpha = 0.3) +
  geom_hline(yintercept = 10, color = 'black') +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50)) +
  scale_x_discrete(labels = c(expression(paste('-P 100',mu,M)),
                              expression(paste('-P 30',mu,M)),
                              expression(paste('+P 100',mu,M)),
                              expression(paste('+P 30',mu,M)))) +
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(color = "black",fill = NA, size = 3),
        axis.text.y = element_text(color = 'black', size = 20),
        axis.title.y = element_text(face = 'bold', color = 'black', size = 24),
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.x = element_text(color = "black", size = 20),
        strip.text = element_text(color = "black", size = 32, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))
p1
ggsave('overall_pi_boxplot.svg', p1, width = 10, height = 8)

p1 <- ggplot(dat,aes(x = Condition, y = Pi_content, color = Bacteria, fill = Bacteria)) +
  geom_hline(yintercept = 10, color = 'grey', size = 3) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), alpha = 0.4) +
  #geom_point(position = position_jitterdodge(), size = 0.7, alpha = 0.3) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50)) +
  scale_x_discrete(labels = c(expression(paste('-P 100',mu,M)),
                              expression(paste('-P 30',mu,M)),
                              expression(paste('+P 100',mu,M)),
                              expression(paste('+P 30',mu,M)))) +
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(color = "black",fill = NA, size = 3),
        axis.text.y = element_text(color = 'black', size = 20),
        axis.title.y = element_text(face = 'bold', color = 'black', size = 24),
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.x = element_text(color = "black", size = 20),
        strip.text = element_text(color = "black", size = 32, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))
p1
ggsave('overall_pi_violin.svg', p1, width = 10, height = 8)

# Test overall bacteril effect
m1 <- aov(log(Pi_content) ~ Bacteria * StartP * EndP + Experiment ,data = binP.all)
summary(m1)
drop1(m1)

m1 <- aov(log(Pi_content) ~ Bacteria  + StartP + EndP + Experiment +
            StartP * Bacteria + EndP * Bacteria + StartP * EndP,data = binP.all)
m1.sum <- summary(m1)
m1.sum
drop1(m1)
TukeyHSD(m1)

m1 <- lm(log(Pi_content) ~ Bacteria  + StartP + EndP + Experiment +
           StartP * Bacteria + EndP * Bacteria + StartP * EndP,data = binP.all)
m1.sum <- summary(m1)
m1.sum
m1.sum$r.squared
m1.sum$adj.r.squared
