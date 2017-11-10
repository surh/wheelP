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
library(ggplot2)

# setwd("~/rhizogenomics/experiments/2017/today7/phenotypes/")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")

# Read data
# Elongation <- read.table("elongation_single_community_test.txt", header = TRUE, sep = "\t")
# Area <- read.table("area_single_community_test.txt", header = TRUE, sep = "\t")
# Pi <- read.table("pi_single_community_test.txt", header = TRUE, sep = "\t")
# Totalroot <- read.table("totalroot_single_community_test.txt", header = TRUE, sep = "\t")

Elongation <- read.table("~/rhizogenomics/experiments/2017/2017-03-01.wheel_phenotypes/elongation_single_community_test.txt", header = TRUE, sep = "\t")
Area <- read.table("~/rhizogenomics/experiments/2017/2017-03-01.wheel_phenotypes/area_single_community_test.txt", header = TRUE, sep = "\t")
Pi <- read.table("~/rhizogenomics/experiments/2017/2017-03-01.wheel_phenotypes/pi_single_community_test.txt", header = TRUE, sep = "\t")
Totalroot <- read.table("~/rhizogenomics/experiments/2017/2017-03-01.wheel_phenotypes/totalroot_single_community_test.txt", header = TRUE, sep = "\t")

## Added version of plots for Tatiana's comment, undecided which one to use yet
# Plot
dat <- list(Pi = Pi, Elongation = Elongation,
            Area = Area, Totalroot = Totalroot)
svglite::svglite("syncom_effects_clustered.svg", width = 6, height = 6)
tab <- plot_combined_effects(dat = dat)
dev.off()

svglite::svglite("syncom_effects_clustered_standardized.svg", width = 6, height = 6)
tab <- plot_combined_effects(dat = dat,standardize = TRUE)
dev.off()

svglite::svglite("syncom_effects_notclustered.svg", width = 6, height = 6)
tab <- plot_combined_effects(dat = dat,cluster = FALSE)
dev.off()

svglite::svglite("syncom_effects_notclustered_standardized.svg", width = 6, height = 6)
tab <- plot_combined_effects(dat = dat,cluster = FALSE, standardize = TRUE)
dev.off()

dat.sc <- dat
rm(Elongation,Area,Pi,Totalroot,dat,tab)

# Read data
# Elongation <- read.table("elongation_block_test.txt", header = TRUE, sep = "\t")
# Area <- read.table("area_block_test.txt", header = TRUE, sep = "\t")
# Pi <- read.table("pi_block_test.txt", header = TRUE, sep = "\t")
# Totalroot <- read.table("totalroot_block_test.txt", header = TRUE, sep = "\t")

Elongation <- read.table("~/rhizogenomics/experiments/2017/2017-03-01.wheel_phenotypes/elongation_block_test.txt", header = TRUE, sep = "\t")
Area <- read.table("~/rhizogenomics/experiments/2017/2017-03-01.wheel_phenotypes/area_block_test.txt", header = TRUE, sep = "\t")
Pi <- read.table("~/rhizogenomics/experiments/2017/2017-03-01.wheel_phenotypes/pi_block_test.txt", header = TRUE, sep = "\t")
Totalroot <- read.table("~/rhizogenomics/experiments/2017/2017-03-01.wheel_phenotypes/totalroot_block_test.txt", header = TRUE, sep = "\t")

# Plot
dat <- list(Pi = Pi, Elongation = Elongation,
            Area = Area, Totalroot = Totalroot)
svglite::svglite("block_effects_clustered.svg", width = 6, height = 5)
tab <- plot_combined_effects(dat = dat)
dev.off()


tab <- scale(tab, center = TRUE)
tab <- as.data.frame(tab)
tab$Block <- row.names(tab)
tab <- reshape2::melt(tab,id.vars = "Block")
tab <- cbind(tab, do.call(rbind, strsplit(x = as.character(tab$variable), split = "[.]")))
colnames(tab)[4:5] <- c("Phen","cond")
tab$Block <- factor(tab$Block,
                    levels = c("P1","P2","P3","I1","I2","I3","N1","N2","N3"))


p1 <- ggplot2::ggplot(tab, ggplot2::aes(x =  variable, y = Block)) +
  ggplot2::geom_tile(ggplot2::aes(fill = value)) +
  ggplot2::scale_fill_gradient2(low = c("#8e0152","#de77ae"),
                                min = "#f7f7f7",
                                high = c("#7fbc41","#276419")) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                 axis.text = ggplot2::element_text(color = "black"))
p1
ggplot2::ggsave("block_effects_sorted.svg", p1, width = 5, height = 4)

#### Plot all as one
# Prepare syncom data
dat.sc$Pi$Phenotype <- 'Pi'
dat.sc$Elongation$Phenotype <- 'Main'
dat.sc$Area$Phenotype <- 'Area'
dat.sc$Totalroot$Phenotype <- 'Net'

dat.sc <- do.call(rbind,dat.sc)
dat.sc$Type <- 'SynCom'
head(dat.sc)

# Prepare block data
dat$Pi$Phenotype <- 'Pi'
dat$Elongation$Phenotype <- 'Main'
dat$Area$Phenotype <- 'Area'
dat$Totalroot$Phenotype <- 'Net'

dat <- do.call(rbind,dat)
dat$Type <- 'Block'
head(dat)

# Combine
Dat <- rbind(dat.sc, dat)
Dat$Phenotype <- factor(Dat$Phenotype, levels = c('Pi','Main','Area','Net'))
Dat$Type <- factor(Dat$Type, levels = c('SynCom','Block'))
Dat$Significant <- ""
Dat$Significant[ Dat$p.value < 0.05 ] <- 'X'
Dat$Condition <- NULL
Dat$Condition[ Dat$StartP == '-Pi,0.5%Suc'  & Dat$EndP == '100 uM,0%Suc' ] <- 'minusP_100uM'
Dat$Condition[ Dat$StartP == '-Pi,0.5%Suc'  & Dat$EndP == '30 uM,0%Suc' ] <- 'minusP_30uM'
Dat$Condition[ Dat$StartP == '+Pi,0.5%Suc'  & Dat$EndP == '100 uM,0%Suc' ] <- 'plusP_100uM'
Dat$Condition[ Dat$StartP == '+Pi,0.5%Suc'  & Dat$EndP == '30 uM,0%Suc' ] <- 'plusP_30uM'
Dat$Condition <- factor(Dat$Condition, levels = c( 'minusP_100uM','minusP_30uM','plusP_100uM','plusP_30uM'))
Dat$SynCom <- factor(Dat$SynCom, levels = c(rev(c('P1P2','P2P3','P3I1','I1I2','I2I3',
                                              'I3N1','N1N2','N2N3','P1N3','P1P3','P1I1',
                                              'P2I1','P2N3','P3N3')),
                                            rev(c("P1","P2","P3","I1","I2","I3","N1","N2","N3"))))
head(Dat)

Dat$Scaled.Effect <- NA
for(c in levels(Dat$Condition)){
  for(t in levels(Dat$Type)){
    for(p in levels(Dat$Phenotype)){
      # c <- levels(Dat$Condition)[1]
      # t <- levels(Dat$Type)[1]
      # p <- levels(Dat$Phenotype)[1]

      s <- subset(Dat, Condition == c & Type == t & Phenotype == p)
      # # cat(nrow(s),"\n")
      # s <- sd(s$Estimate)

      s <- s$Estimate

      s <- sqrt(sum(s^2 / (length(s) - 1)))


      Dat$Scaled.Effect[ Dat$Condition == c & Dat$Type == t & Dat$Phenotype == p ] <- Dat$Estimate[ Dat$Condition == c & Dat$Type == t & Dat$Phenotype == p ] / s
      # Dat$Scaled.Effect[ Dat$Condition == c & Dat$Type == t & Dat$Phenotype == p ] <- scale(Dat$Estimate[ Dat$Condition == c & Dat$Type == t & Dat$Phenotype == p ], center = FALSE)

    }

  }
}

p1 <- ggplot(Dat, aes(x =  Condition, y = SynCom)) +
  # facet_wrap(~ Type + Phenotype, scales = 'free_y', nrow = 1) +
  facet_grid(~ Type + Phenotype) +
  geom_tile(aes(fill = Scaled.Effect)) +
  geom_text(aes(label = Significant)) +
  scale_fill_gradient2(low = c("#8e0152","#de77ae"),
                                min = "#f7f7f7",
                                high = c("#7fbc41","#276419")) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(color = "black"),
        panel.background = element_blank())
p1
ggsave('heatmap_effects_all.svg',p1, width = 9, height = 6)


