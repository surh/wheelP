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
library(edgeR)

# setwd("~/rhizogenomics/experiments/2017/today5//")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")

data(wheelP.rna)
Dat <- wheelP.rna
rm(wheelP.rna)

# Read core
pi.core <- read.table("~/rhizogenomics/data/phosphate/core_pi_up.txt",
                      stringsAsFactors = FALSE)$V1
meja.core <- read.table("~/rhizogenomics/data/phosphate/MeJA_core.txt",
                        stringsAsFactors = FALSE)$V1
bth.core <- read.table("~/rhizogenomics/data/phosphate/BTH_core.txt",
                       stringsAsFactors = FALSE)$V1
bth.meja.core <- read.table("~/rhizogenomics/data/phosphate/BTH_MeJA_core.txt",
                            stringsAsFactors = FALSE)$V1
flg22.acute.core <- read.table("~/rhizogenomics/data/phosphate/fl22_acute_core.txt",
                               stringsAsFactors = FALSE)$V1
syncomP.core <- read.table("~/rhizogenomics/data/phosphate/PBI_RNAseq/2017-03-05.syncom.responsive.genes.txt",
                            stringsAsFactors = FALSE, header = TRUE, sep = "\t")
syncomP.core <- syncomP.core$Gene[ syncomP.core$Fold > 0 ]

syncomP.core.dn <- read.table("~/rhizogenomics/data/phosphate/PBI_RNAseq/2017-03-05.syncom.responsive.genes.txt",
                           stringsAsFactors = FALSE, header = TRUE, sep = "\t")
syncomP.core.dn <- syncomP.core.dn$Gene[ syncomP.core.dn$Fold < 0 ]

bactlowp.core <- read.table("~/rhizogenomics/data/phosphate/PBI_RNAseq/2017-03-05.syncomlowP.responsive.genes.txt",
                           stringsAsFactors = FALSE, header = TRUE, sep = "\t")
bactlowp.core <- bactlowp.core$Gene[ bactlowp.core$Fold > 0 ]

# Gene lengths for RPKM
gene.lengths <- read.table("~/rhizogenomics/data/phosphate/gene_lengths.txt",
                           row.names = 1, header = TRUE)

# Calculate RPKM
Dat.norm <- create_dataset(Tab = rpkm(x = Dat$Tab,
                                      gene.length = gene.lengths[ taxa(Dat), ]),
                           Map = Dat$Map,
                           Tax = Dat$Tax)
# Calculate z-score
Dat.norm <- create_dataset(t(scale(t(Dat$Tab))),Dat$Map)

# Plot cores
p1 <- plot_set(Dat = Dat.norm,gene.set = pi.core,
               groups = c("Bacteria","Phosphate"))
p1
ggsave("corePi.wheel.png",p1,width = 12, height = 4)
ggsave("corePi.wheel.svg",p1,width = 12, height = 4)

p1 <- plot_set(Dat = Dat.norm,gene.set = bth.core,
               groups = c("Bacteria","Phosphate"))
p1
ggsave("coreBTH.wheel.png",p1,width = 12, height = 4)
ggsave("coreBTH.wheel.svg",p1,width = 12, height = 4)

p1 <- plot_set(Dat = Dat.norm,gene.set = meja.core,
               groups = c("Bacteria","Phosphate"))
p1
ggsave("coreMeJA.wheel.png",p1,width = 12, height = 4)
ggsave("coreMeJA.wheel.svg",p1,width = 12, height = 4)

p1 <- plot_set(Dat = Dat.norm,gene.set = bth.meja.core,
               groups = c("Bacteria","Phosphate"))
p1
ggsave("coreBTHMeJA.wheel.png",p1,width = 12, height = 4)
ggsave("coreBTHMeJA.wheel.svg",p1,width = 12, height = 4)

p1 <- plot_set(Dat = Dat.norm,gene.set = flg22.acute.core,
               groups = c("Bacteria","Phosphate"))
p1
ggsave("coreflg22acute.wheel.png",p1,width = 12, height = 4)
ggsave("coreflg22acute.wheel.svg",p1,width = 12, height = 4)

p1 <- plot_set(Dat = Dat.norm,gene.set = syncomP.core,
               groups = c("Bacteria","Phosphate"))
p1
ggsave("coresyncomp.wheel.png",p1,width = 12, height = 4)
ggsave("coresyncomp.wheel.svg",p1,width = 12, height = 4)

p1 <- plot_set(Dat = Dat.norm,gene.set = syncomP.core.dn,
               groups = c("Bacteria","Phosphate"))
p1
ggsave("coresyncompdn.wheel.png",p1,width = 12, height = 4)
ggsave("coresyncompdn.wheel.svg",p1,width = 12, height = 4)

p1 <- plot_set(Dat = Dat.norm,gene.set = bactlowp.core,
               groups = c("Bacteria","Phosphate"))
p1
ggsave("coresyncomlowp.wheel.png",p1,width = 12, height = 4)
ggsave("coresyncomlowp.wheel.svg",p1,width = 12, height = 4)


######### ADDITIVITY #####
# Mostly correlates with other phenotypes
# Calculate index for core
Dat.sub <- remove_taxons(Dat.norm,setdiff(taxa(Dat.norm),pi.core))
Dat.sub$Map$core.index <- collapse_matrix(Dat.sub$Tab,
                                          groups = rep("core.index",
                                                       nrow(Dat.sub$Tab)),
                                          dim = 1,FUN = mean)
plotgg_var(Dat.sub, var.name = "core.index",x = "Phosphate",col = "Bacteria")

# Prepare data
dat <- Dat.sub$Map
dat$Well <- dat$Genotype <- dat$Concentration <- dat$A260.230 <- dat$A260.280 <- NULL
dat$ID <- dat$Sample <- NULL
dat$StartP <- NULL
dat$EndP <- NULL
dat$StartP[ dat$Phosphate %in% c('minusP_30uM','minusP_100uM') ] <- "-Pi,0.5%Suc"
dat$StartP[ dat$Phosphate %in% c('plusP_30uM', 'plusP_100uM') ] <- "+Pi,0.5%Suc"
dat$EndP[ dat$Phosphate %in% c('plusP_100uM','minusP_100uM') ] <- "100 uM,0%Suc"
dat$EndP[ dat$Phosphate %in% c('plusP_30uM','minusP_30uM') ] <- "30 uM,0%Suc"
dat$Extraction <- factor(dat$Extraction)
dat$Plate <- factor(dat$Plate)
dat$Experiment <- factor(dat$Experiment)
dat$StartP <- factor(dat$StartP)
dat$EndP <- factor(dat$EndP)
head(dat)

# Test community
# setwd("~/rhizogenomics/experiments/2017/today6/")
Res.sc <- test_single_community_phenotype(Dat = dat,
                                          var.name = 'core.index',
                                          dir = './pi_core_images/',
                                          bacteria.col = 'Bacteria',
                                          ref.level = 'No Bacteria',
                                          plot = TRUE,f1.extra = "+ Experiment + Plate + Extraction")
summary(qvalue::qvalue(Res.sc$p.value))

# Test block
Res.block <- test_single_block_phenptype(Phen = dat,phenotype = 'core.index',
                                         confounders = c('Experiment', 'Plate', 'Extraction'),
                                         use_abun = FALSE)
summary(qvalue::qvalue(Res.block$p.value))

# Predict from main
Pred <- predict_from_main(dat = Res.sc, Res.main = Res.block)

Res.sc$Measured <- Res.sc$Estimate
Res.sc$Predicted <- Pred$Estimate
Res.sc$SE.pred <- Pred$SE
head(Res.sc)

p1 <- ggplot(Res.sc,aes(x = Measured, y = Predicted)) +
  facet_grid(StartP ~ EndP) +
  geom_point(size = 3) +
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE,xpos = 0, ypos = 2.5, size = 3) +
  geom_errorbarh(aes(xmin = Measured - SE, xmax = Measured + SE)) +
  geom_errorbar(aes(ymin = Predicted - SE.pred, ymax = Predicted + SE.pred)) +
  geom_smooth(method = 'lm') +
  #geom_abline(intercept = 0, slope =1) +
  xlab(label = "Measured (cm)") +
  theme_classic() +
  theme(axis.title = element_text(face = "bold",
                                  color = "black",
                                  size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(face = "bold",
                                  size = 22))
p1
ggsave('corepi.additivity.svg',p1,width = 7, height = 7)

#### THE FOLLOWING SECTION IS INCORRECT!!
#### IT CALCULATES GOODNESS OF FIT OF ADDITIVE ONLY MODEL
#### BUT IT DOES NOT COMPARE IT TO MEASURED COMMUNITY EFFECTS
#### LIKE FOR THE OTHER PHENOTYPES

### WILL FIX, THIS IS LOW PRIORITY SINCE IT MOSTLY
### CORRELATES WITH OTHER PHENOTYPES

#
# f1 <- formula(core.index ~ P1 + P2 + P3 +
#                 I1 + I2 + I3 + N1 + N2 + N3 +
#                 Experiment)
#
# # Fit per conditions
# Res <- NULL
# dat2 <- subset(dat,Phosphate == "minusP_100uM")
# m1 <- lm(f1,
#          data = dat2)
# summary(m1)
# dat2$Pred <- fitted(m1)
# res <- cbind(aggregate(core.index ~ Bacteria,dat2,FUN = mean),
#              aggregate(Pred ~ Bacteria,dat2,FUN = mean))
# res$PrePi <- "minusP"
# res$PosPi <- "100uM"
# Res <- rbind(Res,res)
#
# dat2 <- subset(dat,Phosphate == "minusP_30uM")
# m1 <- lm(f1,
#          data = dat2)
# summary(m1)
# dat2$Pred <- fitted(m1)
# res <- cbind(aggregate(core.index ~ Bacteria,dat2,FUN = mean),
#              aggregate(Pred ~ Bacteria,dat2,FUN = mean))
# res$PrePi <- "minusP"
# res$PosPi <- "30uM"
# Res <- rbind(Res,res)
#
# dat2 <- subset(dat,Phosphate == "plusP_30uM")
# m1 <- lm(f1,
#          data = dat2)
# summary(m1)
# dat2$Pred <- fitted(m1)
# res <- cbind(aggregate(core.index ~ Bacteria,dat2,FUN = mean),
#              aggregate(Pred ~ Bacteria,dat2,FUN = mean))
# res$PrePi <- "plusP"
# res$PosPi <- "30uM"
# Res <- rbind(Res,res)
#
# dat2 <- subset(dat,Phosphate == "plusP_100uM")
# m1 <- lm(f1,
#          data = dat2)
# summary(m1)
# dat2$Pred <- fitted(m1)
# res <- cbind(aggregate(core.index ~ Bacteria,dat2,FUN = mean),
#              aggregate(Pred ~ Bacteria,dat2,FUN = mean))
# res$PrePi <- "plusP"
# res$PosPi <- "100uM"
# Res <- rbind(Res,res)
#
# Res <- Res[ Res$Bacteria != "No Bacteria", ]
# p1 <- ggplot(Res,aes(x = core.index, y = Pred)) +
#   facet_grid(PrePi ~ PosPi) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   theme_bw()
# p1
# ggsave("core.pi.index.png",p1,width = 4,height = 4)
# ggsave("core.pi.index.svg",p1,width = 4,height = 4)
#
