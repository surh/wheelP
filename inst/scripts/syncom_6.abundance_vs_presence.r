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
devtools::document()

setwd("~/rhizogenomics/experiments/2017/today5/")
outdir <- "~/rhizogenomics/experiments/2017/today5/"

data(Elongation)
data(wheelP.mapsplit)

Dat <- obtain_block_abundances(Dat = wheelP.mapsplit,
                               varnames = c("Bacteria","Experiment","Pre.Pi", "Pos.Pi"),
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
Dat$Experiment <- Dat$EXP
Dat$EXP <- NULL

# Rename condition names
Dat$StartP <- Dat$Pre.Pi
Dat$EndP <- Dat$Pos.Pi
Dat$Pre.Pi <- NULL
Dat$Pos.Pi <- NULL

# Collapse elongation to match sequencing data
Phen <- Elongation
Phen <- aggregate(Elongation ~ Bacteria + Experiment + StartP + EndP, data = Phen, FUN = mean)
head(Phen)

# Test phenotype on collapsed data
Res.sc <- test_single_community_phenotype(Dat = Phen,dir = "Elongation/",
                                          var.name = "Elongation",
                                          bacteria.col = "Bacteria", ref.level = 'none',
                                          plot = FALSE, f1.extra = "+ Experiment")
summary(qvalue::qvalue(Res.sc$p.value))

# Now test block effects
Res.abun <- test_single_block_phenptype(Phen = Phen,
                                        abun = Dat,
                                        phenotype = 'Elongation',
                                        variables = c("Bacteria", "Experiment", "StartP", "EndP"),
                                        confounders = c("Experiment"),
                                        use_abun = TRUE)
summary(qvalue::qvalue(Res.abun$p.value))

Res.block <- test_single_block_phenptype(Phen = Phen,
                                         phenotype = 'Elongation',
                                         variables = c("Bacteria", "Experiment", "StartP", "EndP"),
                                         confounders = c("Experiment"),
                                         use_abun = FALSE)
summary(qvalue::qvalue(Res.block$p.value))

# Compare abundance and not abundance based
head(Res.block)
head(Res.abun)

dat <- Res.block
dat$t.value <- dat$p.value <- NULL
dat$Estimate.block <- dat$Estimate
dat$SE.block <- dat$SE
dat$Estimate <- dat$SE <- NULL
dat$Estimate.abun <- Res.abun$Estimate
dat$SE.abun <- Res.abun$SE
head(dat)

p1 <- ggplot(dat,aes(x = Estimate.block, y = Estimate.abun)) +
  facet_grid(StartP ~ EndP) +
  geom_point() +
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, xpos = -2, ypos = 1.5) +
  geom_errorbar(aes(ymin = Estimate.abun - SE.abun, ymax = Estimate.abun + SE.abun)) +
  geom_errorbarh(aes(xmin = Estimate.block - SE.block, xmax = Estimate.block + SE.block)) +
  geom_smooth(method = "lm") +
  theme_classic() +
  theme(axis.title = element_text(face = "bold",
                                  color = "black",
                                  size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(face = "bold",
                                  size = 22))
p1
ggsave("block_effect_comparison.svg", p1, width = 7, height = 7)


# Predict from main effects
Pred.block <- predict_from_main(dat = Res.sc, Res.main = Res.block)

dat <- Res.sc
dat$Measured <- dat$Estimate
dat$Predicted <- Pred.block$Estimate
dat$SE.pred <- Pred.block$SE
head(dat)
p1 <- ggplot(dat,aes(x = Measured, y = Predicted)) +
  facet_grid(StartP ~ EndP) +
  geom_point(size = 3) +
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, xpos = -3, ypos = 1.5) +
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

Pred.abun <- predict_from_main(dat = Res.sc, Res.main = Res.abun)

dat <- Res.sc
dat$Measured <- dat$Estimate
dat$Predicted <- Pred.abun$Estimate
dat$SE.pred <- Pred.abun$SE
head(dat)
p1 <- ggplot(dat,aes(x = Measured, y = Predicted)) +
  facet_grid(StartP ~ EndP) +
  geom_point(size = 3) +
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, xpos = -3, ypos = 1.5) +
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

# Plot together
dat <- cbind(rbind(Res.sc,Res.sc),
             rbind(data.frame(Predicted = Pred.abun$Estimate,
                              SE.pred = Pred.abun$SE, Type = "Abundance"),
                   data.frame(Predicted = Pred.block$Estimate,
                              SE.pred = Pred.block$SE, Type = "Block")))
head(dat)
p1 <- ggplot(dat,aes(x = Estimate, y = Predicted, color = Type)) +
  facet_grid(StartP ~ EndP) +
  geom_point(size = 3) +
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, xpos = -3.5) +
  geom_errorbarh(aes(xmin = Estimate - SE, xmax = Estimate + SE)) +
  geom_errorbar(aes(ymin = Predicted - SE.pred, ymax = Predicted + SE.pred)) +
  geom_smooth(method = 'lm') +
  xlab(label = "Measured (cm)") +
  theme_classic() +
  theme(axis.title = element_text(face = "bold",
                                  color = "black",
                                  size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(face = "bold",
                                  size = 22))
p1
ggsave("abundance_vs_block_predictions.svg", p1, width = 7, height = 7)
