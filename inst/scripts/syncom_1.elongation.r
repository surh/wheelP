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
#library(coefplot)

setwd("~/rhizogenomics/experiments/2017/today11")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")
# source("/home/sur/rhizogenomics/src/trunk/phosphate_code/monoP/functions.r")

Dat <- read.table("~/rhizogenomics/experiments/2015/2015-07-27.wheelP/elongation_full.txt",
                  sep="\t", header = TRUE)
head(Dat)

# Order communities, and rename based on new labels
Dat$Bacteria <- factor(Dat$Bacteria, levels = c("none","G1G2","G2G3","G3N1",
                                                "N1N2","N2N3","N3B1","B1B2",
                                                "B2B3","G1B3","G1N1","G2N1",
                                                "G1G3","G2B3","G3B3"))
levels(Dat$Bacteria) <- c("none","P1P2","P2P3","P3I1",
                          "I1I2","I2I3","I3N1","N1N2",
                          "N2N3","P1N3","P1I1","P2I1",
                          "P1P3","P2N3","P3N3")
# Save data to package
Elongation <- Dat
devtools::use_data(Elongation, pkg = "~/rhizogenomics/src/trunk/phosphate_code/wheelP/",
                   overwrite = TRUE)
rm(Elongation)
###################################
data(Elongation)
Dat <- Elongation

# Plot experimental reproducibility
p1 <- ggplot(Dat,aes(x = Bacteria, y = Elongation, col = Experiment)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x  = element_text(size = 16, angle = 90))
p1
ggsave("elongation_experiment.svg",p1, width = 8, height = 5)

## Fit one at a time
dir <- "elongation_images/"
Res <- test_single_community_phenotype(Dat = Dat, dir = dir, var.name = "Elongation",
                                       bacteria.col = "Bacteria", ref.level = "none",
                                       plot = TRUE, f1.extra = "+Experiment + Plate")
Res
write.table(Res,"elongation_single_community_test.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Plot results of single community analysis
# Prepare data for plotting
singlecoms <- paste(rep(c("P","I","N"),each = 3),rep(1:3,times = 3),sep="")
dat <- Res
dat$Block1 <- substring(dat$SynCom,1,2)
dat$Block2 <- substring(dat$SynCom,3,4)
dat$Block1 <- factor(dat$Block1 , levels = singlecoms)
dat$Block2 <- factor(dat$Block2 , levels = rev(singlecoms))

dat$Significant <- ""
dat$Significant[ dat$p.value < 0.05 ] <- "X"

head(dat)

p1 <- ggplot(dat, aes(x = Block1, y = Block2)) +
  facet_grid(StartP ~ EndP) +
  geom_tile(aes(fill = Estimate), size = 1, col = "black") +
  geom_text(aes(label = Significant), size = 8) +
  scale_fill_gradient2(low = c("#8e0152","#de77ae"),mid = "#f7f7f7",
                       high = c("#7fbc41","#276419"),
                       midpoint = 0,
                       na.value = "#404040",
                       guide = guide_colorbar(title = "SynCom effect\non main root\nelongation (cm)")) +
  scale_x_discrete(limits = singlecoms) +
  scale_y_discrete(limits = rev(singlecoms)) +
  scale_color_manual(values = c("#8c510a","#01665e")) +
  ggtitle("Elongation") +
  theme_classic()
p1
ggsave("elongation_images/heatmap_mix_estimate.svg",p1, width = 8, height = 7)
rm(singlecoms,dat,p1,dir)

## Now we estimate the main effects
Res2 <- test_single_block_phenptype(Phen = Dat,
                                    phenotype = 'Elongation',
                                    confounders = c("Experiment", "Plate"),
                                    use_abun = FALSE)

write.table(Res2,"elongation_block_test.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# compare predictions from main effects
# with measured effects
Pred <- predict_from_main(dat = Res, Res.main = Res2)
head(Pred)

singlecoms <- paste(rep(c("P","I","N"),each = 3),rep(1:3,times = 3),sep="")
Full <- Res
Full$Block1 <- substring(Full$SynCom,1,2)
Full$Block2 <- substring(Full$SynCom,3,4)
Full$Block1 <- factor(Full$Block1 , levels = singlecoms)
Full$Block2 <- factor(Full$Block2 , levels = rev(singlecoms))

Full$Type <- "Measured"
Pred$Type <- "Predicted"
Full <- rbind(Full,Pred)

p1 <- ggplot(Full, aes(x = Block1, y = Block2)) +
  facet_grid(StartP ~ EndP) +
  geom_tile(aes(fill = Estimate), col = "black", size = 1) +
  scale_fill_gradient2(low =  c("#8e0152","#de77ae"),
                       mid = "#f7f7f7",
                       high = c("#7fbc41","#276419"),
                       midpoint = 0,
                       na.value = "#404040",
                       guide = guide_colorbar(title = "SynCom effect\non main root\nelongation (cm)")) +
  scale_x_discrete(limits = singlecoms) +
  scale_y_discrete(limits = rev(singlecoms)) +
  ggtitle(label = "Elongation") +
  theme_classic()
p1
ggsave("elongation_images/heatmap_mix_full.svg",p1, width = 8, height = 7)

Res$Measured <- Res$Estimate
Res$Predicted <- Pred$Estimate
Res$SE.pred <- Pred$SE
p1 <- ggplot(Res,aes(x = Measured, y = Predicted)) +
  facet_grid(StartP ~ EndP) +
  geom_point(size = 3) +
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, xpos = -3.5, ypos = -0.5) +
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
ggsave("elongation_images/estimate_vs_pred.svg",p1, width = 7, height = 7)
dir.create("figuredata/")
fig4.elongation <- p1$data
save(fig4.elongation, file = "figuredata/fig4.elongation.rda")
