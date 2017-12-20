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

library(ggplot2)
library(wheelP)

# setwd("~/rhizogenomics/experiments/2017/today11")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")
# source("/home/sur/rhizogenomics/src/trunk/phosphate_code/monoP/functions.r")

# Read data, same object as Pi
data(Pi)
Dat <- Pi
rm(Pi)

# Plot experimental reproducibility
p1 <- ggplot(Dat,aes(x = Bacteria, y = total.root, col = Experiment)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x  = element_text(size = 16, angle = 90))
p1
ggsave("totalroot_experiment.svg",p1,width = 5, height = 5)

## Fit one at a time
dir <- "totalroot_images/"
Res <- test_single_community_phenotype(Dat = Dat, dir = dir, var.name = "total.root",
                                       bacteria.col = "Bacteria", ref.level = "none",
                                       plot = TRUE, f1.extra = "+ Experiment")
Res
write.table(Res,"totalroot_single_community_test.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)

# Plot results of single community analysis
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
                       guide = guide_colorbar(title = "Difference in\nroot network length\n between SynCom and \nNo Bacteria (cm)")) +
  scale_x_discrete(limits = singlecoms) +
  scale_y_discrete(limits = rev(singlecoms)) +
  scale_color_manual(values = c("#8c510a","#01665e")) +
  ggtitle(label = "Total Root Network") +
  theme_classic()
p1
ggsave("totalroot_images/heatmap_mix_estimate.svg",p1, width = 8, height = 7)
rm(singlecoms,dat,p1,dir)

## Now we estimate the main effects
Res2 <- test_single_block_phenptype(Phen = Dat,
                                    phenotype = 'total.root',
                                    confounders = c("Experiment"),
                                    use_abun = FALSE)

write.table(Res2,"totalroot_block_test.txt",
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
                       guide = guide_colorbar(title = "SynCom effect\ntotal root network\nbetween SynCom\nand no bacteria")) +
  scale_x_discrete(limits = singlecoms) +
  scale_y_discrete(limits = rev(singlecoms)) +
  ggtitle(label = "Total Root Network") +
  theme_classic()
p1
ggsave("totalroot_images/heatmap_mix_full.svg",p1, width = 8, height = 7)

Res$Measured <- Res$Estimate
Res$Predicted <- Pred$Estimate
Res$SE.pred <- Pred$SE
p1 <- ggplot(Res,aes(x = Measured, y = Predicted)) +
  facet_grid(StartP ~ EndP) +
  geom_point(size = 3) +
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE, xpos = -1.7, ypos = 3) +
  geom_errorbarh(aes(xmin = Measured - SE, xmax = Measured + SE)) +
  geom_errorbar(aes(ymin = Predicted - SE.pred, ymax = Predicted + SE.pred)) +
  geom_smooth(method = 'lm') +
  # geom_abline(intercept = 0, slope =1) +
  xlab(label = "Measured (cm)") +
  theme_classic() +
  theme(axis.title = element_text(face = "bold",
                                  color = "black",
                                  size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(face = "bold",
                                  size = 22))
p1
ggsave("totalroot_images/estimate_vs_pred.svg",p1, width = 7, height = 7)
dir.create("figuredata/")
fig4.totalroot <- p1$data
save(fig4.totalroot, file = "figuredata/fig4.totalroot.rda")
