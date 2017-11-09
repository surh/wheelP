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

# setwd("~/rhizogenomics/experiments/2017/today4/")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")
# source("/home/sur/rhizogenomics/src/trunk/phosphate_code/monoP/functions.r")

Dat <- read.table("~/rhizogenomics/data/synthetic/wheel_phosphate/2017-03-01.master.txt",
                  sep="\t", header = TRUE)

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
Pi <- Dat
devtools::use_data(Pi, pkg = "~/rhizogenomics/src/trunk/phosphate_code/wheelP/",
                   overwrite = TRUE)
rm(Pi)

##########
data(Pi)
Dat <- Pi
rm(Pi)

# Plot experimental reproducibility
p1 <- ggplot(Dat,aes(x = Bacteria, y = Pi_content, col = Experiment)) +
  geom_boxplot() +
  scale_y_log10(breaks = c(1,5,10,15,20)) +
  theme_classic() +
  theme(axis.text.x  = element_text(size = 16, angle = 90))
p1
ggsave("pi_experiment.svg",p1, width = 5, height = 5)

## Fit one at a time
dir <- "pi_images/"
Dat$LogPi <- log(Dat$Pi_content)
Res <- test_single_community_phenotype(Dat = Dat, dir = dir, var.name = "LogPi",
                                       bacteria.col = "Bacteria", ref.level = "none",
                                       plot = TRUE, f1.extra = "+ Experiment")
Res
write.table(Res,"pi_single_community_test.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

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
                       guide = guide_colorbar(title = "Fold change in\nPi content between\nSynCom and \nNo Bacteria")) +
  scale_x_discrete(limits = singlecoms) +
  scale_y_discrete(limits = rev(singlecoms)) +
  scale_color_manual(values = c("#8c510a","#01665e")) +
  ggtitle(label = "Pi content") +
  theme_classic()
p1
ggsave("pi_images/heatmap_mix_estimate.svg",p1, width = 8, height = 7)
rm(singlecoms,dat,p1,dir)

## Now we estimate the main effects
Res2 <- test_single_block_phenptype(Phen = Dat,
                                    phenotype = 'LogPi',
                                    confounders = c("Experiment"),
                                    use_abun = FALSE)

# singlecoms <- paste(rep(c("P","I","N"),each = 3),rep(1:3,times = 3),sep="")
#
# Res2 <- NULL
# for(i in 1:2){
#   for(j in 1:2){
#     m1 <- block_effects(Dat = Dat, cond1 = i, cond2 = j,
#                         var.name = "LogPi",
#                         keep.vars = c("Experiment"))
#     m1.sum <- summary(m1)
#
#     res <- data.frame(SynCom = singlecoms, StartP = levels(Dat$StartP)[i],
#                       EndP = levels(Dat$EndP)[j],
#                       Estimate = m1.sum$coefficients[ singlecoms, 1 ],
#                       SE = m1.sum$coefficients[ singlecoms, 2 ],
#                       t.value = m1.sum$coefficients[ singlecoms, 3 ],
#                       p.value = m1.sum$coefficients[ singlecoms, 4 ])
#     Res2 <- rbind(Res2,res)
#     # p1 <- coefplot(m1, intercept = FALSE,innerCI = 1, outerCI = 2,
#     #                coefficients = singlecoms,
#     #                lwdOuter = 0.5, lwdInner = 2.5, pointSize = 4,
#     #                color = "black", zeroColor = "red", zeroType = 1)
#     # p1 <- p1 + theme(axis.text.y = element_text(color = "black"),
#     #                  panel.border = element_rect(fill = NA, color = "black", size = 2),
#     #                  panel.background = element_rect(fill = "white")) +
#     #   ggtitle(paste(levels(Dat$StartP)[i],"=>",levels(Dat$EndP)[j]))
#     # p1
#     #
#     # filename <- paste(dir,"/coefplot_",i,j,".png",sep = "")
#     # ggsave(filename = filename, p1 , width = 3.5, height = 5)
#   }
# }
write.table(Res2,"pi_block_test.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# compare predictions from main effects
# with measured effects
Pred <- predict_from_main(dat = Res, Res.main = Res2)
head(Pred)

# # Format observed effects
# dat <- Res
# dat$Block1 <- substring(dat$SynCom,1,2)
# dat$Block2 <- substring(dat$SynCom,3,4)
# dat$Block1 <- factor(dat$Block1 , levels = singlecoms)
# dat$Block2 <- factor(dat$Block2 , levels = rev(singlecoms))
#
# # Predict from main effects
# Pred <- NULL
# for(i in 1:nrow(Res)){
#   # i <- 1
#   index1 <- Res2$StartP == dat$StartP[i] &
#     Res2$EndP == dat$EndP[i] &
#     Res2$SynCom == dat$Block1[i]
#   index2 <- Res2$StartP == dat$StartP[i] &
#     Res2$EndP == dat$EndP[i] &
#     Res2$SynCom == dat$Block2[i]
#   additiveguess <- Res2$Estimate[index1] + Res2$Estimate[index2]
#
#   res <- data.frame(SynCom = paste(Res2$SynCom[index1],Res2$SynCom[index2],sep = ""),
#                     StartP = dat$StartP[i], EndP = dat$EndP[i],
#                     Estimate = additiveguess, SE = NA, t.value = NA,
#                     p.value = NA, Block1 = Res2$SynCom[index2],
#                     Block2 = Res2$SynCom[index1])
#   Pred <- rbind(Pred,res)
# }

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
                       guide = guide_colorbar(title = "SynCom effect\non Pi content\nelongation (cm)")) +
  scale_x_discrete(limits = singlecoms) +
  scale_y_discrete(limits = rev(singlecoms)) +
  ggtitle(label = "Pi content") +
  theme_classic()
p1
ggsave("pi_images/heatmap_mix_full.svg",p1, width = 8, height = 7)

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
  xlab(label = "Measured (log fold-change)") +
  theme_classic() +
  theme(axis.title = element_text(face = "bold",
                                  color = "black",
                                  size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(face = "bold",
                                  size = 22))
p1
ggsave("pi_images/estimate_vs_pred.svg",p1, width = 7, height = 7)
