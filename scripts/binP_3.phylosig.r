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
library(ggtree)
library(phytools)
# setwd("~/rhizogenomics/experiments/2017/2017-02-28.binP_phylosig//")
data(Features,package = "PEBG")

binP <- read.table("~/rhizogenomics/data/phosphate/exudate_isolate_screening/features/strain_test_results.txt",
                   sep = "\t", header = TRUE)
head(binP)
row.names(binP) <- as.character(binP$Strain)
binP$Strain.1 <- NULL
colnames(binP)[6:9] <- paste(colnames(binP[2:5]),".qval",sep = "")
head(binP)


tested <- matrix(binP$Strain,ncol = 1)
tested[,1] <- gsub(pattern = "^M",replacement = "",x = tested[,1])
binP$Strain <- tested[,1]
row.names(binP) <- tested[,1]

devtools::use_data(binP,
                   pkg = "~/rhizogenomics/src/trunk/phosphate_code/wheelP/",
                   overwrite = TRUE)
write.table(tested, "binP.tested.txt", sep = "\t",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

tree <- read.tree("binP_tree/filter/ref_aligned_pfiltered.tre")
Dat <- binP[ tree$tip.label, ]
Dat <- Dat[ ,2:5 ]
head(Dat)

ladder <- TRUE
pal <- colorRampPalette(colors = c("#8e0152","#de77ae",
                                   "#f7f7f7","#7fbc41","#276419"))
p1 <- ggtree::ggtree(tree, ladderize = ladder)
p1 <- ggtree::gheatmap(p1,as.data.frame(Dat),
                       colnames = FALSE, color = NA,colnames_position = "top") +
  scale_fill_gradient2(low =c("#8e0152","#de77ae"),
                       mid = "#f7f7f7", high = c("#7fbc41","#276419"),
                       midpoint = 0)
#+coord_polar(theta = "y")
df <- ggtree::get_heatmap_column_position(p1, by = "top")
p1 <- p1 + geom_text(data = df, aes(x, y, label = label), angle = 90, nudge_y = 17)
p1
ggsave("growth_curve_features_heatmap.svg",p1, width = 5, height = 7)


RES <- NULL
for(i in 1:ncol(Dat)){
  # i <- 1
  res <- phytools::phylosig(tree = tree,x = Dat[,i],
                            test = TRUE, method = "lambda")
  res <- data.frame(Phenotype = colnames(Dat)[i],
                    lambda = res$lambda, p.value = res$P)
  RES <- rbind(RES,res)

}
head(RES)

p1 <- ggplot(RES,aes(x = Phenotype, y = lambda)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90),
        panel.background = element_blank())
p1
ggsave("pagel.lambda.svg", p1, width = 6, height = 4)

p1 <- ggplot(RES,aes(x = Phenotype, y = p.value)) +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90),
        panel.background = element_blank())
p1
ggsave("pagel.lambda.pval.svg", p1, width = 6, height = 4)

max(RES$p.value)
min(RES$lambda)

#### correlationsw with GC

Dat <- binP[,c(2:5,10)]
head(binP)
head(Dat)
Dat <- cbind(Dat, Features[ row.names(Dat), 2:5])
head(Dat)

Dat <- Dat[ apply(!is.na(Dat),1,all), ]
head(Dat)

group1 <- 1:5
group2 <- 6:9

RES <- NULL
for(g1 in group1){
  for(g2 in group2){
    # g1 <- 1
    # g2 <- 6
    #
    res <- data.frame(AUC = colnames(Dat)[g2], binP = colnames(Dat)[g1],
                      x = Dat[,g2], y = Dat[,g1])

    RES <- rbind(RES,res)

  }

}


head(RES)

p1 <- ggplot(RES, aes(x,y)) +
  facet_grid(binP ~ AUC) +
  geom_point(aes(color = AUC), size = 0.5) +
  geom_smooth(method = "loess") +
  AMOR::theme_blackbox +
  theme(strip.text.y = element_text(angle = 0))
p1
ggsave("auc_vs_binP.svg", p1, width = 8, height = 6)
