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

############# Adding data from flatfiles ######################
data(Strain.auc)
head(Strain.auc)
write.table(row.names(Strain.auc),file = "strains.auc.id.txt",
            row.names = FALSE,col.names = FALSE, quote = FALSE)

Features <- cbind(Strain = row.names(Strain.auc), as.data.frame(log(Strain.auc)))
head(Features)

# Reading data from flatfile. Data is now in package as "Features"
feat <- read.table("~/rhizogenomics/data/phosphate/exudate_isolate_screening/features/PhBacFeature.txt",
                   sep = "\t", header = TRUE, comment.char = "", row.names = 1)
head(feat)

# Homogenizing features with the rest of the info about growth cruves
setdiff(row.names(feat),Features$Strain)
setdiff(Features$Strain,row.names(feat))

feat <- feat[ , -grep("_AUC$", colnames(feat))]
Features <- cbind(Features, feat[row.names(Features),])
write.table(Features, "2017-01-27.features.txt")
devtools::use_data(Features, pkg = "~/rhizogenomics/src/trunk/phosphate_code/wheelP/",
                   overwrite = FALSE)

# Reading tree_file from flatfile
gc.tree <- read.tree("~/rhizogenomics/experiments/2017/2017-02-27.growthcurves_phylosig/growth_curves_tree/filter/ref_aligned_pfiltered.tre")
devtools::use_data(gc.tree, pkg = "~/rhizogenomics/src/trunk/phosphate_code/wheelP/",
                   overwrite = FALSE)
rm(feat,Strain.auc,gc.tree,Features)
####################################################################
# From here it should run without access to the flatfiles
data(Strain.auc)
data(Features)
data(gc.tree)
tree <- gc.tree
rm(gc.tree)

Dat <- Features[ tree$tip.label, -1]
Dat <- as.data.frame(scale(Dat))
Dat[ Dat > 3] <- 3
Dat[ Dat < -3 ] <- -3
head(Dat)
dat.clus <- hclust(dist(t(Dat)),method = "complete")
Dat <- Dat[ ,dat.clus$order ]

# Plot
ladder <- TRUE
pal <- colorRampPalette(colors = c("#8e0152","#de77ae",
                                   "#f7f7f7","#7fbc41","#276419"))
p1 <- ggtree::ggtree(tree, ladderize = ladder)
# p1 <- ggtree::groupOTU(object = p1,cls) +
#   #aes(color = group) +
#   ggtree::geom_tiplab(aes(label = Isolate)) +
#   scale_color_manual(values = c(phyla_colors),
#                      labels = color_labels)
# p1 <- p1 + ggtree::geom_tippoint(aes_string(color = var))
p1 <- ggtree::gheatmap(p1,as.data.frame(Dat),
                       colnames = FALSE, color = NA,colnames_position = "top") +
  scale_fill_gradientn(colours = pal(10),trans = "identity",na.value = "white")
#+coord_polar(theta = "y")
df <- ggtree::get_heatmap_column_position(p1, by = "top")
p1 <- p1 + geom_text(data = df, aes(x, y, label = label), angle = 90, nudge_y = 22)
p1
ggsave("growth_curve_features_heatmap.svg",p1, width = 10, height = 10)

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
