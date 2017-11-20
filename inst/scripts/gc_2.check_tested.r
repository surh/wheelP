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
data("Strain.auc")
data("All.filtered")

# Scale
Strain.auc <- log(Strain.auc)
Strain.auc <- t(scale(t(Strain.auc)))
attr(Strain.auc,which = "scaled:center") <- NULL
attr(Strain.auc,which = "scaled:scale") <- NULL

# Cluster and get groups
area.clus <- hclust(dist(Strain.auc, method = "euclidean"))
k <- 10
groups <- cutree(area.clus,k=k)
group.order <- unique(groups[area.clus$labels[ area.clus$order]])

# Reformat
Strain.auc <- data.frame(Strain.auc)
Strain.auc$Strain <- row.names(Strain.auc)
Strain.auc <- reshape2::melt(data = Strain.auc,
                   id.vars = "Strain",
                   value.name = "Zscore",
                   variable.name = "condition")

# Reorder for plotting
Strain.auc$Group <- factor(groups[ as.character(Strain.auc$Strain) ],
                           levels = rev(group.order))
levels(Strain.auc$Group) <- 1:10
Strain.auc$Strain <- factor(Strain.auc$Strain,
                            levels = rev(area.clus$labels[ area.clus$order ]))
Strain.auc$condition <- factor(Strain.auc$condition,
                               levels = c("minus2plusP",
                                          "plus2minusP",
                                          "minusP","plusP"))
head(Strain.auc)

groups <- apply(reshape2::acast(data = Strain.auc,formula = Strain ~ condition,value.var = "Group"),1,unique)
r.squared <- strain_rsqrd(All.filtered)
row.names(r.squared) <- r.squared$Strain

Dat <- data.frame(Strain = names(groups),Group = groups,
                  r.squared = r.squared[names(groups),"r.squared"] )
head(Dat)

# write.table(Dat, "~/rhizogenomics/experiments/2017/today10/clusters.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(Dat, "clusters.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#b7-112

# tested <- read.table("~/rhizogenomics/data/phosphate/exudate_isolate_screening/features/strain_test_results.txt",
#                      sep = "\t", header = TRUE)
data(binP)
tested <- binP
rm(binP)

head(tested)
tested <- as.character(tested$Strain)
tested <- gsub(pattern = "^M",replacement = "",x = tested)
tested

Dat$Tested <- Dat$Strain %in% tested
table(Dat$Tested)
Dat$Strain <- factor(Dat$Strain, levels = as.character(Dat$Strain[ order(Dat$r.squared, decreasing = TRUE) ]))

head(Dat)


p1 <- ggplot(Dat, aes(x = Strain, y = Tested)) +
  facet_wrap(~ Group, scales = "free", ncol = 2) +
  geom_bar(stat = "identity", col = "black", fill = "black") +
  theme(axis.text.x = element_blank())
p1
ggsave("check_tested.svg",p1,width = 5,height = 5)
