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
library(gplots)
library(ggplot2)
# setwd("~/rhizogenomics/experiments/2017/today4/")

data("All.filtered")

#################
# Test AUC
RES <- test_auc(Dat = All.filtered)
head(RES, 10)

# Find strains that are always different between exudates and no exudates,
# and strains that are different between exudates
Tab <- reshape2::acast(formula = Strain ~ comparison,data = RES,value.var = "significant")
exudate_strains <- row.names(Tab[ Tab[,"minusP-minus2plusP"] == "yes" &
                                    Tab[,"plus2minusP-minusP"] == "yes" &
                                    Tab[,"plusP-minus2plusP"] == "yes" &
                                    Tab[,"plusP-plus2minusP"] == "yes",])
exudate_difference_stains <- row.names(Tab[ Tab[,"plus2minusP-minus2plusP"] == "yes" ,])


# Plot strains that are different in at leaast one pair of conditions
Tab <- reshape2::acast(formula = comparison ~ Strain,data = RES,value.var = "difference")
Tab <- Tab[ ,colnames(Tab) %in% c(exudate_difference_stains,exudate_strains)]
dim(Tab)
Tab[1:5,1:5]

pal <- colorRampPalette(colors = c("#2c7bb6","#ffffbf","#d7191c"))
svglite::svglite("differences_heatmap.svg",width = 10, height = 5)
gplots::heatmap.2(Tab,col = pal(10),trace = "none",margins = c(4,17))
dev.off()

rm(RES,Tab,exudate_difference_stains,exudate_strains,pal)
########################### Area clustering #########################

# Get median per time point and integrate area
Dat.area <- aggregate(OD600 ~ Strain + condition + hrs,
                      data = All.filtered, FUN = median)
Dat.area <- aggregate(OD600 ~ Strain + condition, data = Dat.area, FUN = sum)

# Cluster
Strain.auc <- reshape2::acast(formula =  Strain ~ condition,
                    data = Dat.area, value.var = "OD600")
devtools::use_data(Strain.auc,pkg = "~/rhizogenomics/src/trunk/PGCA/PEBG/",
                   overwrite = FALSE)

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

# Plot clustering
head(Strain.auc)
pal <- colorRampPalette(colors = c("#8e0152","#de77ae","#f7f7f7","#7fbc41","#276419"))
# pal <- colorRampPalette(colors = c("#67a9cf","#ffffbf","#ef8a62"))
# pal <- colorRampPalette(colors = c("#2c7bb6","#ffffbf","#d7191c"))

svglite::svglite("auc_heatmap_gplots.svg",width = 4, height = 10)
gplots::heatmap.2(Strain.auc, Rowv = area.clus$order,
          RowSideColors = rainbow(10)[ groups[ row.names(Strain.auc) ] ],
          cexCol = 0.8, col = pal(10),trace = "none")
dev.off()

hist(1:10,col = rainbow(10)[group.order],breaks = 0:10)

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

# subset(Strain.auc, Strain %in% c("CL21","1","10","Ecoli"))
p1 <- ggplot(Strain.auc,
             aes(x = condition, y = Strain)) +
  facet_grid(Group ~ condition, scales = "free", space = "free") +
  geom_tile(aes(fill = Zscore)) +
  scale_fill_gradientn(colours = pal(10)) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold"),
        strip.text.x = element_text(angle = 90))
p1
ggsave("auc_heatmap_ggplot.svg",p1, width = 4, height = 8)

# Plot groups
for(g in unique(Strain.auc$Group)){
  # g <- unique(Strain.auc$Group)[3]

  Dat <- subset(Strain.auc, Group == g)
  p1 <- ggplot(Dat,aes(x = condition, y = Zscore)) +
    geom_boxplot() +
    AMOR::theme_blackbox +
    theme(axis.text = element_text(size = 12))
  #p1
  filename <- paste("group",g,".svg",sep = "")
  ggsave(filename = filename,plot = p1,width = 5, height = 5)
}
