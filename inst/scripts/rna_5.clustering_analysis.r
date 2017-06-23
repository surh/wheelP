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

# setwd("~/rhizogenomics/experiments/2017/today7//")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")

# Load count table
data(wheelP.rna)
Dat <- wheelP.rna
rm(wheelP.rna)

# Calculate RPKM
gene.lengths <- read.table("~/rhizogenomics/data/phosphate/gene_lengths.txt",
                           row.names = 1, header = TRUE)
Dat.norm <- create_dataset(Tab = rpkm(x = Dat$Tab,
                                      gene.length = gene.lengths[ taxa(Dat), ]),
                           Map = Dat$Map,
                           Tax = Dat$Tax)
sd <- apply(Dat.norm$Tab,1,sd)
Dat.norm <- remove_taxons(Dat.norm,
                          taxons = setdiff(taxa(Dat.norm),
                                           names(sd[sd >= 1])))
Dat.norm <- create_dataset(t(scale(t(Dat.norm$Tab))),Dat.norm$Map,Dat.norm$Tax)
Dat.norm

# Subset genes that have low RPKM standard deviation
# Dat.norm <- remove_taxons(Dat.norm,
#                           taxons = setdiff(taxa(Dat.norm),
#                                            names(sd[sd >= 1])))

# Cluster
sam.dis <- dist(Dat.norm$Tab)
sam.clus <- hclust(sam.dis)
rm(sam.dis)
gc()
plot(sam.clus)
clusters <- cutree(sam.clus, k = 7)
table(clusters)
cluster.order <- unique(clusters[ sam.clus$labels[sam.clus$order ] ])
strain.order <- sam.clus$labels[sam.clus$order ]

# temp2 <- pool_samples.default(Dat.sub$Tab,
#                               groups = interaction(Dat.sub$Map$Bacteria,
#                                                    Dat.sub$Map$Phosphate),
#                               FUN = mean)
# temp2 <- temp2$Tab
# temp2 <- dist(temp2)
# temp2 <- hclust(temp2)
#plot(temp2)

# Prepare for heatmap
p1 <-  heatgg(Dat.norm,facet = ~ Phosphate + Bacteria,trans = "identity")

dat <- p1$data
rm(p1)
gc()
# dat <- subset(dat,Taxon %in% levels(dat$Taxon)[1:100])
dat$Abundance[ dat$Abundance > 3 ] <- 3
dat$Abundance[ dat$Abundance < -3 ] <- -3
dat$Cluster <- factor(clusters[ dat$Taxon ], levels = cluster.order)
dat$Taxon <- factor(dat$Taxon,levels = strain.order)
levels(dat$Cluster) <- paste("c",1:length(levels(dat$Cluster)),sep = "")

clusters <- tapply(X = dat$Cluster,dat$Taxon, function(x) unique(as.character(x)))
clusters <- data.frame(Gene = names(clusters), Cluster = clusters)
clusters$Cluster <- factor(clusters$Cluster,levels = levels(dat$Cluster))
write.table(clusters,file = "cluster_assignment.txt",sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

p1 <- ggplot(dat, aes(x = Taxon, y = Sample)) +
  facet_grid(Phosphate + Bacteria ~ Cluster, scales = "free",space = "free") +
  geom_tile(aes(fill = Abundance)) +
  # scale_fill_gradient2(na.value = NA) +
  scale_fill_gradient2(low =  c("#8e0152","#de77ae"),
                       mid = "#f7f7f7",
                       high = c("#7fbc41","#276419"),
                       midpoint = 0,
                       na.value = "#404040",
                       guide = guide_colorbar(title = "z-score")) +

  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        panel.background = element_blank(),
        panel.spacing.y = unit(0,"lines"),
        axis.ticks = element_blank())
# p1
ggsave("heatmap.png",p1,width = 10,height = 10)
rm(p1,dat,sam.clus,cluster.order,gene.lengths,sd)
gc()

## heatmap pooled
Dat.norm <- pool_samples(Dat.norm$Tab,
                         groups = interaction(Dat.norm$Map$Bacteria,Dat.norm$Map$Phosphate),
                         FUN = mean)
Map <- data.frame(do.call(rbind,strsplit(samples(Dat.norm),"[.]")))
colnames(Map) <- c("Bacteria","Phosphate")
Map$Bacteria <- factor(Map$Bacteria, levels = levels(Dat$Map$Bacteria))
Map$Phosphate <- factor(Map$Phosphate, levels = levels(Dat$Map$Phosphate))
row.names(Map) <- samples(Dat.norm)
Dat.norm <- create_dataset(Tab = Dat.norm$Tab,Map = Map)


# Prepare for heatmap
p1 <-  heatgg(Dat.norm,facet = ~ Phosphate + Bacteria,trans = "identity")

dat <- p1$data
rm(p1)
gc()
# dat <- subset(dat,Taxon %in% levels(dat$Taxon)[1:100])
dat$Abundance[ dat$Abundance > 3 ] <- 3
dat$Abundance[ dat$Abundance < -3 ] <- -3
dat$Cluster <- clusters[ as.character(dat$Taxon), "Cluster"]
dat$Taxon <- factor(dat$Taxon,levels = strain.order)

p1 <- ggplot(dat, aes(x = Taxon, y = SAMPLEID)) +
  facet_grid(Phosphate + Bacteria ~ Cluster, scales = "free",space = "free") +
  geom_tile(aes(fill = Abundance)) +
  # scale_fill_gradient2(na.value = NA) +
  scale_fill_gradient2(low =  c("#8e0152","#de77ae"),
                       mid = "#f7f7f7",
                       high = c("#7fbc41","#276419"),
                       midpoint = 0,
                       na.value = "#404040",
                       guide = guide_colorbar(title = "z-score")) +

  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        panel.background = element_blank(),
        panel.spacing.y = unit(0,"lines"),
        axis.ticks = element_blank())
# p1
ggsave("heatmap_pooled.png",p1,width = 10,height = 10)

