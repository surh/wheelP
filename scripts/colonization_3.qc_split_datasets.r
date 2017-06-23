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
library(AMOR)

# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")
# setwd("/home/sur/rhizogenomics/experiments/2017/today4")

data("wheelP.full")
Dat <- wheelP.full
rm(wheelP.full)


p1 <- phylogram(collapse_by_taxonomy(Dat,level = 1),facet = ~ Fraction,nrow.legend = 7)
p1
ggsave("taxatype_phylogram.svg", p1, width = 8, height = 4)


p1 <- plotgg_var(Dat,var.name = "Depth", x = "Fraction") +
  scale_y_log10()
p1
ggsave("depth.svg", p1, width = 5, height = 5)

# Normalize by blank or no bacteria
blank.norm <- round(rowMeans(subset(Dat,Fraction == "BLANK")$Tab))
nobact.norm <- round(rowMeans(subset(Dat, Bacteria == "No Bacteria")$Tab))

#Dat$Tab <- t(t(Dat$Tab) - blank.norm)
Dat$Tab <- t(t(Dat$Tab) - nobact.norm)
Dat$Tab[ Dat$Tab < 0 ] <- 0
Dat$Map$Depth.norm <- colSums(Dat$Tab)
Dat <- clean(Dat)
p1 <- plotgg_var(Dat,var.name = "Depth.norm", x = "Fraction") +
  scale_y_log10()
p1
ggsave("depth.norm.svg", p1, width = 5, height = 5)

# Remove contams
Dat <- remove_taxons(Dat,taxons = taxa(Dat)[ grep("contaminant",taxa(Dat)) ] )
Dat$Map$Usable <- colSums(Dat$Tab)
Dat <- clean(Dat)
p1 <- plotgg_var(Dat,var.name = "Usable", x = "Fraction") +
  scale_y_log10()
p1
ggsave("usable.svg", p1, width = 5, height = 5)

# Clean non-relevant samples
Dat <- subset(Dat, Fraction != "ABSENT",clean = TRUE,drop = TRUE)
Dat <- subset(Dat, Fraction != "NA.",clean = TRUE,drop = TRUE)
#Dat <- subset(Dat, Fraction != "BLANK",clean = TRUE,drop = TRUE)
Dat <- clean(Dat)

plotgg_var(Dat,var.name = "Usable", x = "Inoculated", col = "Bacteria") +
  scale_y_log10()

plotgg_var(Dat,var.name = "Usable", x = "Inoculated", col = "Plate") +
  scale_y_log10()


Dat$Tax$Block <- NA
Dat$Tax$Block[ Dat$Tax$ID %in% c("CL21","217","CL59","278",
                                 "CL96","CL28","113","215","160")] <- "P1"
# Missing CL4"
Dat$Tax$Block[ Dat$Tax$ID %in% c("499","340","137A","137B","385",
                                 "CL81", "CL71","474","371")] <- "P2"
# missing 229
Dat$Tax$Block[ Dat$Tax$ID %in% c("210","28","146","186","74",
                                 "111A","111B","121")] <- "P3"
# missing CL17
Dat$Tax$Block[ Dat$Tax$ID %in% c("77","CL87","122","10","CL155",
                                 "329","CL149","143")] <- "I1"
# missing CL151
Dat$Tax$Block[ Dat$Tax$ID %in% c("275","283","CL151","30A","316",
                                 "211","295","131","36")] <- "I2"
# missing 240
Dat$Tax$Block[ Dat$Tax$ID %in% c("370","20","317","363","33",
                                 "168","165")] <- "I3"
# Missing 350
Dat$Tax$Block[ Dat$Tax$ID %in% c("233","339","375","123","CL89",
                                 "259","214","72")] <- "N1"
# Missing CL9
Dat$Tax$Block[ Dat$Tax$ID %in% c("494","345A","345B","224","69",
                                 "218","CL11","376","CL32") ] <- "N2"
# Missing 88, 59
Dat$Tax$Block[ Dat$Tax$ID %in% c("384","267","13A","13B","1B",
                                 "CL52","351")] <- "N3"
Dat$Tax$Block <- factor(Dat$Tax$Block,levels = c("P1","P2","P3",
                                                 "I1","I2","I3",
                                                 "N1","N2","N3"))

Dat$Tax[ is.na(Dat$Tax$Block), ]
table(Dat$Tax$Block)

# Heatmap
p1 <- heatgg(Dat)
head(p1$data)
dat <- p1$data
dat$Block <- Dat$Tax[ as.character(dat$Taxon), "Block" ]
dat$Taxon <- factor(dat$Taxon,
                    levels = as.character(Dat$Tax$ID[ order(Dat$Tax$Block)]))

p1 <- ggplot(dat,aes(x = Taxon, y = SAMPLEID)) +
  facet_grid(Bacteria ~ ., scales = "free") +
  geom_tile(aes(fill = Abundance)) +
  # scale_fill_continuous(trans = "log2") +
  scale_fill_gradientn(colours = c("darkblue","blue","yellow","orange","red"),
                       trans = "log2") +
  # scale_fill_gradient2(high = c("red","orange"),low = c("darkblue","blue"),
  #                      mid = "yellow",midpoint = 256, trans = "log2") +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
p1
ggsave("log2_abundance.png", p1, width = 12, height = 8)

p1 <- p1 + facet_grid(Bacteria ~ Block, scales = "free")
ggsave("log2_abundance_faceted.png", p1, width = 12, height = 8)


dat2 <- aggregate(Abundance ~ Block + Bacteria, data = p1$data,FUN = sum)
p1 <- ggplot(dat2,aes(x = Block, y = Bacteria)) +
  #facet_grid(Bacteria ~ Block, scales = "free") +
  geom_tile(aes(fill = Abundance)) +
  scale_fill_continuous(trans = "log2",low = "white",high = "darkred") +
  theme(strip.text.y = element_text(angle = 0))
p1
ggsave("log2_abundance_byblock.png", p1, width = 6, height = 6)


# read dist
dist.seq <- read.table("ref.substr.trim350.aln.pim", comment.char = "#")
dist.seq$V1 <- NULL
row.names(dist.seq) <- as.character(dist.seq$V2)
dist.seq$V2 <- NULL
colnames(dist.seq) <- row.names(dist.seq)

p1 <- plot_cross_mapping(Dat$Tax,dist.seq,id = 98)
p1
ggsave("cross_mapping98.svg",p1,width = 12, height = 4)
p1 <- plot_cross_mapping(Dat$Tax,dist.seq,id = 100)
p1
ggsave("cross_mapping100.svg",p1,width = 12, height = 4)

# Write counts per taxonomic level per group
write.qiime(pool_samples(collapse_by_taxonomy(Dat,level = 7),groups = "Bacteria",FUN = sum)$Tab,
            file = "genus_mapall.txt")
write.qiime(pool_samples(collapse_by_taxonomy(Dat,level = 6),groups = "Bacteria",FUN = sum)$Tab,
            file = "family_mapall.txt")
write.qiime(pool_samples(collapse_by_taxonomy(Dat,level = 5),groups = "Bacteria",FUN = sum)$Tab,
            file = "order_mapall.txt")

### Compare to syncomP ##############

sp.wheel <- read.am("syncomP.tables/otutab.sp.wheelmap.trim350.txt",
                    taxonomy = "taxonomy", format = "qiime",simplify = FALSE)
sp.wheel$Map <- data.frame(ID = colnames(sp.wheel$Tab), row.names = colnames(sp.wheel$Tab))
sp.wheel <- heatgg(normalize(sp.wheel,norm = colSums(sp.wheel$Tab)),
                   trans = "sqrt")$data
sp.wheel$Sample.group <- "SyncomP"
sp.wheel$Map <- "wheelP"

wheel <- heatgg(normalize(Dat, norm = colSums(Dat$Tab)), trans = "sqrt")$data
wheel$Sample.group <- "wheelP"
wheel$Map <- "wheelP"
wheel <- wheel[ ,c("ID","SAMPLEID","Taxon","Abundance","Sample.group","Map")]
sp.wheel <- sp.wheel[ ,c("ID","SAMPLEID","Taxon","Abundance","Sample.group","Map")]

dat <- rbind(wheel,sp.wheel)
dat$Taxon <- factor(dat$Taxon, levels = Dat$Tax$ID[ order(Dat$Tax$Block)])

p1 <- ggplot(dat,aes(x = Taxon, y = SAMPLEID)) +
  facet_grid(Sample.group ~ ., scales = "free", space = "free") +
  geom_tile(aes(fill = Abundance)) +
  #scale_fill_continuous(trans = "log2") +
  scale_fill_gradient2(low = c("darkblue","blue"),mid = "yellow",
                       high = c("orange","red"), midpoint = 1,
                       trans = "log2") +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, color = "black"))
# p1
ggsave("heatmap_syncomP_wheelP.png",p1,width = 10, height = 10, dpi = 200)

p1 <- ggplot(data.frame(acast(data = aggregate(Abundance ~ Taxon + Sample.group,
                                               data = dat,
                                               FUN = mean),
                              formula = Taxon ~ Sample.group,
                              value.var = "Abundance")),
              aes(x = SyncomP, y = wheelP)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point()
p1
ggsave("cor.syncomp.wheelp.svg",p1,width = 4,height = 4)
cor.test(p1$data$SyncomP,p1$data$wheelP)
cor.test(log10(p1$data$SyncomP),log10(p1$data$wheelP))

######################




p1 <- heatgg(normalize(Dat,"Usable"))
head(p1$data)
dat <- p1$data
dat$Block <- Dat$Tax[ as.character(dat$Taxon), "Block" ]


p1 <- heatgg(subset(normalize(Dat,"Usable"), Fraction == "Inoculum"))
head(p1$data)
dat <- p1$data
dat$Block <- Dat$Tax[ as.character(dat$Taxon), "Block" ]


p1 <- ggplot(dat,aes(x = Taxon, y = SAMPLEID)) +
  facet_grid(Bacteria ~ Block, scales = "free") +
  geom_tile(aes(fill = Abundance)) +
  scale_fill_continuous(trans = "log2") +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90))
p1
ggsave("norm_abundance_inoc.svg", p1, width = 6, height = 6)

dat2 <- aggregate(Abundance ~ Block + Bacteria, data = p1$data,FUN = mean)
p1 <- ggplot(dat2,aes(x = Block, y = Bacteria)) +
  #facet_grid(Bacteria ~ Block, scales = "free") +
  geom_tile(aes(fill = Abundance)) +
  scale_fill_continuous(trans = "log2",low = "white",high = "darkred") +
  theme(strip.text.y = element_text(angle = 0))
p1
ggsave("log2_abundance_byblock.png", p1, width = 6, height = 6)





set.seed(429247)
# keep <- sample(x = unique(wheel$ID),
#                size = length(levels(sp.wheel$ID)),
#                replace = FALSE)
# wheel <- droplevels(subset(wheel, ID %in% keep))
wheel <- wheel[ , c("ID","SAMPLEID","Taxon","Abundance","Sample.group","Map")]
sp.wheel <- sp.wheel[ , c("ID","SAMPLEID","Taxon","Abundance","Sample.group","Map")]
dat <- rbind(wheel, sp.wheel)
dat$Taxon <- factor(dat$Taxon, levels = Dat$Tax$ID[ order(Dat$Tax$Block)])

p1 <- ggplot(dat, aes( x = Taxon, y = SAMPLEID)) +
  facet_grid(Sample.group ~ ., scales = "free") +
  geom_tile(aes(fill = Abundance)) +
  scale_fill_gradient(trans = "log2") +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_blank())
p1












aggregate(Abundance ~ G1, data = subset(dat, Block == "G1" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ G2, data = subset(dat, Block == "G2" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ G3, data = subset(dat, Block == "G3" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ N1, data = subset(dat, Block == "N1" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ N2, data = subset(dat, Block == "N2" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ N3, data = subset(dat, Block == "N3" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ B1, data = subset(dat, Block == "B1" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ B2, data = subset(dat, Block == "B2" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ B3, data = subset(dat, Block == "B3" & Bacteria != "No Bacteria"),FUN = mean)




aggregate(Abundance ~ B1, data = subset(dat, Block == "B2" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ B2, data = subset(dat, Block == "B2" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ B3, data = subset(dat, Block == "B2" & Bacteria != "No Bacteria"),FUN = mean)

aggregate(Abundance ~ G1, data = subset(dat, Block == "B2" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ G2, data = subset(dat, Block == "B2" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ G3, data = subset(dat, Block == "B2" & Bacteria != "No Bacteria"),FUN = mean)

aggregate(Abundance ~ N1, data = subset(dat, Block == "B2" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ N2, data = subset(dat, Block == "B2" & Bacteria != "No Bacteria"),FUN = mean)
aggregate(Abundance ~ N3, data = subset(dat, Block == "B2" & Bacteria != "No Bacteria"),FUN = mean)




