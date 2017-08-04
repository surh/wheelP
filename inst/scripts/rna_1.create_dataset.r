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

# setwd("~/rhizogenomics/experiments/2017/today5//")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")

# Read data
Tab <- read.table("~/rhizogenomics/data/synthetic/wheel_phosphate/RNAseq/SC_gene_counts.txt",
                  sep = "\t", row.names = 1, header = TRUE)
Map <- read.table("~/rhizogenomics/data/synthetic/wheel_phosphate/RNAseq/2016-03-14.map.txt",
                  header = TRUE, stringsAsFactors = FALSE)
row.names(Map) <- as.character(Map$ID)
head(Map)

Map$Well <- factor(Map$Well)
Map$Plate <- factor(Map$Plate)
Map$Genotype <- factor(Map$Genotype)
Map$Bacteria[ Map$Bacteria == "-" ] <- "No Bacteria"
Map$Bacteria[ Map$Bacteria == "SC1" ] <- "P1P2"
Map$Bacteria[ Map$Bacteria == "SC2" ] <- "P2P3"
Map$Bacteria[ Map$Bacteria == "SC3" ] <- "P3I1"
Map$Bacteria[ Map$Bacteria == "SC4" ] <- "I1I2"
Map$Bacteria[ Map$Bacteria == "SC5" ] <- "I2I3"
Map$Bacteria[ Map$Bacteria == "SC6" ] <- "I3N1"
Map$Bacteria[ Map$Bacteria == "SC7" ] <- "N1N2"
Map$Bacteria[ Map$Bacteria == "SC8" ] <- "N2N3"
Map$Bacteria[ Map$Bacteria == "SC9" ] <- "P1N3"
Map$Bacteria[ Map$Bacteria == "SC10" ] <- "P1I1"
Map$Bacteria[ Map$Bacteria == "SC11" ] <- "P2I1"
Map$Bacteria[ Map$Bacteria == "SC12" ] <- "P1P3"
Map$Bacteria[ Map$Bacteria == "SC13" ] <- "P2N3"
Map$Bacteria[ Map$Bacteria == "SC14" ] <- "P3N3"
Map$Bacteria <- factor(Map$Bacteria,
                       levels = c("No Bacteria","P1P2","P2P3","P3I1","I1I2","I2I3",
                                  "I3N1","N1N2","N2N3","P1N3","P1P3","P1I1","P2I1",
                                  "P2N3","P3N3"))

Map$P1 <- 0
Map$P1[ grep("P1",Map$Bacteria) ] <- 1
Map$P2 <- 0
Map$P2[ grep("P2",Map$Bacteria) ] <- 1
Map$P3 <- 0
Map$P3[ grep("P3",Map$Bacteria) ] <- 1
Map$I1 <- 0
Map$I1[ grep("I1",Map$Bacteria) ] <- 1
Map$I2 <- 0
Map$I2[ grep("I2",Map$Bacteria) ] <- 1
Map$I3 <- 0
Map$I3[ grep("I3",Map$Bacteria) ] <- 1
Map$N1 <- 0
Map$N1[ grep("N1",Map$Bacteria) ] <- 1
Map$N2 <- 0
Map$N2[ grep("N2",Map$Bacteria) ] <- 1
Map$N3 <- 0
Map$N3[ grep("N3",Map$Bacteria) ] <- 1

Map$Phosphate[ Map$Phosphate == "-P30uM" ] <- "minusP_30uM"
Map$Phosphate[ Map$Phosphate == "+P30uM" ] <- "plusP_30uM"
Map$Phosphate[ Map$Phosphate == "-P100uM" ] <- "minusP_100uM"
Map$Phosphate[ Map$Phosphate == "+P100uM" ] <- "plusP_100uM"
Map$Phosphate <- factor(Map$Phosphate,
                        levels = c("plusP_100uM","minusP_100uM","plusP_30uM","minusP_30uM"))
Map$Replicate <- factor(Map$Replicate)
Map$Experiment <- factor(Map$Experiment)
Map$Extraction <- factor(Map$Extraction)
Map <- Map[ colnames(Tab), ]

# Create Dataset
Dat <- create_dataset(Tab = Tab,Map = Map)
Dat <- clean(Dat)

Dat.pca <- PCA(Dat, cor = TRUE)
p1 <- plotgg(Dat.pca,col = "Phosphate", point_size = 1)
p1
ggsave("PCA_phosphate.svg",p1,width = 4, height = 4)

p1 <- plotgg(Dat.pca,col = "Bacteria", point_size = 1)
p1
ggsave("PCA_Bacteria.svg",p1,width = 4, height = 4)

p1 <- plotgg(Dat.pca,col = "Replicate", point_size = 1)
p1
ggsave("PCA_Replicate.svg",p1,width = 4, height = 4)

p1 <- plotgg(Dat.pca,col = "Experiment", point_size = 1)
p1
ggsave("PCA_Experiment.svg",p1,width = 4, height = 4)

p1 <- plotgg(Dat.pca,col = "Extraction", point_size = 1)
p1
ggsave("PCA_Extraction.svg",p1,width = 4, height = 4)

p1 <- plotgg(Dat.pca,col = "Concentration", point_size = 1)
p1
ggsave("PCA_Concentration.svg",p1,width = 4, height = 4)

p1 <- plotgg(Dat.pca,col = "A260.280", point_size = 1)
p1
ggsave("PCA_Concentration.svg",p1,width = 4, height = 4)

p1 <- plotgg(Dat.pca,col = "A260.230", point_size = 1)
p1
ggsave("PCA_A260.230.svg",p1,width = 4, height = 4)

# save(Dat,file = "wheel_dataset.rdat")
wheelP.rna <- Dat
devtools::use_data(wheelP.rna,
                   pkg = "~/rhizogenomics/github/wheelP/",
                   overwrite = TRUE)
