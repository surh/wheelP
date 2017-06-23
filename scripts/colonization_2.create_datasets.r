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

# devtools::document(pkg = "~/rhizogenomics/src/trunk/phosphate_code/wheelP/")
# setwd("~/rhizogenomics/experiments/2017/today4/")

# Read count tables
Tab1 <- read.am("~/rhizogenomics/experiments/2017/2017-02-06.wheel_colonization/otu_table_150903.txt",
                format = "qiime",taxonomy = "taxonomy")
# Tab1 <- remove_samples(Tab1,samples = "NA.")
Tab1

Tab2 <- read.am("~/rhizogenomics/experiments/2017/2017-02-06.wheel_colonization/otu_table_160201.txt",
                format = "qiime",taxonomy = "taxonomy")
# Tab2 <- remove_samples(Tab2,samples = "NA.")
Tab2

Tab3 <- read.am("~/rhizogenomics/experiments/2017/2017-02-06.wheel_colonization/otu_table_160909.txt",
                format = "qiime",taxonomy = "taxonomy")
# Tab3 <- remove_samples(Tab3,samples = "NA.")
Tab3

Tab4 <- read.am("~/rhizogenomics/experiments/2017/2017-02-06.wheel_colonization/otu_table_170118.txt",
                format = "qiime",taxonomy = "taxonomy")
# Tab4 <- remove_samples(Tab4,samples = "NA.")
Tab4

Dat <- combine_datasets(list(Tab1,Tab2,Tab3,Tab4))
Dat

# Read mapping files
Map1 <- read.table("~/rhizogenomics/experiments/2016/2016-12-06.wheelP_demultref/R150903_demultiplex_map.txt",
                   sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
row.names(Map1) <- Map1$rnaID
head(Map1)

Map2 <- read.table("~/rhizogenomics/experiments/2016/2016-12-06.wheelP_demultref/R160201_demultiplex_map.txt",
                   sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
row.names(Map2) <- Map2$rnaID
head(Map2)

Map3 <- read.table("~/rhizogenomics/experiments/2016/2016-12-06.wheelP_demultref/R160909_demultiplex_map.txt",
                   sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
absent <- Map3[ Map3$rnaID == "ABSENT", ]
Map3 <- Map3[ Map3$rnaID != "ABSENT", ]
row.names(Map3) <- as.character(Map3$rnaID)
head(Map3)

Map4 <- read.table("~/rhizogenomics/experiments/2016/2016-12-06.wheelP_demultref/R170118_demultiplex_map.txt",
                   sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
row.names(Map4) <- Map4$rnaID
head(Map4)

# give priority to run 4 over run 2, which had very few reads (less than 10 for
# most samples)
redundant <- intersect(row.names(Map2),row.names(Map4))
Map2 <- Map2[ setdiff(row.names(Map2),redundant), ]

Map <- rbind(Map1,Map2,Map3,Map4)
row.names(Map)

## Add ABSENT and NA
setdiff(samples(Dat),row.names(Map))
absent <- absent[1,]
absent
absent$ID <- "ABSENT"
absent$Well <- "ABSENT"
absent$Sample_ID <- "ABSENT"
absent$Sample_Name <- "ABSENT"
absent$Sample_Plate <- "ABSENT"
absent$Sample_Well <- "ABSENT"
absent$I7_Index_ID <- "ABSENT"
absent$index <- "ABSENT"
absent
na <- absent
na
na[ na == "ABSENT" ] <- "NA."
na$Plate <- "NA."
na$Run  <- "NA."
na$Barcode2 <- "NA."
na$Frameshift <- "NA."
na
add <- rbind(absent,na)
row.names(add) <- add$rnaID
add

Map <- rbind(Map,add)
row.names(Map)

## Add block variables
Map$Genotype <- factor(Map$Genotype,
                       levels = c("Col-0", "Inoculum",
                                  "T0Agar","BLANK","ABSENT","NA." ))
Map$Bacteria[ Map$Bacteria == "B3G1" ] <- "G1B3"
colnames(Map)
Map$P1 <- 0
Map$P1[ grep("G1",Map$Bacteria) ] <- 1
Map$P2 <- 0
Map$P2[ grep("G2",Map$Bacteria) ] <- 1
Map$P3 <- 0
Map$P3[ grep("G3",Map$Bacteria) ] <- 1
Map$I1 <- 0
Map$I1[ grep("N1",Map$Bacteria) ] <- 1
Map$I2 <- 0
Map$I2[ grep("N2",Map$Bacteria) ] <- 1
Map$I3 <- 0
Map$I3[ grep("N3",Map$Bacteria) ] <- 1
Map$N1 <- 0
Map$N1[ grep("B1",Map$Bacteria) ] <- 1
Map$N2 <- 0
Map$N2[ grep("B2",Map$Bacteria) ] <- 1
Map$N3 <- 0
Map$N3[ grep("B3",Map$Bacteria) ] <- 1
colnames(Map)
Map$Bacteria <- factor(Map$Bacteria,levels = c("No Bacteria",
                                               "G1G2",
                                               "G2G3",
                                               "G1G3",
                                               "G1N1",
                                               "G1N2",
                                               "G1N3",
                                               "G2N1",
                                               "G3N1",
                                               "N1N2",
                                               "N1N3",
                                               "N2N3",
                                               "N1B1",
                                               "N3B1",
                                               "G1B1",
                                               "G1B2",
                                               "G2B1",
                                               "G2B2",
                                               "G3B1",
                                               "G3B2",
                                               "N1B2",
                                               "B1B2",
                                               "B1B3",
                                               "B2B3",
                                               "N1B3",
                                               "G1B3",
                                               "G2B3",
                                               "G3B3",
                                               "BLANK",
                                               "ABSENT",
                                               "NA."))
levels(Map$Bacteria) <- c("No Bacteria",
                          "P1P2",
                          "P2P3",
                          "P1P3",
                          "P1I1",
                          "P1I2",
                          "P1I3",
                          "P2I1",
                          "P3I1",
                          "I1I2",
                          "I1I3",
                          "I2I3",
                          "I1N1",
                          "I3N1",
                          "P1N1",
                          "P1N2",
                          "P2N1",
                          "P2N2",
                          "P3N1",
                          "P3N2",
                          "I1N2",
                          "N1N2",
                          "N1N3",
                          "N2N3",
                          "I1N3",
                          "P1N3",
                          "P2N3",
                          "P3N3",
                          "BLANK",
                          "ABSENT",
                          "NA.")


Map$Pre.Pi <- NA
Map$Pre.Pi[ grep("^-P",Map$Phosphate) ] <- "-Pi,0.5%Suc"
Map$Pre.Pi[ grep("^[+]P",Map$Phosphate) ] <- "+Pi,0.5%Suc"
Map$Pre.Pi[ grep("BLANK",Map$Phosphate) ] <- "BLANK"
Map$Pre.Pi[ grep("Inoculum",Map$Phosphate) ] <- "Inoculum"
Map$Pre.Pi[ grep("ABSENT",Map$Phosphate) ] <- "ABSENT"
Map$Pre.Pi[ grep("NA.",Map$Phosphate) ] <- "NA."
Map$Pos.Pi <- NA
Map$Pos.Pi[ grep("30uM$",Map$Phosphate) ] <- "30 uM,0%Suc"
Map$Pos.Pi[ grep("100uM$",Map$Phosphate) ] <- "100 uM,0%Suc"
Map$Pos.Pi[ grep("BLANK",Map$Phosphate) ] <- "BLANK"
Map$Pos.Pi[ grep("Inoculum",Map$Phosphate) ] <- "Inoculum"
Map$Pos.Pi[ grep("ABSENT",Map$Phosphate) ] <- "ABSENT"
Map$Pos.Pi[ grep("NA.",Map$Phosphate) ] <- "NA."
Map$Pre.Pi <- factor(Map$Pre.Pi, levels = c("+Pi,0.5%Suc","-Pi,0.5%Suc",
                                            "BLANK","Inoculum","ABSENT",
                                            "NA."))
Map$Pos.Pi <- factor(Map$Pos.Pi, levels = c("100 uM,0%Suc", "30 uM,0%Suc",
                                            "BLANK","Inoculum","ABSENT",
                                            "NA."))
colnames(Map)
Map$Replicate[ Map$Fraction == "BLANK" ] <- "BLANK"
colnames(Map)
Map$Fraction <- factor(Map$Fraction, levels = c("Agar","Root","Inoculum",
                                                "T0Agar","BLANK","ABSENT",
                                                "NA."))
colnames(Map)

Map$Inoculated <- as.character(Map$Bacteria)
Map$Inoculated[ !(Map$Inoculated %in% c("No Bacteria", "BLANK", "ABSENT", "NA.")) ] <- "+Bacteria"
Map$Inoculated <- factor(Map$Inoculated,
                         levels = c("No Bacteria","+Bacteria","BLANK", "ABSENT", "NA."))

Map[ rowSums(is.na(Map)) > 0, ]


setdiff(samples(Dat),row.names(Map))
Map <- Map[ samples(Dat), ]
Map <- droplevels(Map)

Dat <- create_dataset(Tab = Dat$Tab,Map = Map,Tax = Dat$Tax)
Dat$Map$Depth <- colSums(Dat$Tab)
wheelP.full <- Dat
wheelP.full
devtools::use_data(wheelP.full,pkg = "~/rhizogenomics/src/trunk/phosphate_code/wheelP/",
                   overwrite = TRUE)

