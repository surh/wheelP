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

#setwd("~/rhizogenomics/experiments/2017/today3/")

rna <- read.table("~/rhizogenomics/data/synthetic/wheel_phosphate/RNAseq/2016-03-14.map.txt",
                  sep = "\t", header = TRUE)
rna[1:5,1:5]
rna$Plate <- paste("Plate", rna$Plate, sep = "")
# Reading master demultiplex file for 16S
wheel <- read.table("~/rhizogenomics/data/synthetic/wheel_phosphate/RNAseq/wheel.txt",
                    sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE)
wheel[1:5,1:5]

# Match old IDs from RNA
ids <- as.character(rna$ID)
names(ids) <- paste(rna$Plate,rna$Well,sep = ".")
ids
wheel$rnaID <- ids[ paste(wheel$Plate,wheel$Well,sep = ".") ]
wheel$rnaID

# Remove space in phosphate
wheel$Phosphate <- gsub("^ ",replacement = "", wheel$Phosphate)

# Change community names
wheel$Bacteria[ wheel$Bacteria == "-" ] <- "No Bacteria"
wheel$Bacteria[ wheel$Bacteria == " -" ] <- "No Bacteria"
wheel$Bacteria[ wheel$Bacteria == "SC1" ] <- "P1P2"
wheel$Bacteria[ wheel$Bacteria == "SC2" ] <- "P2P3"
wheel$Bacteria[ wheel$Bacteria == "SC3" ] <- "P3I1"
wheel$Bacteria[ wheel$Bacteria == "SC4" ] <- "I1I2"
wheel$Bacteria[ wheel$Bacteria == "SC5" ] <- "I2I3"
wheel$Bacteria[ wheel$Bacteria == "SC6" ] <- "I3N1"
wheel$Bacteria[ wheel$Bacteria == "SC7" ] <- "N1N2"
wheel$Bacteria[ wheel$Bacteria == "SC8" ] <- "N2N3"
wheel$Bacteria[ wheel$Bacteria == "SC9" ] <- "P1N3"
wheel$Bacteria[ wheel$Bacteria == "SC10" ] <- "P1I1"
wheel$Bacteria[ wheel$Bacteria == "SC11" ] <- "P2I1"
wheel$Bacteria[ wheel$Bacteria == "SC12" ] <- "P1P3"
wheel$Bacteria[ wheel$Bacteria == "SC13" ] <- "P2N3"
wheel$Bacteria[ wheel$Bacteria == "SC14" ] <- "P3N3"
ftable(Plate ~ Bacteria, wheel)

# Check problematic samples
subset(wheel, SampleID == "SC0073")
subset(wheel, SampleID == "SC0074")
subset(wheel, SampleID == "SC0075")
subset(wheel, SampleID == "SC0073-75")

# Rename agar sample missing one of the three
wheel$SampleID [ wheel$SampleID == "SC0112-113" ] <- "SC0112-114"

# Renname samples with one zero to avoid ambiguity
wheel$SampleID [ wheel$SampleID == "SC01" ] <- "SC1001"
wheel$SampleID [ wheel$SampleID == "SC02" ] <- "SC1002"
wheel$SampleID [ wheel$SampleID == "SC03" ] <- "SC1003"
wheel$SampleID [ wheel$SampleID == "SC01-3" ] <- "SC1001-1003"

# Create dataset to match agar samples to root samples
dat <- subset(wheel, Fraction == "Agar")
dat <- data.frame(Start = NA, End = NA,
                  Experiment = dat$Experiment,
                  Bacteria = dat$Bacteria,
                  Phosphate = dat$Phosphate,
                  SampleID = dat$SampleID,
                  ID = paste(dat$SampleID,
                             dat$Experiment,
                             dat$Bacteria,
                             dat$Phosphate,
                             sep = "."))
row.names(dat) <- as.character(dat$ID)
head(dat)
dat[ , c("Start", "End") ] <- do.call(rbind,
                                      lapply(strsplit(gsub(pattern = "^SC",
                                                           replacement = "",
                                                           as.character(dat$SampleID)),
                                                      "-"),
                                             as.numeric))
dat$AgarGroup <- paste("AG",1:nrow(dat), sep = "")
head(dat)

# Get sample number for roots
wheel$root_sample_num <- as.numeric(gsub(pattern = "^SC", replacement = "", x = wheel$SampleID))
head(wheel)
wheel$AgarGroup <- NA

# Add agar group for root samples
for(i in 1:nrow(wheel)){
  # i <- 1
  # i <- 251
  # i <- 61
  # wheel[ i, ]

  if(wheel$Fraction[i] == "Agar"){
    mydat <- subset(dat, SampleID == wheel$SampleID[i])
    mydat <- subset(mydat, Experiment == wheel$Experiment[i])
    mydat <- subset(mydat, Phosphate ==  wheel$Phosphate[i])
    mydat <- subset(mydat, Bacteria ==  wheel$Bacteria[i])
    #mydat
  }else if(wheel$Fraction[i] == "Root"){
    mydat <- subset(dat, Start <= wheel$root_sample_num[i] & End >= wheel$root_sample_num[i])
    mydat <- subset(mydat, Experiment == wheel$Experiment[i])
    mydat <- subset(mydat, Phosphate ==  wheel$Phosphate[i])
    mydat <- subset(mydat, Bacteria ==  wheel$Bacteria[i])
    # mydat
  }

  if(nrow(mydat) != 1){
    stop("ERROR")
  }

  wheel$AgarGroup[i] <- mydat$AgarGroup
}

colSums(is.na(wheel))
rowSums(is.na(wheel))
ftable(Fraction ~ is.na(root_sample_num), wheel)
wheel$root_sample_num <- NULL

wheel$rnaID[ is.na(wheel$rnaID) ] <- paste("SC0", 433:(432 + table(is.na(wheel$rnaID))["TRUE"]), sep = "")
max(table((wheel$rnaID)))
colSums(is.na(wheel))

# Read validation and add columns
validation <- read.table("~/rhizogenomics/data/synthetic/wheel_phosphate/RNAseq/validation.txt",
                         sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE)
head(validation)
validation$rnaID <- NA
validation$AgarGroup <- NA
all(colnames(wheel) == colnames(validation))

# Fix naming mistakes
validation$SampleID[ validation$SampleID == "BLANK " ] <- "BLANK"
validation$SampleID[ validation$SampleID == "G3B1 TO" ] <- "G3B1 T0"
validation$SampleID[ validation$SampleID == "BIB2 INOC" ] <- "B1B2 INOC"

# Read metadata
meta <- read.table("~/rhizogenomics/data/synthetic/wheel_phosphate/RNAseq/validation_legend.txt",
                   sep = "\t", header = TRUE)
row.names(meta) <- as.character(meta$No.)
head(meta)

# Add phosphate
meta$Phosphate <- "-P30uM"
validation$Phosphate[ validation$SampleID %in% row.names(meta) ] <- "-P30uM"
validation$Phosphate[ grep("T0$", validation$SampleID) ] <- "-P30uM"
validation$Phosphate[ grep("INOC$", validation$SampleID) ] <- "Inoculum"
validation$Phosphate[ validation$SampleID  == "BLANK" ] <- "BLANK"
validation$Phosphate[ validation$SampleID  == "ABSENT" ] <- "ABSENT"
validation[ is.na(validation$Phosphate), ]

# Add rnaID
start <- max(as.numeric(sub("^SC", "", wheel$rnaID))) + 1
validation$rnaID <- paste("SC0", start:(start+nrow(validation) - 1),sep = "")
validation$rnaID[ validation$SampleID  == "ABSENT" ] <- "ABSENT"
validation[ is.na(validation$rnaID), ]

# Add genotype
validation$Genotype[ validation$SampleID %in% row.names(meta) ] <- "Col-0"
validation$Genotype[ grep("T0$", validation$SampleID) ] <- "T0Agar"
validation$Genotype[ grep("INOC$", validation$SampleID) ] <- "Inoculum"
validation$Genotype[ validation$SampleID  == "BLANK" ] <- "BLANK"
validation$Genotype[ validation$SampleID  == "ABSENT" ] <- "ABSENT"
validation[ is.na(validation$Genotype), ]

# Add genotype
validation$Genotype[ validation$SampleID %in% row.names(meta) ] <- "Col-0"
validation$Genotype[ grep("T0$", validation$SampleID) ] <- "T0Agar"
validation$Genotype[ grep("INOC$", validation$SampleID) ] <- "Inoculum"
validation$Genotype[ validation$SampleID  == "BLANK" ] <- "BLANK"
validation$Genotype[ validation$SampleID  == "ABSENT" ] <- "ABSENT"
validation[ is.na(validation$Genotype), ]

# Correct experiment
validation$Experiment <- paste("Validation", validation$Replicate, sep = "")
validation$Experiment[ validation$Experiment == "ValidationNA" ] <- NA
validation$Experiment[ validation$SampleID  == "BLANK" ] <- "BLANK"
validation$Experiment[ validation$SampleID  == "ABSENT" ] <- "ABSENT"
validation[ is.na(validation$Genotype), ]

# Correct Fraction
validation$Fraction[ validation$Fraction == "agar" ] <- "Agar"
validation$Fraction[ validation$Fraction == "root" ] <- "Root"
validation$Fraction[ grep("T0$", validation$SampleID) ] <- "T0Agar"
validation$Fraction[ grep("INOC$", validation$SampleID) ] <- "Inoculum"
validation$Fraction[ validation$SampleID  == "BLANK" ] <- "BLANK"
validation$Fraction[ validation$SampleID  == "ABSENT" ] <- "ABSENT"
validation[ is.na(validation$Fraction), ]

# Add bacteria
validation$Bacteria <- as.character(meta[ validation$SampleID, "SynCom" ])
validation$Bacteria[ validation$Bacteria == "NB" ] <- "No Bacteria"
validation$Bacteria[ validation$SampleID  == "BLANK" ] <- "BLANK"
validation$Bacteria[ validation$SampleID  == "ABSENT" ] <- "ABSENT"
validation$Bacteria[ grep("^G3B1", validation$SampleID) ] <- "G3B1"
validation$Bacteria[ grep("^B1B3", validation$SampleID) ] <- "B1B3"
validation$Bacteria[ grep("^B1B2", validation$SampleID) ] <- "B1B2"
validation$Bacteria[ grep("^N1B2", validation$SampleID) ] <- "N1B2"
validation$Bacteria[ grep("^G3B2", validation$SampleID) ] <- "G3B2"
validation$Bacteria[ grep("^B2B3", validation$SampleID) ] <- "B2B3"
validation$Bacteria[ grep("^G2B2", validation$SampleID) ] <- "G2B2"
validation$Bacteria[ grep("^G1G3", validation$SampleID) ] <- "G1G3"
validation$Bacteria[ grep("^G1B2", validation$SampleID) ] <- "G1B2"
validation$Bacteria[ grep("^G3N1", validation$SampleID) ] <- "G3N1"
validation$Bacteria[ grep("^N1B3", validation$SampleID) ] <- "N1B3"
validation$Bacteria[ grep("^G1B1", validation$SampleID) ] <- "G1B1"
validation$Bacteria[ grep("^N1B1", validation$SampleID) ] <- "N1B1"
validation$Bacteria[ grep("^G3B3", validation$SampleID) ] <- "G3B3"
validation$Bacteria[ grep("^G2B1", validation$SampleID) ] <- "G2B1"
validation$Bacteria[ grep("^G2B3", validation$SampleID) ] <- "G2B3"
validation$Bacteria[ grep("^N1N3", validation$SampleID) ] <- "N1N3"
validation$Bacteria[ grep("^G1N3", validation$SampleID) ] <- "G1N3"
validation$Bacteria[ grep("^G1B3", validation$SampleID) ] <- "G1B3"
validation$Bacteria[ grep("^G1N2", validation$SampleID) ] <- "G1N2"

# validation$Bacteria[ grep("^G3B1", validation$Bacteria) ] <- "P3N1"
# validation$Bacteria[ grep("^B1B3", validation$Bacteria) ] <- "N1N3"
# validation$Bacteria[ grep("^B1B2", validation$Bacteria) ] <- "N1N2"
# validation$Bacteria[ grep("^N1B2", validation$Bacteria) ] <- "I1N2"
# validation$Bacteria[ grep("^G3B2", validation$Bacteria) ] <- "P3N2"
# validation$Bacteria[ grep("^B2B3", validation$Bacteria) ] <- "N2N3"
# validation$Bacteria[ grep("^G2B2", validation$Bacteria) ] <- "P2N2"
# validation$Bacteria[ grep("^G1G3", validation$Bacteria) ] <- "P1P3"
# validation$Bacteria[ grep("^G1B2", validation$Bacteria) ] <- "P1N2"
# validation$Bacteria[ grep("^G3N1", validation$Bacteria) ] <- "P3I1"
# validation$Bacteria[ grep("^N1B3", validation$Bacteria) ] <- "I1N3"
# validation$Bacteria[ grep("^G1B1", validation$Bacteria) ] <- "P1N1"
# validation$Bacteria[ grep("^N1B1", validation$Bacteria) ] <- "I1N1"
# validation$Bacteria[ grep("^G3B3", validation$Bacteria) ] <- "P3N3"
# validation$Bacteria[ grep("^G2B1", validation$Bacteria) ] <- "P2N1"
# validation$Bacteria[ grep("^G2B3", validation$Bacteria) ] <- "P2N3"
# validation$Bacteria[ grep("^N1N3", validation$Bacteria) ] <- "I1I3"
# validation$Bacteria[ grep("^G1N3", validation$Bacteria) ] <- "P1I3"
# validation$Bacteria[ grep("^G1B3", validation$Bacteria) ] <- "P1N3"
# validation$Bacteria[ grep("^G1N2", validation$Bacteria) ] <- "P1I2"
validation[ is.na(validation$Bacteria), ]
unique(validation$Bacteria)
validation$Bacteria <- factor(validation$Bacteria)
levels(validation$Bacteria)
levels(validation$Bacteria) <- c("ABSENT","N1N2","N1N3","N2N3",
                                 "BLANK","P1N1","P1N2","P1N3",
                                 "P1P3","P1I2","P1I3","P2N1",
                                 "P2N2","P2N3","P3N1","P3N2",
                                 "P3N3","P3I1","I1N1","I1N2",
                                 "I1N3","I1I3","No Bacteria")
validation$Bacteria <- as.character(validation$Bacteria)

# add agar group
validation$AgarGroup <- as.character(meta[ validation$SampleID, "AgarGroup" ])
validation$AgarGroup[!is.na(validation$AgarGroup)] <- paste("V",validation$Replicate[ !is.na(validation$AgarGroup) ], validation$AgarGroup[!is.na(validation$AgarGroup)],sep = "")
validation$AgarGroup[ validation$Fraction == "Inoculum" ] <- validation$rnaID[ validation$Fraction == "Inoculum" ]
validation$AgarGroup[ validation$Fraction == "BLANK" ] <- validation$rnaID[ validation$Fraction == "BLANK" ]
validation$AgarGroup[ validation$Fraction == "T0Agar" ] <- validation$rnaID[ validation$Fraction == "T0Agar" ]
validation$AgarGroup[ validation$Fraction == "ABSENT" ] <- validation$rnaID[ validation$Fraction == "ABSENT" ]
validation[ is.na(validation$AgarGroup), ]

# Correct replicate ORDER MATTERS
validation$Replicate[ validation$Fraction %in% c("Agar", "Root")] <- as.numeric(validation$SampleID[ validation$Fraction %in% c("Agar", "Root")]) %% 3
validation$Replicate[ validation$Replicate == 0 ] <- 3
validation$Replicate[ validation$Fraction == "BLANK" ] <- validation$rnaID[ validation$Fraction == "BLANK" ]
validation$Replicate[ validation$Fraction == "ABSENT" ] <- validation$rnaID[ validation$Fraction == "ABSENT" ]
validation[ is.na(validation$Replicate), ]

colSums(is.na(validation))
rowSums(is.na(validation))

## Combine
Dat <- rbind(wheel, validation)
head(Dat)

### Add blocks
Dat$P1 <- 0
Dat$P1[ grep("P1",Dat$Bacteria) ] <- 1
Dat$P2 <- 0
Dat$P2[ grep("P2",Dat$Bacteria) ] <- 1
Dat$P3 <- 0
Dat$P3[ grep("P3",Dat$Bacteria) ] <- 1
Dat$I1 <- 0
Dat$I1[ grep("I1",Dat$Bacteria) ] <- 1
Dat$I2 <- 0
Dat$I2[ grep("I2",Dat$Bacteria) ] <- 1
Dat$I3 <- 0
Dat$I3[ grep("I3",Dat$Bacteria) ] <- 1
Dat$N1 <- 0
Dat$N1[ grep("N1",Dat$Bacteria) ] <- 1
Dat$N2 <- 0
Dat$N2[ grep("N2",Dat$Bacteria) ] <- 1
Dat$N3 <- 0
Dat$N3[ grep("N3",Dat$Bacteria) ] <- 1
colnames(Dat)

### Reorder bacteria
unique(Dat$Bacteria)
Dat$Bacteria[ Dat$Bacteria == "N3P1" ] <- "P1N3"
Dat$Bacteria <- factor(Dat$Bacteria, levels = c("No Bacteria","P1P2","P2P3",
                                                "P1P3","P1I1","P1I2","P1I3",
                                                "P2I1","P3I1","I1I2","I1I3",
                                                "I2I3","I1N1","I3N1","P1N1",
                                                "P1N2","P2N1","P2N2","P3N1",
                                                "P3N2",
                                                "I1N2","N1N2","N1N3","N2N3",
                                                "I1N3","P1N3","P2N3",
                                                "P3N3","BLANK",
                                                "ABSENT"))

### Reorder Genotype
unique(Dat$Genotype)
Dat$Genotype <- factor(Dat$Genotype,
                       levels = c("Col-0", "Inoculum",
                                  "T0Agar","BLANK","ABSENT"))

### Create & order Pre/Pos Pi
Dat$Pre.Pi <- NA
Dat$Pre.Pi[ grep("^-P",Dat$Phosphate) ] <- "-Pi,0.5%Suc"
Dat$Pre.Pi[ grep("^[+]P",Dat$Phosphate) ] <- "+Pi,0.5%Suc"
Dat$Pre.Pi[ grep("BLANK",Dat$Phosphate) ] <- "BLANK"
Dat$Pre.Pi[ grep("Inoculum",Dat$Phosphate) ] <- "Inoculum"
Dat$Pre.Pi[ grep("ABSENT",Dat$Phosphate) ] <- "ABSENT"
# Dat$Pre.Pi[ grep("NA.",Dat$Phosphate) ] <- "NA."
Dat$Pos.Pi <- NA
Dat$Pos.Pi[ grep("30uM$",Dat$Phosphate) ] <- "30 uM,0%Suc"
Dat$Pos.Pi[ grep("100uM$",Dat$Phosphate) ] <- "100 uM,0%Suc"
Dat$Pos.Pi[ grep("BLANK",Dat$Phosphate) ] <- "BLANK"
Dat$Pos.Pi[ grep("Inoculum",Dat$Phosphate) ] <- "Inoculum"
Dat$Pos.Pi[ grep("ABSENT",Dat$Phosphate) ] <- "ABSENT"
# Dat$Pos.Pi[ grep("NA.",Dat$Phosphate) ] <- "NA."
Dat$Pre.Pi <- factor(Dat$Pre.Pi, levels = c("+Pi,0.5%Suc","-Pi,0.5%Suc",
                                            "BLANK","Inoculum","ABSENT"))
Dat$Pos.Pi <- factor(Dat$Pos.Pi, levels = c("100 uM,0%Suc", "30 uM,0%Suc",
                                            "BLANK","Inoculum","ABSENT"))
### Order fraction
unique(Dat$Fraction)
Dat$Fraction <- factor(Dat$Fraction, levels = c("Agar","Root","Inoculum",
                                                "T0Agar","BLANK","ABSENT"))

### Add inoculated
Dat$Inoculated <- as.character(Dat$Bacteria)
Dat$Inoculated[ !(Dat$Inoculated %in% c("No Bacteria", "BLANK", "ABSENT", "NA.")) ] <- "+Bacteria"
Dat$Inoculated <- factor(Dat$Inoculated,
                         levels = c("No Bacteria","+Bacteria","BLANK", "ABSENT"))



colnames(Dat)
table(Dat$Well, useNA = "ifany")
table(Dat$Genotype, useNA = "ifany")
table(Dat$Bacteria, useNA = "ifany")
table(Dat$Phosphate, useNA = "ifany")
table(Dat$Replicate, useNA = "ifany")
table(Dat$Experiment, useNA = "ifany")
table(Dat$Plate, useNA = "ifany")
table(Dat$Run, useNA = "ifany")
table(Dat$Barcode2, useNA = "ifany")
table(Dat$Frameshift, useNA = "ifany")
table(Dat$Sample_ID, useNA = "ifany")
table(Dat$Sample_Name, useNA = "ifany")
table(Dat$Sample_Plate, useNA = "ifany")
table(Dat$Sample_Well, useNA = "ifany")
table(Dat$I7_Index_ID, useNA = "ifany")
table(Dat$index, useNA = "ifany")
table(Dat$Fraction, useNA = "ifany")
table(Dat$rnaID, useNA = "ifany")
sort(table(Dat$rnaID, useNA = "ifany"))
table(table(Dat$rnaID, useNA = "ifany"))
table(Dat$AgarGroup, useNA = "ifany")
sort(table(Dat$AgarGroup, useNA = "ifany"))
table(table(Dat$AgarGroup, useNA = "ifany"))

# Add re-run
rerun <- subset(Dat, Plate == "Plate3" | Plate == "Plate4" )
rerun$Run <- "R170118"
Dat <- rbind(Dat,rerun)
Dat$ID <- paste("WP",1:nrow(Dat), sep = "")

table(Dat$rnaID, useNA = "ifany")
sort(table(Dat$rnaID, useNA = "ifany"))
table(table(Dat$rnaID, useNA = "ifany"))
table(Dat$AgarGroup, useNA = "ifany")
sort(table(Dat$AgarGroup, useNA = "ifany"))
table(table(Dat$AgarGroup, useNA = "ifany"))

write.table(Dat,"metadadata_full.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE,
            quote = FALSE)
Map.colonization <- Dat
devtools::use_data(Map.colonization,
                   pkg = "~/rhizogenomics/github/wheelP/",
                   overwrite = TRUE)
head(Dat)

rm(list = ls())
Tax <- read.table("~/rhizogenomics/data/synthetic/wheel_phosphate/ref_wheel/strainmap.txt",
                  sep = "\t")
head(Tax)
Taxonomy <- read.table("~/rhizogenomics/data/synthetic/wheel_phosphate/ref_wheel/taxonomy.txt",
                       sep = "\t")
head(Taxonomy)


row.names(Tax) <- as.character(Tax$V1)
colnames(Tax) <- c("ID","ID2","Type","Block","MonoID")
row.names(Taxonomy) <- as.character(Taxonomy$V1)
Tax$Taxonomy <- Taxonomy[ row.names(Tax), "V2"]

# Fix CL18 type
Tax$Type[ Tax$Type == "N3" ] <- "good"
Tax$Type[ Tax$Type == "Unknown" ] <- "good"

Tax <- Tax[ ,c("ID","Taxonomy","Type", "Block", "MonoID") ]
Tax$Block
levels(Tax$Block) <- c("N1","N2","N3","contaminant","P1","P2","P3","I1","I2","I3")
Tax$Block <- factor(Tax$Block, levels = c("P1","P2","P3",
                                          "I1","I2","I3",
                                          "N1","N2","N3",
                                          "contaminant"))
Tax$Block
Tax <- droplevels(Tax)
levels(Tax$Type)
levels(Tax$Type) <- c("contaminant","Negative","Positive","Indifferent")
Tax$Type <- factor(Tax$Type,
                   levels = c("Indifferent", "Positive",
                              "Negative", "contaminant"))

Tax.colonization <- Tax
devtools::use_data(Tax.colonization,
                   pkg = "~/rhizogenomics/github/wheelP/",
                   overwrite = TRUE)
write.table(Tax.colonization,
            file = "Taxonomy.txt",sep = "\t", col.names = NA, row.names = TRUE,
            quote = FALSE)


