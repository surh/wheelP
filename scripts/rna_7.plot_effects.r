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

date()

indir <- "~/rhizogenomics/experiments/2017/2017-03-07.wheel_rna/killdevil/"

data(wheelP.rna)
Dat <- wheelP.rna
rm(wheelP.rna)

# Calculate RPKM
gene.lengths <- read.table("~/rhizogenomics/data/phosphate/gene_lengths.txt",
                           row.names = 1, header = TRUE)
core.pi <- read.table("~/rhizogenomics/data/phosphate/core_pi_up.txt")
Dat.norm <- create_dataset(Tab = rpkm(x = Dat$Tab,
                                      gene.length = gene.lengths[ taxa(Dat), ]),
                           Map = Dat$Map,
                           Tax = Dat$Tax)


filename <- "positive_vs_negative.txt"
filename <- paste(indir,"/",filename,sep = "")
Res <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
head(Res)

Res <- subset(Res, Gene %in% core.pi$V1 & FDR < 0.05)
Res


filename <- "block_main_effects.txt"
filename <- paste(indir,"/",filename,sep = "")
Res <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
head(Res)

gene.set <- core.pi$V1
gene.set <- Res$Gene[ Res$Coef == "P1" & Res$FDR < 0.05 & Res$logFC > 0 ]

Res <- subset(Res, Gene %in% gene.set)

Res$Block <- factor(Res$Coef, levels = c("P1","P2","P3",
                                         "I1","I2","I3",
                                         "N1","N2","N3"))
Res$Increase <- Res$logFC > 0
ftable(Increase ~ Block, Res)


p1 <- ggplot(Res, aes(x = Block, y = logFC)) +
  # geom_boxplot() +
  geom_violin()
p1
