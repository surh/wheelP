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

library(edgeR)
library(wheelP)
# library(AMOR)
# library(GO.db)
# library(org.At.tair.db)
# library(metacoder)

# setwd("~/rhizogenomics/experiments/2017/today7//")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")

date()

# Load previously fit model
data(m1.wheel)
m1 <- m1.wheel
rm(m1.wheel)

# Test data
# load("../today5/m1.test.rdat")
# m1 <- m1.wheel
# rm(m1.wheel)

# Compare bacteria vs no baceria
contrast <- c(0,1,1,1,1,1,1,1,1,1,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0)
qlf <- glmQLFTest(m1,contrast = contrast)
qlf.tt <- topTags(qlf,n = length(row.names(qlf$coefficients)))
qlf.tt <- qlf.tt$table
qlf.tt$Gene <- row.names(qlf.tt)
row.names(qlf.tt) <- NULL
hist(qlf.tt$PValue)
write.table(qlf.tt, "bacteria_vs_nobac.txt", sep ="\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

# gene.set <- row.names(subset(qlf.tt$table, FDR < 0.01 & logFC > 0 ))
# p1 <- plot_set(Dat = Dat.norm,
#                gene.set = gene.set,
#                groups = c("Bacteria", "Phosphate"))
# p1
# ggplot(p1$data,aes(x = Phosphate, y = z.score, color = Bacteria)) +
#   geom_boxplot()


# find genes DE in all syncom
RES <- NULL
for(i in 2:10){
  # i <- 2
  contrast <- c(0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,0,0,
                0)
  contrast[i] <- 1
  qlf <- glmQLFTest(m1,contrast = contrast)
  qlf.tt <- topTags(qlf,n = length(row.names(qlf$coefficients)))
  qlf.tt <- qlf.tt$table
  qlf.tt$Gene <- row.names(qlf.tt)
  row.names(qlf.tt) <- NULL
  qlf.tt$Coef <- colnames(coef(m1))[contrast == 1]
  RES <- rbind(RES,qlf.tt)
}
write.table(RES, "block_main_effects.txt", sep ="\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

RES <- read.table("block_main_effects.txt", header = TRUE)
# head(RES)
# RES2 <- NULL
# for(gene in unique(RES$Gene)){
#   # gene <- unique(RES$Gene)[1]
#   dat <- subset(RES, Gene == gene)
#
#   if(all(dat$logFC < 0 | dat$logFC > 0) & all(dat$PValue < 0.01)){
#     res <- data.frame(Gene = gene, logFC = mean(dat$logFC))
#     RES2 <- rbind(RES2,res)
#   }
# }
# RES2

rm(RES,RES2,dat,gene,res,i)

# Compare positive vs negative
contrast <- c(0,1,1,1,0,0,0,-1,-1,-1,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0)
qlf <- glmQLFTest(m1,contrast = contrast)
qlf.tt <- topTags(qlf,n = length(row.names(qlf$coefficients)))
qlf.tt <- qlf.tt$table
qlf.tt$Gene <- row.names(qlf.tt)
row.names(qlf.tt) <- NULL
hist(qlf.tt$PValue)
write.table(qlf.tt, "positive_vs_negative.txt", sep ="\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)


# Compare negatives in both pre treatment
contrast <- c(0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,1,-1,0,1,-1,0,1,
              -1)
qlf <- glmQLFTest(m1,contrast = contrast)
topTags(qlf,n = 20)
qlf.tt <- topTags(qlf,n = length(row.names(qlf$coefficients)))
qlf.tt <- qlf.tt$table
qlf.tt$Gene <- row.names(qlf.tt)
row.names(qlf.tt) <- NULL
hist(qlf.tt$PValue)
write.table(qlf.tt, "negativeplusp_vs_negativeminusP.txt", sep ="\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

# Compare N1 vs N3 in 30uM
contrast <- c(0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,1,1,0,0,0,0,-1,
              -1)
qlf <- glmQLFTest(m1,contrast = contrast)
topTags(qlf,n = 20)
qlf.tt <- topTags(qlf,n = length(row.names(qlf$coefficients)))
qlf.tt <- qlf.tt$table
qlf.tt$Gene <- row.names(qlf.tt)
row.names(qlf.tt) <- NULL
hist(qlf.tt$PValue)
write.table(qlf.tt, "N1in30uM_vs_N3in30uM.txt", sep ="\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

# Compare N2 vs N3 in 30uM
contrast <- c(0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,1,1,0,-1,
              -1)
qlf <- glmQLFTest(m1,contrast = contrast)
topTags(qlf,n = 20)
qlf.tt <- topTags(qlf,n = length(row.names(qlf$coefficients)))
qlf.tt <- qlf.tt$table
qlf.tt$Gene <- row.names(qlf.tt)
row.names(qlf.tt) <- NULL
hist(qlf.tt$PValue)
write.table(qlf.tt, "N2in30uM_vs_N3in30uM.txt", sep ="\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

# Compare N1 vs N2 in 30uM
contrast <- c(0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,1,1,0,-1,-1,0,0,
              0)
qlf <- glmQLFTest(m1,contrast = contrast)
topTags(qlf,n = 20)
qlf.tt <- topTags(qlf,n = length(row.names(qlf$coefficients)))
qlf.tt <- qlf.tt$table
qlf.tt$Gene <- row.names(qlf.tt)
row.names(qlf.tt) <- NULL
hist(qlf.tt$PValue)
write.table(qlf.tt, "N1in30uM_vs_N2in30uM.txt", sep ="\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

# Compare I3 in plusP_100uM to I2 in plusP_100uM
contrast <- c(0,0,0,0,0,1,-1,0,0,0,
              0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,1,0,0,-1,
              0,0,0,0,0,0,0,0,0,0,
              0)
qlf <- glmQLFTest(m1,contrast = contrast)
topTags(qlf,n = 20)
qlf.tt <- topTags(qlf,n = length(row.names(qlf$coefficients)))
qlf.tt <- qlf.tt$table
qlf.tt$Gene <- row.names(qlf.tt)
row.names(qlf.tt) <- NULL
hist(qlf.tt$PValue)
write.table(qlf.tt, "I2inMP100uM_vs_I3inMP100uM.txt", sep ="\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)


date()
