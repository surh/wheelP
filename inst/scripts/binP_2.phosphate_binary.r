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
# setwd("~/rhizogenomics/experiments/2017/today8/")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")
# source("/home/sur/rhizogenomics/src/trunk/phosphate_code/monoP/functions.r")

# Load data
data(binP.all)

# Get strains
Strains <- levels(binP.all$Treatment)
Strains <- Strains[ Strains != "No Bacteria" ]

# test
strain.test <- test_all_strains(Strains = Strains, treatment_col = "Treatment",
                                Master = binP.all, plot = FALSE)

# Plot results
p1 <- plot_res(test = strain.test, qval_thres = 0.05)
p1
ggsave("strain_res_qval0.05.svg",p1,width = 5, height = 30)
p1 <- plot_res(test = strain.test, qval_thres = 0.1)
p1
ggsave("strain_res_qval0.1.svg",p1,width = 5, height = 30)

# Write output
write.table(strain.test$RES,file = "bin_lm.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
write.table(strain.test$SE,file = "bin_lm_se.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
write.table(strain.test$PVAL,file = "bin_lm_pval.txt",sep="\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
write.table(strain.test$QVAL,file = "bin_lm_qval.txt",sep="\t",row.names = FALSE,col.names = TRUE, quote = FALSE)

# Test for enrichment of phosphate modulators per condition
set.seed(2810)
perm_binaryP_enrich(strain.test)
set.seed(2810)
perm_binaryP_enrich(strain.test,min.val = -Inf, max.val = 0)
set.seed(2810)
perm_binaryP_enrich(strain.test,min.val = 0, max.val = Inf)

# Make data.frame of individual strain test results
binP <- data.frame(strain.test$RES,strain.test$SE,strain.test$QVAL)
head(binP)
row.names(binP) <- as.character(binP$Strain)
binP$Strain.1 <- binP$Strain.2 <- NULL
colnames(binP)[6:9] <- paste(colnames(binP[2:5]),".se",sep = "")
colnames(binP)[10:13] <- paste(colnames(binP[2:5]),".qval",sep = "")
head(binP)

# Get groups for
# Good ones are positive in at least one of the plusP pre-treatments
good <- binP$Strain[ (binP$plusP_100uM > 0 & binP$plusP_100uM.qval < 0.1) |
                       (binP$plusP_30uM > 0  & binP$plusP_30uM.qval < 0.1) ]

# Verybad are bad in at least three conditions
verybad <- binP$Strain[ rowSums(binP[,2:5] < 0 & binP[,10:13] < 0.1) >= 3 ]

# Bad are bad in at least two conditions with higher significance
bad <- binP$Strain[ rowSums(binP[,2:5] < 0 & binP[,10:13] < 0.05) >= 2 ]

# Indifferent strains are not significant in any condtions
# Originally picked a random, but when I changed the factor names for homogenizing
# the strain names accross datasets, it seems like the factor order changed and so
# the same seed produces a different result. Since the experiment is already done,
# I am just manually indicating which strains were used.
# set.seed(313)
# med <- sample(x = med,size = 26,replace = FALSE)
med <- c("77","CL87","122","CL17","10","CL155",
         "329","CL149","143","275","283",
         "CL151","30","316","211","295","131",
         "36","370","20","317","240","363","33",
         "168","165")

# Calculate mean
binP$Mean <- rowMeans(binP[,2:5])

# Assign groups
binP$Group <- "none"
binP$Group[ binP$Strain %in% good ] <- "Positive"
binP$Group[ binP$Strain %in% setdiff(union(bad,verybad),c("R219","pf"))] <- "Negative"
binP$Group[ binP$Strain %in% med] <- "Indifferent"
binP$Group[ binP$Mean == min(binP$Mean[ binP$Group == "none" ])] <- "Negative" # add 1 to have same sizes

# Order data
binP$Group <- factor(binP$Group,levels = c("none","Negative","Indifferent","Positive"))
binP <- binP[ order(binP$Group, binP$Mean,decreasing = TRUE ),]
table(binP$Group)
binP$Functional.class <- binP$Group
binP$Group <- NULL
head(binP)

# Add block
binP$Block <- NA
binP$Block[ binP$Functional.class == "Positive" ] <- do.call(c,apply(data.frame(i = c('P1','P2','P3'),
                                                                                t = c(9,9,8)),1,function(x) rep(x[1], x[2])))



binP$Block[ binP$Functional.class == "Indifferent" ] <- do.call(c,apply(data.frame(i = c('I1','I2','I3'),
                                                                              t = c(9,9,8)),1,function(x) rep(x[1], x[2])))


binP$Block[ binP$Functional.class == "Negative" ] <- do.call(c,apply(data.frame(i = c('N1','N2','N3'),
                                                                                   t = c(9,9,8)),1,function(x) rep(x[1], x[2])))
binP$Block <- factor(binP$Block, levels = c('P1','P2','P3','I1','I2','I3','N1','N2','N3'))
ftable(Block ~ Functional.class, data = binP)
head(binP,30)

# Save dataset in package
# confirmed it matches original data
devtools::use_data(binP,pkg = "~/rhizogenomics/src/trunk/phosphate_code/wheelP/",
                   overwrite = TRUE)
write.table(binP,'supp_table_binary.txt', sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(binP,'~/rhizogenomics/experiments/2017/today10/supp_table_binary.txt', sep = "\t", quote = FALSE, row.names = FALSE)



