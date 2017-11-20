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
library(Biostrings)

date()

# Feature data
data(Features)
head(Features)

# Taxonomy
Tax <- read.table('~/rhizogenomics/experiments/2017/today10/growth_curve_ref_taxonomy_RDP.txt',
                  stringsAsFactors = FALSE, comment.char = "", quote = "", sep = ';')
head(Tax, 20)
Tax$V2 <- Tax$V4 <- Tax$V6 <- Tax$V8 <- Tax$V10 <- Tax$V12 <- Tax$V14 <- NULL
row.names(Tax) <- as.character(Tax$V1)
Tax <- Tax[ row.names(Features), ]
colnames(Tax) <- c('Strain.tax', 'Kingdom', 'Phylum','Class','Order','Family','Genus')
head(Tax, 20)

# Ref seqs
seqs <- readDNAStringSet("~/rhizogenomics/experiments/2017/2017-02-27.growthcurves_phylosig/growth_curves_tree/ref.fasta")
setdiff(names(seqs), Features$Strain)
length(setdiff(Features$Strain,names(seqs)))
seqs <- as.character(seqs)
seqs <- data.frame(Strain.seqs = row.names(Features), Seq = seqs[row.names(Features)], row.names = row.names(Features))
head(seqs)
sum(is.na(seqs$Seq))

# Clusters
Clus <- read.table("~/rhizogenomics/experiments/2017/today10/clusters.txt", header = TRUE)
head(Clus)
colnames(Clus)[1] <- 'Strain.clus'
row.names(Clus) <- as.character(Clus$Strain.clus)
Clus <- Clus[ row.names(Features), ]
head(Clus)

## Combine
Dat <- cbind(Features, Clus, Tax, seqs)
head(Dat)
# write.table(Dat, '~/rhizogenomics/experiments/2017/today10/supp_table_growth_curves.txt', sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Dat, 'supp_table_growth_curves.txt', sep = "\t", quote = FALSE, row.names = FALSE)
date()
