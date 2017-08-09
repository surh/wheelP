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

setwd("~/rhizogenomics/experiments/2017/today4/")

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

#########################################
# filename <- "positive_vs_negative.txt"
# filename <- paste(indir,"/",filename,sep = "")
# Res <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
# head(Res)
#
# Res <- subset(Res, Gene %in% core.pi$V1 & FDR < 0.05)
# Res

#####################################3
# Plot effects of blocks
outdir <- "rna/"
dir.create("rna/")

filename <- "block_main_effects.txt"
filename <- paste(indir,"/",filename,sep = "")
Res <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
head(Res)

# Get matrix of effects
dat <- acast(Gene ~ Coef, data = Res, fun.aggregate = sum,value.var = "logFC")
fdr <- acast(Gene ~ Coef, data = Res, fun.aggregate = sum,value.var = "FDR")

# Turn non-significant to zero and bound coefficients for easier visualization.
dat[ fdr >= 0.05 ] <- 0
dat[ dat > 3 ] <- 3
dat[ dat < -3 ] <- -3

# Remove never significant genes
dat <- dat[ rowSums((dat == 0)) < ncol(dat), ]

# Reorder and count significant
dat <- dat[ ,c("P1","P2","P3","I1","I2","I3","N1","N2","N3") ]
colSums(dat > 0)
colSums(dat < 0)
colSums(dat != 0)

nrow(dat)

filename <- paste(outdir,"/","heatmap_significant_block_effects.png",sep = "")
png(filename,width = 2000,height = 2000, res = 600)
gplots::heatmap.2(dat,scale = "none",trace = "none",
                  col = colorRampPalette(colors = c("red","white","blue"))(11),
                  labRow = "", distfun = dist,dendrogram = "row",Colv = FALSE,
                  #RowSideColors = c("white","black")[ row.names(dat) %in% as.character(core.pi$V1) + 1 ],
                  margins = c(5,5))
dev.off()
#heatmap(dat)

ndiffs <- function(x){
  x.t <- t(x)
  res <- apply(x,1,function(v) colSums(v != x.t))
  return(res)
}

nident <- function(x){
  x.t <- t(x)
  res <- apply(x,1,function(v) colSums(v == x.t))
  return(res)
}

# mat <- matrix(c(1,1,3,2,1,1,4,2,0),ncol = 3)
# ndiffs(mat)
# nident(mat)

numup <- function(x){
  x[ x <= 0] <- 0
  x[ x > 0 ] <- 1
  x.t <- t(x)
  n <- rowSums(x)
  apply(x,1,function(v) colSums(v + x.t == 2) / n)
}

numdown <- function(x){
  x <- 1*(x < 0)
  x.t <- t(x)
  n <- rowSums(x)
  apply(x,1,function(v) colSums(v + x.t == 2) / n)
}

ppv.up <- 100*numup(t(dat))
ppv.dn <- 100*numdown(t(dat))
ppv.up <- melt(ppv.up,varnames = c("Reference.Block","Validation.Block"),value.name = "PPV")
ppv.dn <- melt(ppv.dn,varnames = c("Reference.Block","Validation.Block"),value.name = "PPV")

p1 <- ggplot(ppv.up,aes(x = Reference.Block,y = Validation.Block)) +
  geom_tile(aes(fill = log2(PPV))) +
  # geom_tile(aes(fill = PPV)) +
  scale_fill_gradient2(low = "#d01c8b",mid = "white", high = "#4dac26",midpoint = 4,na.value = "#404040")
p1

p1 <- ggplot(ppv.dn,aes(x = Reference.Block,y = Validation.Block)) +
  geom_tile(aes(fill = log2(PPV))) +
  # geom_tile(aes(fill = PPV)) +
  scale_fill_gradient2(low = "#d01c8b",mid = "white", high = "#4dac26",midpoint = 4,na.value = "#404040")
p1

# Write gene lists
for( block in colnames(dat)){
  set <- dat[,block]

  up <- names(set)[ which(set > 0) ]
  filename <- paste(outdir,"/",block,"_up.txt",sep = "")
  write.table(up,file = filename,quote = FALSE, row.names = FALSE, col.names = FALSE)

  dn <- names(set)[ which(set < 0) ]
  filename <- paste(outdir,"/",block,"_dn.txt",sep = "")
  write.table(dn,file = filename,quote = FALSE, row.names = FALSE, col.names = FALSE)
}

rm(dat,filename,ppv.up,ppv.dn,outdir,block,set)
##################################################3
# plot go
blocks <- c("P1","P2","P3","I1","I2","I3","N1","N2","N3")
godir <- "rna/block_go/"
GO <- NULL
for(block in blocks){
  # block <- blocks[1]

  filename <- paste(godir,"/",block,".up.txt",sep = "")
  cat(filename,"\n")
  up <- read.table(filename, sep = "\t",stringsAsFactors = FALSE, quote = "")

  filename <- paste(godir,"/",block,".dn.txt",sep = "")
  cat(filename,"\n")
  dn <- read.table(filename, sep = "\t",stringsAsFactors = FALSE, quote = "")

  up$Block <- block
  up$Direction <- "up"
  dn$Block <- block
  dn$Direction <- "down"

  # up <- subset(up,V6 < 0.05)
  # dn <- subset(dn,V6 < 0.05)

  GO <- rbind(GO,up,dn)
}
names(GO) <- c("GO.name", "GO.term.size", "GO.ID","GO.found","p.value","FDR","genes","Block","Direction")
head(GO)

GO$Block <- factor(GO$Block, levels = blocks)

GO$log10.FDR <- GO$FDR
GO$log10.FDR[ GO$log10.FDR > 0.05 ] <- 0
GO$log10.FDR <- -log10(GO$log10.FDR)
GO$log10.FDR[ GO$log10.FDR == Inf ] <- NA

dat <- acast(GO.name ~ Direction + Block, fun.aggregate = function(x){sum (x < 0.05)},value.var = "FDR",data = GO)
dat <- dat[ order(rowSums(dat), decreasing = TRUE), ]
head(dat,20)

dat <- dat[ order(rowSums(dat[,1:9]), decreasing = TRUE), ]
head(dat[,1:9],20)

dat <- dat[ order(rowSums(dat[,10:18]), decreasing = TRUE), ]
head(dat[,10:18],20)

gos <- c("DEFENSE_RESPONSE","INNATE_IMMUNE_RESPONSE","IMMUNE_RESPONSE","RESPONSE_TO_OTHER_ORGANISM","RESPONSE_TO_SALICYLIC_ACID_STIMULUS","SYSTEMIC_ACQUIRED_RESISTANCE","RESPONSE_TO_JASMONIC_ACID_STIMULUS","CELLULAR_RESPONSE_TO_PHOSPHATE_STARVATION",
         "REGULATION_OF_PHOSPHATE_METABOLIC_PROCESS","RESPONSE_TO_CHEMICAL_STIMULUS","GLUCOSINOLATE_BIOSYNTHETIC_PROCESS","GLUCOSINOLATE_METABOLIC_PROCESS",
         "ABSCISIC_ACID_METABOLIC_PROCESS","ABSCISIC_ACID_MEDIATED_SIGNALING_PATHWAY","CELLULAR_RESPONSE_TO_ABSCISIC_ACID_STIMULUS","RESPONSE_TO_ABSCISIC_ACID_STIMULUS")

gos <- c("DEFENSE_RESPONSE",
         "SALICYLIC_ACID_BIOSYNTHETIC_PROCESS",
         "JASMONIC_ACID_BIOSYNTHETIC_PROCESS",
         "GLUCOSINOLATE_BIOSYNTHETIC_PROCESS",
         "ROOT_DEVELOPMENT",
         "CELLULAR_RESPONSE_TO_IRON_ION_STARVATION",
         "CELLULAR_RESPONSE_TO_CALCIUM_ION_STARVATION",
         "CELLULAR_RESPONSE_TO_MAGNESIUM_STARVATION",
         "CELLULAR_RESPONSE_TO_NITROGEN_STARVATION",
         "CELLULAR_RESPONSE_TO_PHOSPHATE_STARVATION",
         "CELLULAR_RESPONSE_TO_STARVATION",
         "CELLULAR_RESPONSE_TO_SUCROSE_STARVATION",
         "CELLULAR_RESPONSE_TO_SULFATE_STARVATION",
         "CELLULAR_RESPONSE_TO_SULFUR_STARVATION",
         "RESPONSE_TO_STARVATION")
gos <- c("DEFENSE_RESPONSE",
         "SALICYLIC_ACID_BIOSYNTHETIC_PROCESS",
         "JASMONIC_ACID_BIOSYNTHETIC_PROCESS",
         "GLUCOSINOLATE_BIOSYNTHETIC_PROCESS",
         "ROOT_DEVELOPMENT",
         "CELLULAR_RESPONSE_TO_IRON_ION_STARVATION",
         "CELLULAR_RESPONSE_TO_NITROGEN_STARVATION",
         "CELLULAR_RESPONSE_TO_PHOSPHATE_STARVATION")

dat <- subset(GO,GO.name %in% gos)
dat$GO.name <- factor(dat$GO.name, levels = gos)
dat$Direction <- factor(dat$Direction, levels = c("down","up"))

df_grid <- expand.grid(Direction = levels(dat$Direction), Block = levels(dat$Block), GO.name = levels(dat$GO.name))
dat <- merge(dat,df_grid,by = c("Direction","Block","GO.name"), all = TRUE)

# subset(dat,Block == "N1")
# head(dat)

# dat$DirBlock <- interaction(dat$Block,dat$Direction,drop = FALSE, sep = ".")
p1 <- ggplot(dat, aes(y = GO.name, x = Block )) +
  facet_grid(~ Direction, scales = "free_x", drop = FALSE) +
  geom_tile(aes(fill = log10.FDR), col = "black", size = 2) +
  scale_fill_gradient2(low = "yellow",mid = "orange",
                       high = "red",midpoint = 60,na.value = "white",
                       breaks = c(10,30,60,90,120)) +
  guides(fill = guide_legend(title = "-log10(FDR)")) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, color = "black", face = "bold"),
        axis.text.y = element_text(color = "black", face = "bold"))
p1
ggsave("rna/go_table.svg",p1,width = 10,height = 4)

# gene.set <- core.pi$V1
# gene.set <- Res$Gene[ Res$Coef == "P1" & Res$FDR < 0.05 & Res$logFC > 0 ]
#
#
# Res <- subset(Res, Gene %in% gene.set)
#
# Res$Block <- factor(Res$Coef, levels = c("P1","P2","P3",
#                                          "I1","I2","I3",
#                                          "N1","N2","N3"))
# Res$Increase <- Res$logFC > 0
# ftable(Increase ~ Block, Res)
#
#
# p1 <- ggplot(Res, aes(x = Block, y = logFC)) +
#   # geom_boxplot() +
#   geom_violin()
# p1
