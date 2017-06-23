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
setwd("~/rhizogenomics/experiments/2017/today3/")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")

# READ tables
# dir <- "../today3/tables_mapsplit_trim350_0.985/150903/"
dir <- "tables_mapsplit_trim350_0.985/150903/"
Dat1 <- read_mapsplit_dir(dir = dir)

# dir <- "../today3/tables_mapsplit_trim350_0.985/160201/"
dir <- "tables_mapsplit_trim350_0.985/160201/"
Dat2 <- read_mapsplit_dir(dir = dir)

# dir <- "../today3/tables_mapsplit_trim350_0.985/160909/"
dir <- "tables_mapsplit_trim350_0.985/160909/"
Dat3 <- read_mapsplit_dir(dir = dir)

# dir <- "../today3/tables_mapsplit_trim350_0.985/170118/"
dir <- "tables_mapsplit_trim350_0.985/170118/"
Dat4 <- read_mapsplit_dir(dir = dir)

# Combine and import metadata and taxonomy
DAT <- combine_datasets(list(Dat1,Dat2,Dat3,Dat4))
colSums(DAT$Tab)["ABSENT"]
DAT <- remove_samples(DAT,"ABSENT")
data(Map.colonization)
data(Tax.colonization)

# Clean map
Map.colonization <- subset(Map.colonization, rnaID != "ABSENT")
# Give priority to run 4 over run 2
Map.colonization <- subset(Map.colonization,
                           !(rnaID %in% names(which(table(Map.colonization$rnaID) > 1)) &
                               Run == "R160201"))
Map.colonization <- droplevels(Map.colonization)
row.names(Map.colonization) <- as.character(Map.colonization$rnaID)

setdiff(samples(DAT),row.names(Map.colonization))
setdiff(taxa(DAT),row.names(Tax.colonization))



DAT <- create_dataset(Tab = DAT$Tab,
                      Map = droplevels(Map.colonization[ samples(DAT), ]),
                      Tax = droplevels(Tax.colonization)[ taxa(DAT), ])

DAT <- subset(DAT, Fraction != "BLANK",drop = TRUE, clean = TRUE)
DAT <- subset(DAT, Fraction != "ABSENT",drop = TRUE, clean = TRUE)
DAT$Map$Depth <- colSums(DAT$Tab)
DAT


ftable(Fraction ~ Bacteria, DAT$Map)

wheelP.mapsplit <- DAT
devtools::use_data(wheelP.mapsplit,
                   pkg = "~/rhizogenomics/src/trunk/phosphate_code/wheelP/",
                   overwrite = TRUE)
rm(DAT,Dat1,Dat2,Dat3,Dat4,dir,wheelP.full,wheelP.mapsplit)

##################
# data("wheelP.mapsplit")
# Dat <- wheelP.mapsplit
# rm(wheelP.mapsplit)
#
# p1 <- ggplot(Dat$Map, aes(x = Depth, y = Depthsplit)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   scale_y_log10() +
#   scale_x_log10()
# p1
# ggsave("depth_bymapping.svg", p1, width = 4, height = 4)
# ggsave("depth_bymapping.png", p1, width = 4, height = 4)
#
#
# summary(p1$data$Depth - p1$data$Depthsplit)
# p1 <- ggplot(Dat$Map, aes(x = Depth - Depthsplit)) +
#   geom_histogram(bins = 20) +
#   scale_x_log10() +
#   AMOR::theme_blackbox
# p1
# ggsave("depthdiff.svg", p1, width = 4, height = 4)
# ggsave("depthdiff.png", p1, width = 4, height = 4)
#
#
# p1 <- phylogram(collapse_by_taxonomy(Dat,level = 1),
#                 facet = ~ Fraction,nrow.legend = 7) +
#   theme(axis.text.x = element_blank(),
#         strip.text.x = element_text(angle = 90))
# p1
# ggsave("taxatype_mapsplit_phylogram.svg", p1, width = 8, height = 4)
# ggsave("taxatype_mapsplit_phylogram.png", p1, width = 8, height = 4)
#
# ### Depth
# p1 <- plotgg_var(subset(Dat, Fraction %in% c("Agar","Root")),
#                  var.name = "Depth", x = "Fraction", col = "Bacteria") +
#   scale_y_log10()
# p1
# ggsave("depth.svg", p1, width = 8, height = 4)
# ggsave("depth.png", p1, width = 8, height = 4)
#
# p1 <- plotgg_var(subset(Dat, Fraction %in% c("Agar","Root")),
#                  var.name = "Depthsplit", x = "Fraction", col = "Bacteria") +
#   scale_y_log10()
# p1
# ggsave("depthsplit.svg", p1, width = 8, height = 4)
# ggsave("depthsplit.png", p1, width = 8, height = 4)
#
# # Remove contams
# Dat <- remove_taxons(Dat,taxons = taxa(Dat)[ grep("contaminant",taxa(Dat)) ] )
# Dat$Map$Usablesplit <- colSums(Dat$Tab)
# Dat <- clean(Dat)
# p1 <- plotgg_var(Dat,var.name = "Usablesplit", x = "Fraction") +
#   scale_y_log10()
# p1
# p1 <- plotgg_var(subset(Dat, Fraction %in% c("Agar","Root")),
#                  var.name = "Usablesplit", x = "Fraction", col = "Bacteria") +
#   scale_y_log10()
# p1
# ggsave("Usablesplit.svg", p1, width = 8, height = 4)
# ggsave("Usablesplit.png", p1, width = 8, height = 4)
#
#
# p1 <- ggplot(subset(Dat$Map, Fraction %in% c("Agar","Root")),
#              aes(x = Fraction, col = Bacteria, y = 100*Depthsplit/Depth)) +
#   geom_boxplot()
# p1
#
# p1 <- ggplot(subset(Dat$Map, Fraction %in% c("Agar","Root")),
#              aes(x = Fraction, col = Bacteria, y = 100*Usablesplit/Depthsplit)) +
#   geom_boxplot()
# p1
#
#
# ## Make tax blocks
# Dat$Tax$Block <- NA
# Dat$Tax$Block[ Dat$Tax$ID %in% c("CL21","217","CL59","278",
#                                  "CL96","CL28","113","215","160")] <- "P1"
# # Missing CL4
# Dat$Tax$Block[ Dat$Tax$ID %in% c("499","340","137A","137B","385",
#                                  "CL81", "CL71","474","371")] <- "P2"
# # missing 229
# Dat$Tax$Block[ Dat$Tax$ID %in% c("210","28","146","186","74",
#                                  "111A","111B","121")] <- "P3"
# # missing CL17
# Dat$Tax$Block[ Dat$Tax$ID %in% c("77","CL87","122","10","CL155",
#                                  "329","CL149","143")] <- "I1"
# # missing CL151
# # Recovered 30B with mapsplit
# Dat$Tax$Block[ Dat$Tax$ID %in% c("275","283","CL151","30A","30B","316",
#                                  "211","295","131","36")] <- "I2"
# # missing 240
# Dat$Tax$Block[ Dat$Tax$ID %in% c("370","20","317","363","33",
#                                  "168","165")] <- "I3"
# # Missing 350, recovered it with mapsplit
# Dat$Tax$Block[ Dat$Tax$ID %in% c("233","339","375","123","CL89",
#                                  "259","214","72","350")] <- "N1"
# # Missing CL9
# Dat$Tax$Block[ Dat$Tax$ID %in% c("494","345A","345B","224","69",
#                                  "218","CL11","376","CL32") ] <- "N2"
# # Missing 88, 59
# # Recovered 59 with mapsplit
# Dat$Tax$Block[ Dat$Tax$ID %in% c("384","267","13A","13B","1B",
#                                  "CL52","351","59")] <- "N3"
# Dat$Tax$Block <- factor(Dat$Tax$Block,levels = c("P1","P2","P3",
#                                                  "I1","I2","I3",
#                                                  "N1","N2","N3"))
#
# # Thre straines (30B, 350, 59) were
# # recovfered with mapsplit vs map all (0.985 usearch trim 350)
# Dat$Tax[ is.na(Dat$Tax$Block), ]
# table(Dat$Tax$Block)
#
#
#
# ### Heatmap
# p1 <- heatgg(Dat)
# head(p1$data)
# dat <- p1$data
# dat$Block <- Dat$Tax[ as.character(dat$Taxon), "Block" ]
#
#
# p1 <- ggplot(dat,aes(x = Taxon, y = SAMPLEID)) +
#   facet_grid(Bacteria ~ Block, scales = "free") +
#   geom_tile(aes(fill = Abundance)) +
#   scale_fill_continuous(trans = "log2") +
#   theme(strip.text.y = element_text(angle = 0),
#         axis.text.y = element_blank(),
#         axis.text.x = element_text(angle = 90, color = "black"))
# p1
# ggsave("log2_counts_mapsplit_heatmap.png",p1,width = 10,height = 7)
#
# dat2 <- aggregate(Abundance ~ Block + Bacteria, data = p1$data,FUN = sum)
# p1 <- ggplot(dat2,aes(x = Block, y = Bacteria)) +
#   #facet_grid(Bacteria ~ Block, scales = "free") +
#   geom_tile(aes(fill = Abundance)) +
#   scale_fill_continuous(trans = "log2",low = "white",high = "darkred") +
#   theme(strip.text.y = element_text(angle = 0))
# p1
# ggsave("log2_counts_mapsplit_byblock.png", p1, width = 6, height = 6)
#
# dat <- subset(dat, Fraction == "Inoculum")
#
# p1 <- ggplot(dat,aes(x = Taxon, y = SAMPLEID)) +
#   facet_grid(Bacteria ~ Block, scales = "free") +
#   geom_tile(aes(fill = Abundance)) +
#   scale_fill_continuous(trans = "log2") +
#   theme(strip.text.y = element_text(angle = 0),
#         axis.text.y = element_blank(),
#         axis.text.x = element_text(angle = 90, color = "black"))
# p1
# ggsave("log2_counts_mapsplit_inoculum_heatmap.png",p1,width = 10,height = 7)
#
#
#
# p1 <- heatgg(normalize(Dat,norm = "Usablesplit"))
# head(p1$data)
# dat <- p1$data
# dat$Block <- Dat$Tax[ as.character(dat$Taxon), "Block" ]
#
#
# p1 <- ggplot(dat,aes(x = Taxon, y = SAMPLEID)) +
#   facet_grid(Bacteria ~ Block, scales = "free") +
#   geom_tile(aes(fill = Abundance)) +
#   scale_fill_continuous(trans = "log2") +
#   theme(strip.text.y = element_text(angle = 0),
#         axis.text.y = element_blank(),
#         axis.text.x = element_text(angle = 90, color = "black"))
# p1
# ggsave("log2_percent_mapsplit_heatmap.png",p1,width = 10,height = 7)
#
# dat2 <- aggregate(Abundance ~ Block + Bacteria, data = p1$data,FUN = mean)
# p1 <- ggplot(dat2,aes(x = Block, y = Bacteria)) +
#   #facet_grid(Bacteria ~ Block, scales = "free") +
#   geom_tile(aes(fill = Abundance)) +
#   scale_fill_continuous(trans = "log2",low = "white",high = "darkred") +
#   theme(strip.text.y = element_text(angle = 0))
# p1
# ggsave("log2_percent_mapsplit_byblock.png", p1, width = 6, height = 6)
#
# dat <- subset(dat, Fraction == "Inoculum")
#
# p1 <- ggplot(dat,aes(x = Taxon, y = SAMPLEID)) +
#   facet_grid(Bacteria ~ Block, scales = "free") +
#   geom_tile(aes(fill = Abundance)) +
#   scale_fill_continuous(trans = "log2") +
#   theme(strip.text.y = element_text(angle = 0),
#         axis.text.y = element_blank(),
#         axis.text.x = element_text(angle = 90, color = "black"))
# p1
# ggsave("log2_percent_mapsplit_inoculum_heatmap.png",p1,width = 10,height = 7)
#
#
# ## Compare taxonomy
# compare_mapping_tax <- function(Dat, level, file){
#   # level <- 7
#   # file <- "genus_mapall.txt"
#
#   tax <- pool_samples(collapse_by_taxonomy(Dat,level = level, FUN = sum),
#                       groups = "Bacteria",FUN = sum)
#   tax <- tax$Tab
#   taxall <- read.am(file = file,
#                     format = "qiime",simplify = TRUE)
#   taxall <- taxall[ row.names(tax), colnames(tax) ]
#   colnames(taxall) <- paste(colnames(taxall),".all",sep="")
#   # tax <- as.data.frame(cbind(tax,taxall))
#
#   RES <- NULL
#   for(group in colnames(tax)){
#     #group <- "P3N2"
#     mytaxon1 <- tax[ , group]
#     mytaxon2 <- taxall[,paste(group,".all",sep = "")]
#
#     # barplot(sort(mytaxon2[ mytaxon1 == 0 ]))
#     # sort(mytaxon2[ mytaxon1 == 0 ])
#
#     prop <- mytaxon2[ mytaxon1 == 0]
#     prop <- prop[ which.max(prop) ] / sum(prop)
#     names(prop) <- NULL
#
#     res <- data.frame(group = group,
#                       split.count = sum(mytaxon1),
#                       all.correct = sum(mytaxon2[ mytaxon1 > 0 ]),
#                       all.contam = sum(mytaxon2[ mytaxon1 == 0 ]),
#                       all.all = sum(mytaxon2),
#                       most.common.contam = names(sort(mytaxon2[ mytaxon1 == 0 ],
#                                                       decreasing = T)[1]),
#                       prop.top.cotam = prop)
#     # res
#     RES <- rbind(RES, res)
#   }
#
#   return(RES)
# }
#
#
#
# res <- compare_mapping_tax(Dat = Dat,level = 7,file = "genus_mapall.txt")
# write.table(res,"genus_comparison.csv",sep = "\t", col.names = TRUE,
#             quote = FALSE, row.names = FALSE)
# # Dat$Tax$ID[ grep(pattern = "Acinetobacter", x = as.character(Dat$Tax$Taxonomy)) ]
# # Dat$Tax$ID[ grep(pattern = "Burkholderiaceae; Burkholderia", x = as.character(Dat$Tax$Taxonomy)) ]
# # Dat$Tax$ID[ grep(pattern = "Micrococcaceae; Renibacterium", x = as.character(Dat$Tax$Taxonomy)) ]
# # res$most.common.contam[ grep(pattern = "P2", x = as.character(res$group)) ]
#
# # res <- compare_mapping_tax(Dat = Dat,level = 6,file = "family_mapall.txt")
# # res <- compare_mapping_tax(Dat = Dat,level = 5,file = "order_mapall.txt")
# p1 <- ggplot(res,aes(x = group, y = split.count / all.correct)) +
#   geom_hline(yintercept = 1) +
#   geom_point() +
#   scale_y_log10(breaks = c(0.1, 0.5, 1/.3, 1/1.2, 1/1.1, 1,
#                            1.1, 1.2, 1.3, 2, 10)) +
#   # scale_x_log10() +
#   # geom_abline(intercept = 0, slope = 1)
#   geom_hline(yintercept = 1) +
#   theme(axis.text.x = element_text(angle = 90))
# p1
# ggsave("split.count.vs.all.correct.svg", p1, width = 4, height = 4)
# summary(res$split.count/res$all.correct)
#
# p1 <- ggplot(res,aes(x = group, y = split.count / all.contam)) +
#   geom_hline(yintercept = 1) +
#   geom_point() +
#   scale_y_log10(breaks = c(1/10, 1/8, 1/4, 1/2, 1,
#                            1, 2, 4, 10, 30)) +
#   # scale_x_log10() +
#   # geom_abline(intercept = 0, slope = 1)
#   geom_hline(yintercept = 1) +
#   theme(axis.text.x = element_text(angle = 90))
# p1
# ggsave("split.count.vs.all.contam.svg", p1, width = 4, height = 4)
# summary(res$split.count/res$all.contam)
#
# p1 <- ggplot(res,aes(x = group, y = 100*all.contam / all.all)) +
#   geom_hline(yintercept = 1) +
#   geom_bar(fill = "black", stat = "identity") +
#   scale_y_log10(breaks = 100*c(0.02, 0.1, 0.25, 0.5, 0.75)) +
#   # scale_x_log10() +
#   # geom_abline(intercept = 0, slope = 1)
#   # geom_hline(yintercept = 1) +
#   theme(axis.text.x = element_text(angle = 90))
# p1
# ggsave("perc.all.contam.svg", p1, width = 4, height = 4)
# summary(100*res$all.contam / res$all.all)
#
#
# p1 <- ggplot(res,aes(x = 100*all.contam / all.all,
#                      y = prop.top.cotam)) +
#   geom_point()
# p1
# ggsave("perc.all.contam.vs.top.contam.svg", p1, width = 4, height = 4)
#
#
# #### Percent map
#
# mapsplit <- read.table("percent_map/mapsplit.txt", sep = "\t")
# mapall <- read.table("percent_map/mapall.txt", sep = "\t")
#
# head(mapsplit)
# mapsplit <- cbind(aggregate(V2 ~ V1, data = mapsplit, FUN = sum),
#                   aggregate(V3 ~ V1, data = mapsplit, FUN = sum))
# all(mapsplit[,1] == mapsplit[,3])
# mapsplit <- mapsplit[,-3]
# colnames(mapsplit) <- c("ID","Total","Hit")
# head(mapsplit)
# mapsplit$Percent <- 100 * mapsplit$Hit / mapsplit$Total
# row.names(mapsplit) <- mapsplit$ID
#
# head(mapall)
# mapall <- cbind(aggregate(V2 ~ V1, data = mapall, FUN = sum),
#                 aggregate(V3 ~ V1, data = mapall, FUN = sum))
# all(mapall[,1] == mapall[,3])
# mapall <- mapall[,-3]
# colnames(mapall) <- c("ID","Total","Hit")
# head(mapall)
# mapall$Percent <- 100 * mapall$Hit / mapall$Total
# row.names(mapall) <- mapall$ID
#
#
# Dat$Map$Percent.mapsplit <- mapsplit[ samples(Dat), "Percent" ]
# Dat$Map$Percent.mapall <- mapall[ samples(Dat), "Percent" ]
#
# p1 <- ggplot(Dat$Map, aes(x = Percent.mapall, y = Percent.mapsplit)) +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1)
# p1
#
# summary(Dat$Map$Percent.mapall - Dat$Map$Percent.mapsplit)
# p1 <- ggplot(Dat$Map, aes(x = Percent.mapall - Percent.mapsplit)) +
#   geom_histogram(bins = 20)
# p1
#
# p1 <- plotgg_var(subset(Dat, Fraction %in% c("Agar","Root", "Inoculum")),
#                  var.name = "Percent.mapsplit", x = "Fraction", col = "Bacteria")
# p1
#
# m1 <- aov(Percent.mapsplit ~ Fraction * Bacteria,
#           data = subset(Dat$Map, Fraction %in% c("Agar","Root", "Inoculum")))
# summary(m1)
# TukeyHSD(m1)
#
# m1 <- lm(Percent.mapsplit ~ Fraction + G1 + G2 + G3 + N1 + N2 + N3 + B1 + B2 + B3,
#          data = subset(Dat$Map, Fraction %in% c("Agar","Root")))
# summary(m1)
#
#
# Dat$Map$Bacteria.block1 <- substr(as.character(Dat$Map$Bacteria),1,2)
# Dat$Map$Bacteria.block2 <- substr(as.character(Dat$Map$Bacteria),3,4)
# Dat$Map$Bacteria.block1 <- factor(Dat$Map$Bacteria.block1,
#                                   levels = c("G1","G2","G3","N1","N2","N3","B1","B2","B3"))
# Dat$Map$Bacteria.block2 <- factor(Dat$Map$Bacteria.block2,
#                                   levels = c("G1","G2","G3","N1","N2","N3","B1","B2","B3"))
#
# p1 <- ggplot(subset(Dat$Map, Fraction %in% c("Agar","Root")),
#              aes(x = Fraction, y = Percent.mapsplit)) +
#   facet_grid(Bacteria.block1 ~ Bacteria.block2, drop = FALSE) +
#   #geom_tile(aes(fill = Percent.mapsplit))
#   geom_violin(draw_quantiles = 0.5, scale = "count", adjust = 2) +
#   geom_hline(yintercept = 27.85)
# p1
#
#
#
# p1 <- plotgg_var(subset(Dat$Map, Fraction %in% c("Agar","Root")),
#                  var.name = "Percent.mapsplit", x = "Fraction", col = "Bacteria")
# p1
#
# p1 <- plotgg_var(subset(Dat$Map, Fraction %in% c("Agar","Root")),
#                  var.name = "Percent.mapsplit", x = "Fraction", col = "Experiment")
# p1
#
# p1 <- plotgg_var(subset(Dat$Map, Fraction %in% c("Agar","Root")),
#                  var.name = "Percent.mapsplit", x = "Fraction", col = "Plate")
# p1
#
# p1 <- plotgg_var(subset(Dat$Map, Fraction %in% c("Agar","Root")),
#                  var.name = "Percent.mapsplit", x = "Fraction", col = "Run")
# p1
#
# p1 <- plotgg_var(subset(Dat$Map, Fraction %in% c("Agar","Root")),
#                  var.name = "Percent.mapsplit", x = "Fraction", col = "Barcode2")
# p1
#
# p1 <- plotgg_var(subset(Dat$Map, Fraction %in% c("Agar","Root")),
#                  var.name = "Percent.mapsplit", x = "Fraction", col = "Frameshift")
# p1
#
# p1 <- plotgg_var(subset(Dat$Map, Fraction %in% c("Agar","Root")),
#                  var.name = "Percent.mapsplit", x = "Fraction", col = "Well")
# p1
#
#
# variables(Dat)
#
# # SC0314 G1G2
# p1 <- phylogram(collapse_by_taxonomy(Dat,level = 6))
# head(p1$data)
# dat <- droplevels(subset(p1$data, Sample == "SC0314" & Abundance > 0))
#
# p1 <- ggplot(dat, aes(x = Sample, y = Abundance, fill = Taxon)) +
#   geom_bar(stat = "identity")
# p1
#
#
# write.table(subset(Dat$Tax, Block == "G1" | Block == "G2"),
#             file = "g1g2.table.txt",sep = "\t", quote = FALSE,
#             col.names = TRUE, row.names = FALSE)

