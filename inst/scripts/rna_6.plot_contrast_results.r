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
library(metacoder)
library(edgeR)
# library(org.At.tair.db)
# library(GO.db)

date()

setwd("~/rhizogenomics/experiments/2017/today9/")
# devtools::document("~/rhizogenomics/src/trunk/phosphate_code/wheelP/")
indir <- "~/rhizogenomics/experiments/2017/2017-03-07.wheel_rna/killdevil/"

data(wheelP.rna)
Dat <- wheelP.rna
rm(wheelP.rna)

# Calculate RPKM
gene.lengths <- read.table("~/rhizogenomics/data/phosphate/gene_lengths.txt",
                           row.names = 1, header = TRUE)
Dat.norm <- create_dataset(Tab = rpkm(x = Dat$Tab,
                                      gene.length = gene.lengths[ taxa(Dat), ]),
                           Map = Dat$Map,
                           Tax = Dat$Tax)

# The annotation file from tair can be downloaded at:
# ftp://ftp.arabidopsis.org/home/tair/Ontologies/Gene_Ontology/
Annot <- read.table("~/rhizogenomics/data/tair/2017-11-15.ATH_GO_GOSLIM.txt",
                    sep = "\t", header = FALSE, quote = '')
head(Annot)

###### Bacteria no bacteria #######
filename <- "bacteria_vs_nobac.txt"
filename <- paste(indir,"/",filename,sep = "")
Res <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
head(Res)

# dat <- Res
# output_folder <- "bacteria_vs_nobac/"
# output_format <- "svg"
# # min_p_value <-0.000001
# type <- GO.db::GOBPPARENTS
# prefix <- "gobp"

data <- metacoder_plot_go(dat = Res,output_folder = "bacteria_vs_nobac/",output_format = "svg",
                          min_fdr = 0.01,type = GO.db::GOBPPARENTS,prefix = "gobp",n.supertaxa = 9,
                          num.changed = 3)

# data <- wheelP::parse_tax_and_plot(file = "bacteria_vs_nobac/gobp_go_res.txt", col = 5,
#                                    output_file = "bacteria_vs_nobac/test.svg",
#                                    n.supertaxa = 9,
#                                    num.changed = 3,
#                                    min_fdr = 0.000001)


#### Try to plot heatmap with GO terms


# GO:0043207  response to external biotic stimulus

# GO:0006952  defense response

Res <- droplevels(subset(Res,FDR < 1e-5))
head(Res)
# Res


gos <- droplevels(subset(Annot,V1 %in% Res$Gene & V8 == 'P'))
gos <- lapply(levels(gos$V6),function(x,gos, Annot){
  d <- unique(as.character(subset(gos,V6 == x)$V1))
  d <- data.frame(GO = x, Count = length(d))
  a <- unique(subset(Annot, V6 == x)$V5)
  d$Annotation <- as.character(a)
  return(d)
}, gos = gos, Annot = Annot)
gos <- do.call(rbind,gos)
gos <- gos[ order(gos$Count, decreasing = TRUE), ]
# gos
head(gos, 20)

# selected_gos <- c('GO:0009751','GO:0009407','GO:0051707')
selected_gos <- c('GO:0009627','GO:0042742','GO:0009697','GO:0015706','GO:0009617','GO:0009862','GO:0010363',
                  'GO:0000165','GO:0009867','GO:0031348','GO:0006952','GO:0009751')

dat <- lapply(selected_gos, function(x){

  genes <- subset(Annot,V6 == x)
  genes <- unique(as.character(genes$V1))

  data.frame(logFC = 10*(Res$Gene %in% genes),
             logCPM = NA, F = NA,
             PValue = NA, FDR = NA,
             Gene = Res$Gene,
             Type = x)
})
dat <- do.call(rbind,dat)
Res$Type <- 'Gene'
dat <- rbind(Res,dat)
dat$Gene <- factor(dat$Gene, levels = as.character(unique(Res$Gene[ order(Res$logFC) ])))
dat$logFC[ dat$logFC > 10 ] <- 10
dat$logFC[ dat$logFC < -10 ] <- -10
head(dat)

p1 <- ggplot(dat,aes(x = Type, y = Gene, fill = logFC)) +
  facet_grid(~ Type, scales = "free_x") +
  geom_tile() +
  scale_fill_gradient2(low =  c("#8e0152","#de77ae"),
                       mid = "#f7f7f7",
                       high = c("#7fbc41","#276419"),
                       midpoint = 0,
                       na.value = "#404040",
                       guide = guide_colorbar(title = "logFC")) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, color = 'black'),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank())
p1




rm(data)
#############  Block main effects ##########
filename <- "block_main_effects.txt"
filename <- paste(indir,"/",filename,sep = "")
Res <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
Res$Coef <- factor(Res$Coef,levels = c("P1","P2","P3",
                                       "I1","I2","I3",
                                       "N1","N2","N3"))
head(Res)
p1 <- ggplot2::ggplot(Res,ggplot2::aes(x = PValue)) +
  ggplot2::facet_wrap(~ Coef, ncol = 3) +
  ggplot2::geom_histogram(bins = 20)
p1
paste(ou)
ggplot2::ggsave("bloc_main_effects_pvals.svg",p1, width = 4, height = 4)

RES2 <- NULL
for(gene in unique(Res$Gene)){
  # gene <- unique(Res$Gene)[1]
  dat <- subset(Res, Gene == gene)

  pos.group <- sum(dat$logFC > 0)
  neg.group <- sum(dat$logFC <=0 )
  major.group <- which.max(c(pos = pos.group, neg = neg.group))

  if(major.group == 2 & all(subset(dat,logFC <= 0)$FDR < 0.05) & neg.group >= 7){
    res <- data.frame(Gene = gene, logFC = mean(dat$logFC[dat$logFC <= 0]),
                      pos = pos.group, neg = neg.group)
    RES2 <- rbind(RES2,res)
  }else if(major.group == 1 & all(subset(dat,logFC > 0)$FDR < 0.05) & pos.group >= 7){
    res <- data.frame(Gene = gene, logFC = mean(dat$logFC[dat$logFC > 0]),
                      pos = pos.group, neg = neg.group)
    RES2 <- rbind(RES2,res)
  }

  # if(all(dat$logFC < 0 | dat$logFC > 0) & all(dat$PValue < 0.01)){
  #   res <- data.frame(Gene = gene, logFC = mean(dat$logFC))
  #   RES2 <- rbind(RES2,res)
  # }
}
RES2

# p1 <- plotgg_taxon(Dat.norm, taxon = RES2$Gene[1],x = "Phosphate", col = "Bacteria")
# p1
# p1 <- plotgg_taxon(Dat.norm, taxon = RES2$Gene[2],x = "Phosphate", col = "Bacteria")
# p1

rm(Res,RES2,res,dat,gene,major.group,neg.group,pos.group,p1)

#############  Block main effects ##########
filename <- "positive_vs_negative.txt"
filename <- paste(indir,"/",filename,sep = "")
Res <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
head(Res)

data <- metacoder_plot_go(dat = Res,output_folder = "positive_vs_negative/",
                          output_format = "svg",
                          min_fdr = 0.01,type = GO.db::GOBPPARENTS,
                          prefix = "gobp",n.supertaxa = 9,
                          num.changed = 3)

################ Compare negatives in both pre treatment #############
filename <- "negativeplusp_vs_negativeminusP.txt"
filename <- paste(indir,"/",filename,sep = "")
Res <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
head(Res)

# data <- metacoder_plot_go(dat = Res,output_folder = "negativeplusP_vs_negativeminusP/",
#                           output_format = "svg",
#                           min_fdr = 0.01,type = GO.db::GOBPPARENTS,
#                           prefix = "gobp",n.supertaxa = 9,
#                           num.changed = 3)


################ N1 vs N3 in 30uM #############
filename <- "N1in30uM_vs_N3in30uM.txt"
filename <- paste(indir,"/",filename,sep = "")
Res <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
head(Res)

data <- metacoder_plot_go(dat = Res,output_folder = "N1in30uM_vs_N3in30uM/",
                          output_format = "svg",
                          min_fdr = 0.01,type = GO.db::GOBPPARENTS,
                          prefix = "gobp",n.supertaxa = 9,
                          num.changed = 3)

################ N2 vs N3 in 30uM #############
filename <- "N2in30uM_vs_N3in30uM.txt"
filename <- paste(indir,"/",filename,sep = "")
Res <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
head(Res)

data <- metacoder_plot_go(dat = Res,output_folder = "N2in30uM_vs_N3in30uM/",
                          output_format = "svg",
                          min_fdr = 0.05,type = GO.db::GOBPPARENTS,
                          prefix = "gobp",n.supertaxa = 9,
                          num.changed = 3)

################ N1 vs N2 in 30uM #############
filename <- "N1in30uM_vs_N2in30uM.txt"
filename <- paste(indir,"/",filename,sep = "")
Res <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
head(Res)

data <- metacoder_plot_go(dat = Res,output_folder = "N1in30uM_vs_N2in30uM/",
                          output_format = "svg",
                          min_fdr = 0.01,type = GO.db::GOBPPARENTS,
                          prefix = "gobp",n.supertaxa = 9,
                          num.changed = 3)

################ (I2 vs I3) in minusP_100uM #############
filename <- "I2inMP100uM_vs_I3inMP100uM.txt"
filename <- paste(indir,"/",filename,sep = "")
Res <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
head(Res)

data <- metacoder_plot_go(dat = Res,output_folder = "I2inMP100uM_vs_I3inMP100uM/",
                          output_format = "svg",
                          min_fdr = 0.1,type = GO.db::GOBPPARENTS,
                          prefix = "gobp",n.supertaxa = 9,
                          num.changed = 3)


####### Plot some genes
# IPS1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G09922",x = "Phosphate",col = "Bacteria")
# PHO2
p2 <- plotgg_taxon(Dat.norm,taxon = "AT2G33770",x = "Phosphate",col = "Bacteria")

p3 <- ggplot(data.frame(Phosphate = p1$data$Phosphate,
                        Bacteria = p1$data$Bacteria,
                        IPS1 = p1$data$Abundance,
                        PHO2 = p2$data$Abundance),
             aes(x = IPS1, y = PHO2)) +
  #facet_wrap(~ Phosphate, ncol = 2) +
  geom_point(aes(col = Bacteria)) +
  scale_y_log10(breaks = c(1,10,15,30,50,100)) +
  scale_x_log10() +
  theme_blackbox
p3
ggsave("IPS1_vs_PHO2.svg",p3,width = 4, height = 4)

# IPS1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G09922",x = "Phosphate",col = "Bacteria")
# PHO1
p2 <- plotgg_taxon(Dat.norm,taxon = "AT3G23430",x = "Phosphate",col = "Bacteria")

p3 <- ggplot(data.frame(Phosphate = p1$data$Phosphate,
                        Bacteria = p1$data$Bacteria,
                        IPS1 = p1$data$Abundance,
                        PHO1 = p2$data$Abundance),
             aes(x = IPS1, y = PHO1)) +
  #facet_wrap(~ Phosphate, ncol = 2) +
  geom_point(aes(col = Bacteria)) +
  # scale_y_log10() +
  scale_x_log10() +
  theme_blackbox
p3

# IPS1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G09922",x = "Phosphate",col = "Bacteria")
# NLA
p2 <- plotgg_taxon(Dat.norm,taxon = "AT1G02860",x = "Phosphate",col = "Bacteria")

p3 <- ggplot(data.frame(Phosphate = p1$data$Phosphate,
                        Bacteria = p1$data$Bacteria,
                        IPS1 = p1$data$Abundance,
                        NLA = p2$data$Abundance),
             aes(x = IPS1, y = NLA)) +
  geom_point(aes(col = Bacteria)) +
  scale_y_log10(breaks = c(5,10,20,30)) +
  scale_x_log10() +
  theme_blackbox
p3
ggsave("IPS1_vs_NLA.svg",p3,width = 4, height = 4)

# IPS1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G09922",x = "Phosphate",col = "Bacteria")
# MYC2
p2 <- plotgg_taxon(Dat.norm,taxon = "AT1G32640",x = "Phosphate",col = "Bacteria")

p3 <- ggplot(data.frame(Phosphate = p1$data$Phosphate,
                        Bacteria = p1$data$Bacteria,
                        IPS1 = p1$data$Abundance,
                        MYC2 = p2$data$Abundance),
             aes(x = IPS1, y = MYC2)) +
  geom_point(aes(col = Bacteria)) +
  scale_y_log10(breaks = c(10,20,50,100,200)) +
  scale_x_log10() +
  theme_blackbox
p3

# IPS1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G09922",x = "Phosphate",col = "Bacteria")
# SUR1
p2 <- plotgg_taxon(Dat.norm,taxon = "AT2G20610",x = "Phosphate",col = "Bacteria")

p3 <- ggplot(data.frame(Phosphate = p1$data$Phosphate,
                        Bacteria = p1$data$Bacteria,
                        IPS1 = p1$data$Abundance,
                        SUR1 = p2$data$Abundance),
             aes(x = IPS1, y = SUR1)) +
  geom_point(aes(col = Bacteria)) +
  scale_y_log10() +
  scale_x_log10() +
  theme_blackbox
p3
ggsave("IPS1_vs_SUR1.svg",p3,width = 4, height = 4)

# IPS1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G09922",x = "Phosphate",col = "Bacteria")
# VSP2
p2 <- plotgg_taxon(Dat.norm,taxon = "AT5G24770",x = "Phosphate",col = "Bacteria")

p3 <- ggplot(data.frame(Phosphate = p1$data$Phosphate,
                        Bacteria = p1$data$Bacteria,
                        IPS1 = p1$data$Abundance,
                        VSP2 = p2$data$Abundance),
             aes(x = IPS1, y = VSP2)) +
  geom_point(aes(col = Bacteria)) +
  #facet_wrap(~ Phosphate, ncol = 2) +
  scale_y_log10() +
  scale_x_log10() +
  theme_blackbox
p3
ggsave("IPS1_vs_VSP2.svg",p3,width = 4, height = 4)


# IPS1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G09922",x = "Phosphate",col = "Bacteria")
# WRKY54
p2 <- plotgg_taxon(Dat.norm,taxon = "AT2G40750",x = "Phosphate",col = "Bacteria")

p3 <- ggplot(data.frame(Phosphate = p1$data$Phosphate,
                        Bacteria = p1$data$Bacteria,
                        IPS1 = p1$data$Abundance,
                        WRKY54 = p2$data$Abundance),
             aes(x = IPS1, y = WRKY54)) +
  geom_point(aes(col = Bacteria)) +
  scale_y_log10() +
  scale_x_log10() +
  theme_blackbox
p3

# SUR1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT2G20610",x = "Phosphate",col = "Bacteria")
# ARGOS
p2 <- plotgg_taxon(Dat.norm,taxon = "AT3G59900",x = "Phosphate",col = "Bacteria") +
  scale_y_log10()
p2

p3 <- ggplot(data.frame(Phosphate = p1$data$Phosphate,
                        Bacteria = p1$data$Bacteria,
                        SUR1 = p1$data$Abundance,
                        ARGOS = p2$data$Abundance),
             aes(x = SUR1, y = ARGOS)) +
  #facet_wrap(~ Phosphate, ncol = 2) +
  geom_point(aes(col = Bacteria)) +
  #geom_smooth(method = "lm") +
  scale_y_log10() +
  scale_x_log10() +
  theme_blackbox
p3
ggsave("SUR1_vs_ARGOS.svg",p3,width = 4, height = 4)

# IPS1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G09922",x = "Phosphate",col = "Bacteria")
# ARGOS
p2 <- plotgg_taxon(Dat.norm,taxon = "AT3G59900",x = "Phosphate",col = "Bacteria")

p3 <- ggplot(data.frame(Phosphate = p1$data$Phosphate,
                        Bacteria = p1$data$Bacteria,
                        IPS1 = p1$data$Abundance,
                        ARGOS = p2$data$Abundance),
             aes(x = IPS1, y = ARGOS)) +
  #facet_wrap(~ Phosphate, ncol = 2) +
  geom_point(aes(col = Bacteria)) +
  scale_y_log10() +
  scale_x_log10() +
  theme_blackbox
p3
ggsave("IPS1_vs_ARGOS.svg",p3,width = 4, height = 4)

# PHO2
p1 <- plotgg_taxon(Dat.norm,taxon = "AT2G33770",x = "Phosphate",col = "Bacteria")
# IPS1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G09922",x = "Phosphate",col = "Bacteria")
# PHO1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G23430",x = "Phosphate",col = "Bacteria")
# NLA
p1 <- plotgg_taxon(Dat.norm,taxon = "AT1G02860",x = "Phosphate",col = "Bacteria")
# PHT1;4
p1 <- plotgg_taxon(Dat.norm,taxon = "AT2G38940",x = "Phosphate",col = "Bacteria")
# PAD4
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G52430",x = "Phosphate",col = "Bacteria")
# NPR1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT1G64280",x = "Phosphate",col = "Bacteria")
# MYC2
p1 <- plotgg_taxon(Dat.norm,taxon = "AT1G32640",x = "Phosphate",col = "Bacteria")
# SUR1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT2G20610",x = "Phosphate",col = "Bacteria")
# WRKY54
p1 <- plotgg_taxon(Dat.norm,taxon = "AT2G40750",x = "Phosphate",col = "Bacteria")
# ARGOS
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G59900",x = "Phosphate",col = "Bacteria")
# VSP2
p1 <- plotgg_taxon(Dat.norm,taxon = "AT5G24770",x = "Phosphate",col = "Bacteria")
# DWF1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G19820",x = "Phosphate",col = "Bacteria")
# CBB3
p1 <- plotgg_taxon(Dat.norm,taxon = "AT5G05690",x = "Phosphate",col = "Bacteria")
# ARL
p1 <- plotgg_taxon(Dat.norm,taxon = "AT2G44080",x = "Phosphate",col = "Bacteria")
# EIN3
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G20770",x = "Phosphate",col = "Bacteria")
# XPL1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G18000",x = "Phosphate",col = "Bacteria")
# PMEI4
p1 <- plotgg_taxon(Dat.norm,taxon = "AT4G25250",x = "Phosphate",col = "Bacteria")
# EIN2
p1 <- plotgg_taxon(Dat.norm,taxon = "AT5G03280",x = "Phosphate",col = "Bacteria")
# BAK1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT4G33430",x = "Phosphate",col = "Bacteria")
p1
# SERK1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT1G71830",x = "Phosphate",col = "Bacteria")
p1
# SERK2
p1 <- plotgg_taxon(Dat.norm,taxon = "AT1G34210",x = "Phosphate",col = "Bacteria")
p1
# SERK4
p1 <- plotgg_taxon(Dat.norm,taxon = "AT2G13790",x = "Phosphate",col = "Bacteria")
p1
# SERK5
p1 <- plotgg_taxon(Dat.norm,taxon = "AT2G13800",x = "Phosphate",col = "Bacteria")
p1
# GAI
p1 <- plotgg_taxon(Dat.norm,taxon = "AT1G14920",x = "Phosphate",col = "Bacteria")
p1
# RGL
p1 <- plotgg_taxon(Dat.norm,taxon = "AT1G66350",x = "Phosphate",col = "Bacteria")
p1
# GID1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G05120",x = "Phosphate",col = "Bacteria")
p1
# RGA
p1 <- plotgg_taxon(Dat.norm,taxon = "AT2G01570",x = "Phosphate",col = "Bacteria")
p1
# SCR
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G54220",x = "Phosphate",col = "Bacteria")
p1
# AUX1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT2G38120",x = "Phosphate",col = "Bacteria")
p1

# LSD1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT1G62830",x = "Phosphate",col = "Bacteria")
p1


# Plot I2I3 & N2N3 plusP_30uM
# ARGOS
p1 <- plotgg_taxon(subset(Dat.norm, Phosphate == "plusP_30uM"),
                   taxon = "AT3G59900",x = "Phosphate",col = "Bacteria") +
  ggtitle(label = "ARGOS (AT3G59900)")
p1
ggsave("ARGOS_plusP_30uM.svg",p1,width = 5, height = 4)
# SERK2
p1 <- plotgg_taxon(subset(Dat.norm, Phosphate == "plusP_30uM"),
                   taxon = "AT1G34210",x = "Phosphate",col = "Bacteria") +
  ggtitle(label = "SERK2 (AT1G34210)")
p1
ggsave("SERK2_plusP_30uM.svg",p1,width = 5, height = 4)
#
p1 <- plotgg_taxon(subset(Dat.norm, Phosphate == "plusP_30uM"),
                   taxon = "AT5G48380",x = "Phosphate",col = "Bacteria") +
  ggtitle(label = " ()")
p1

# IPS1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G09922",x = "Phosphate",col = "Bacteria")
# ARGOS
p2 <- plotgg_taxon(Dat.norm,taxon = "AT3G59900",x = "Phosphate",col = "Bacteria")

temp <- data.frame(Phosphate = p1$data$Phosphate,
                   Bacteria = p1$data$Bacteria,
                   IPS1 = p1$data$Abundance,
                   ARGOS = p2$data$Abundance)
temp <- subset(temp, Phosphate == "minusP_30uM")

p3 <- ggplot(temp,
             aes(x = IPS1, y = ARGOS)) +
  #facet_wrap(~ Phosphate, ncol = 2) +
  geom_point(aes(col = Bacteria)) +
  geom_point(data = subset(temp, Bacteria %in% c("I2I3","N2N3")),
             shape = 1, col = "black", size = 2) +
  scale_y_log10() +
  scale_x_log10() +
  theme_blackbox
p3

# IPS1
p1 <- plotgg_taxon(Dat.norm,taxon = "AT3G09922",x = "Phosphate",col = "Bacteria")
# ARL
p2 <- plotgg_taxon(Dat.norm,taxon = "AT2G44080",x = "Phosphate",col = "Bacteria")

temp <- data.frame(Phosphate = p1$data$Phosphate,
                   Bacteria = p1$data$Bacteria,
                   IPS1 = p1$data$Abundance,
                   ARL = p2$data$Abundance)
temp <- subset(temp, Phosphate == "minusP_30uM")

p3 <- ggplot(temp,
             aes(x = IPS1, y = ARL)) +
  #facet_wrap(~ Phosphate, ncol = 2) +
  geom_point(aes(col = Bacteria)) +
  geom_point(data = subset(temp, Bacteria %in% c("I2I3","N2N3")),
             shape = 1, col = "black", size = 2) +
  scale_y_log10() +
  scale_x_log10() +
  theme_blackbox
p3






ggplot(data.frame(Phosphate = p1$data$Phosphate,
                  Bacteria = p1$data$Bacteria,
                  PHT1.4 = p1$data$Abundance,
                  NPR1 = p2$data$Abundance),
       aes(x = PHT1.4, y = NPR1)) +
  geom_point(aes(col = Bacteria)) +
  # scale_y_log10() +
  scale_x_log10() +
  theme_blackbox



AnnotationDbi::mapIds(x = org.At.tair.db::org.At.tair.db,
                      keys = Res$Gene[1:17],column = "SYMBOL",keytype = "TAIR")


date()
