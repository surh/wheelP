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

#' Test the effect of all single communities
#'
#' @export
test_single_community_phenotype <- function(Dat, dir, var.name,  bacteria.col = "Bacteria",
                                            ref.level = "none", plot = FALSE,
                                            f1.extra = "+ Experiment + Plate"){

  # bacteria.col <- "Bacteria"
  # var.name <- "Pi_content"
  # ref.level <- "none"
  # plot <- TRUE
  # f1.extra <- "+ Experiment"
  #
  if(plot){
    dir.create(dir)
  }
  communities <- levels(Dat[,bacteria.col])
  communities <- communities[ communities != ref.level ]
  Res <- NULL
  for(syncom in communities){
    #syncom <- communities[1]

    # Get data
    experiments <- unique(Dat$Experiment[ Dat[,bacteria.col] == syncom ])
    dat <- subset(Dat, Experiment %in% experiments)
    dat <- dat[ dat[,bacteria.col ]%in% c(syncom,ref.level), ]

    dat[,bacteria.col] <- relevel(droplevels(dat[,bacteria.col]), ref = ref.level)
    if("Plate" %in% colnames(dat)){
      dat$Plate <- factor(dat$Plate)
    }

    if(plot == TRUE){
      p1 <- ggplot(dat,aes_string(x = var.name,fill = bacteria.col)) +
        facet_grid(StartP ~ EndP) +
        geom_density(alpha = 0.3) +
        ggtitle(syncom) +
        theme_classic()
      #p1
      filename <- paste(dir,"/",syncom,"_density.svg",sep = "")
      ggsave(filename,p1, width = 4, height = 4)
    }

    # Test
    f1 <- paste(var.name, " ~ ", bacteria.col, f1.extra, sep = "")
    f1 <- formula(f1)

    m1 <- lm(f1,
             data = subset(dat, StartP == "-Pi,0.5%Suc" & EndP == "100 uM,0%Suc"))
    m1.sum <- summary(m1)
    res <- data.frame(SynCom = syncom,StartP = levels(dat$StartP)[1], EndP = levels(dat$EndP)[1],
                      Estimate = m1.sum$coefficients[2,1], SE = m1.sum$coefficients[2,2],
                      t.value = m1.sum$coefficients[2,3], p.value = m1.sum$coefficients[2,4])
    Res <- rbind(Res,res)

    m1 <- lm(f1,
             data = subset(dat, StartP == "+Pi,0.5%Suc" & EndP == "100 uM,0%Suc"))
    m1.sum <- summary(m1)
    res <- data.frame(SynCom = syncom,StartP = levels(dat$StartP)[2], EndP = levels(dat$EndP)[1],
                      Estimate = m1.sum$coefficients[2,1], SE = m1.sum$coefficients[2,2],
                      t.value = m1.sum$coefficients[2,3], p.value = m1.sum$coefficients[2,4])
    Res <- rbind(Res,res)
    m1 <- lm(f1,
             data = subset(dat, StartP == "-Pi,0.5%Suc" & EndP == "30 uM,0%Suc"))
    m1.sum <- summary(m1)
    res <- data.frame(SynCom = syncom,StartP = levels(dat$StartP)[1], EndP = levels(dat$EndP)[2],
                      Estimate = m1.sum$coefficients[2,1], SE = m1.sum$coefficients[2,2],
                      t.value = m1.sum$coefficients[2,3], p.value = m1.sum$coefficients[2,4])
    Res <- rbind(Res,res)
    m1 <- lm(f1,
             data = subset(dat, StartP == "+Pi,0.5%Suc" & EndP == "30 uM,0%Suc"))
    m1.sum <- summary(m1)
    res <- data.frame(SynCom = syncom,StartP = levels(dat$StartP)[2], EndP = levels(dat$EndP)[2],
                      Estimate = m1.sum$coefficients[2,1], SE = m1.sum$coefficients[2,2],
                      t.value = m1.sum$coefficients[2,3], p.value = m1.sum$coefficients[2,4])
    Res <- rbind(Res,res)
  }

  return(Res)
}
