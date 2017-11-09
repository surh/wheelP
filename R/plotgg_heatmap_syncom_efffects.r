#' Plot syncom effects
#'
#' DEPRECATED
plotgg_heatmap_syncom_effects <- function(Res, cond1, cond2){
  #cond1 <- 2
  #cond2 <- 2

  singlecoms <- paste(rep(c("G","N","B"),each = 3),rep(1:3,times = 3),sep="")
  dat <- subset(Res, StartP == levels(Res$StartP)[cond1] & EndP == levels(Res$EndP)[cond2])
  dat$SynCom1 <- substring(dat$SynCom,1,2)
  dat$SynCom2 <- substring(dat$SynCom,3,4)
  dat$SynCom1 <- factor(dat$SynCom1 , levels = singlecoms)
  dat$SynCom2 <- factor(dat$SynCom2 , levels = rev(singlecoms))

  p1 <- ggplot(dat, aes(x = SynCom1, y = SynCom2)) +
    geom_tile(aes(fill = Estimate, col = p.value < 0.05), size = 3, space = 2) +
    scale_fill_gradient2(low = "#d01c8b",mid = "white", high = "#4dac26",midpoint = 0,na.value = "#404040") +
    scale_color_manual(values = c("#8c510a","#01665e")) +
    ggtitle(paste(levels(dat$StartP)[cond1],"=>",levels(dat$EndP)[cond2])) +
    theme(panel.background = element_blank(),
          axis.title = element_text(face = "bold", size = 12),
          axis.text = element_text(color = "black", size = 10))

  p1
}
