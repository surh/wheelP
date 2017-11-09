#' Predict syncom from main effects
predict_from_main <- function(dat, Res.main){
  singlecoms <- paste(rep(c("P","I","N"),each = 3),rep(1:3,times = 3),sep="")

  #dat <- Res.sc
  dat$Block1 <- substring(dat$SynCom,1,2)
  dat$Block2 <- substring(dat$SynCom,3,4)
  dat$Block1 <- factor(dat$Block1 , levels = singlecoms)
  dat$Block2 <- factor(dat$Block2 , levels = rev(singlecoms))

  Pred <- NULL
  for(i in 1:nrow(dat)){
    index1 <- Res.main$StartP == dat$StartP[i] &
      Res.main$EndP == dat$EndP[i] &
      Res.main$SynCom == dat$Block1[i]
    index2 <- Res.main$StartP == dat$StartP[i] &
      Res.main$EndP == dat$EndP[i] &
      Res.main$SynCom == dat$Block2[i]
    additiveguess <- Res.main$Estimate[index1] + Res.main$Estimate[index2]
    seguess <- Res.main$SE[index1] + Res.main$SE[index2]

    res <- data.frame(SynCom = paste(Res.main$SynCom[index1],Res.main$SynCom[index2],sep = ""),
                      StartP = dat$StartP[i], EndP = dat$EndP[i],
                      Estimate = additiveguess, SE = seguess, t.value = NA,
                      p.value = NA, Block1 = Res.main$SynCom[index2],
                      Block2 = Res.main$SynCom[index1])
    Pred <- rbind(Pred,res)
  }

  return(Pred)
}
