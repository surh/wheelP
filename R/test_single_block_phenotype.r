#' Test the effect of single block on phenotype
#'
#' If abundance is not used, it uses simply 1/0 for
#' presence/absence
#'
#' @export
test_single_block_phenptype <- function(Phen, abun, phenotype,
                                        variables = c("Bacteria", "Experiment", "StartP", "EndP"),
                                        confounders = c("Experiment"), use_abun = TRUE){
  singlecoms <- paste(rep(c("P","I","N"),each = 3),rep(1:3,times = 3),sep="")


  create.design <- TRUE
  if(use_abun){
    Phen <- merge_phen_and_abun(Phen = Phen, abun = abun,
                                columns = singlecoms,
                                varnames = variables)

    # Equalizing norm to compare with presence/absence model
    Phen[,singlecoms] <- Phen[,singlecoms] * 2

    create.design <- FALSE
  }

  Res.block <- NULL
  for(i in 1:2){
    for(j in 1:2){
      m1 <- block_effects(Dat = Phen, cond1 = i, cond2 = j,
                          var.name = phenotype,
                          keep.vars = confounders,
                          create.design = create.design)
      m1.sum <- summary(m1)

      res <- data.frame(SynCom = singlecoms, StartP = levels(Dat$StartP)[i],
                        EndP = levels(Dat$EndP)[j],
                        Estimate = m1.sum$coefficients[ singlecoms, 1 ],
                        SE = m1.sum$coefficients[ singlecoms, 2 ],
                        t.value = m1.sum$coefficients[ singlecoms, 3 ],
                        p.value = m1.sum$coefficients[ singlecoms, 4 ])
      Res.block <- rbind(Res.block,res)
    }
  }

  return(Res.block)
}
