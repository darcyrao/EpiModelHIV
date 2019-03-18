
# MSM -----------------------------------------------------------------

#' @title Initialization Module for the WHAMP model
#'
#' @description This function initializes the master \code{dat} object on which
#'              data are stored, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param_msm}}.
#' @param init An \code{EpiModel} object of class \code{\link{init_msm}}.
#' @param control An \code{EpiModel} object of class \code{\link{control_msm}}.
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @return
#' This function returns the updated \code{dat} object with the initialized values
#' for demographics and disease-related variables.
#'
#' @export
#' @keywords module msm
#'
initialize_msm_whamp <- function(x, param, init, control, s) {

  # Master data list
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control

  dat$attr <- list()
  dat$stats <- list()
  dat$stats$nwstats <- list()
  dat$temp <- list()
  dat$epi <- list()

  ## Network simulation ##
  nw <- list()
  for (i in 1:3) {
    nw[[i]] <- simulate(x[[i]]$fit)
    nw[[i]] <- remove_bad_roles_msm(nw[[i]])
  }

  ## Build initial edgelists
  dat$el <- list()
  dat$p <- list()
  for (i in 1:2) {
    dat$el[[i]] <- as.edgelist(nw[[i]])
    attributes(dat$el[[i]])$vnames <- NULL
    p <- tergmLite::stergm_prep(nw[[i]], x[[i]]$formation, x[[i]]$coef.diss$dissolution,
                                x[[i]]$coef.form, x[[i]]$coef.diss$coef.adj, x[[i]]$constraints)
    p$model.form$formula <- NULL
    p$model.diss$formula <- NULL
    dat$p[[i]] <- p
  }
  dat$el[[3]] <- as.edgelist(nw[[3]])
  attributes(dat$el[[3]])$vnames <- NULL
  p <- tergmLite::ergm_prep(nw[[3]], x[[3]]$formation, x[[3]]$coef.form, x[[3]]$constraints)
  p$model.form$formula <- NULL
  dat$p[[3]] <- p


  # Network parameters
  dat$nwparam <- list()
  for (i in 1:3) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }


  ## Nodal attributes ##

  # Degree terms
  dat$attr$deg.pers <- get.vertex.attribute(x[[1]]$fit$network, "deg.pers")
  dat$attr$deg.main <- get.vertex.attribute(x[[2]]$fit$network, "deg.main")

  # Race
  dat$attr$race..wa <- get.vertex.attribute(nw[[1]], "race..wa")
  num.H..wa <- dat$init$num.H..wa
  num.B..wa <- dat$init$num.B..wa
  num.O..wa <- dat$init$num.O..wa
  
  num <- num.H..wa + num.B..wa + num.O..wa
  ids.H..wa <- which(dat$attr$race..wa == "H")
  ids.B..wa <- which(dat$attr$race..wa == "B")
  ids.O..wa <- which(dat$attr$race..wa == "O")

   
  # Region
  dat$attr$region <- get.vertex.attribute(nw[[1]], "region")
  num.KC <- dat$init$num.KC
  num.OW <- dat$init$num.OW
  num.EW <- dat$init$num.EW
  ids.KC <- which(dat$attr$region == "KC")
  ids.OW <- which(dat$attr$region == "OW")
  ids.EW <- which(dat$attr$region == "EW")
    
  dat$attr$active <- rep(1, num)
  dat$attr$uid <- 1:num
  dat$temp$max.uid <- num

  # Age
  dat$attr$sqrt.age <- get.vertex.attribute(nw[[1]], "sqrt.age")
  dat$attr$age <- dat$attr$sqrt.age^2

  # Risk group
  dat$attr$riskg <- get.vertex.attribute(nw[[3]], "riskg")

  # UAI group
  p1 <- dat$param$cond.pers.always.prob * dat$param$cond.rr
  p2 <- dat$param$cond.inst.always.prob * dat$param$cond.rr
  rho <- dat$param$cond.always.prob.corr
  uai.always <- bindata::rmvbin(num, c(p1, p2), bincorr = (1 - rho) * diag(2) + rho)
  dat$attr$cond.always.pers <- uai.always[, 1]
  dat$attr$cond.always.inst <- uai.always[, 2]

  # Arrival and departure
  dat$attr$arrival.time <- rep(1, num)

  # Circumcision
  circ <- rep(NA, num)
  circ[ids.H..wa] <- sample(apportion_lr(num.H..wa, 0:1, 1 - param$circ.H.prob))
  circ[ids.B..wa] <- sample(apportion_lr(num.B..wa, 0:1, 1 - param$circ.B.prob))
  circ[ids.O..wa] <- sample(apportion_lr(num.O..wa, 0:1, 1 - param$circ.O.prob))
  dat$attr$circ <- circ

  # PrEP Attributes
  dat$attr$prepClass <- rep(NA, num)
  dat$attr$prepElig <- rep(NA, num)
  dat$attr$prepStat <- rep(0, num)
  dat$attr$prepStartTime <- rep(NA, num)
  dat$attr$prepLastRisk <- rep(NA, num)
  dat$attr$prepDisc <- rep(NA, num)
  dat$attr$spontDisc <- rep(0, num) # add attribute to track whether people discontinued PrEP spontaneously
  # dat$attr$prepLastStiScreen <- rep(NA, num)

  # Role class
  role.class <- get.vertex.attribute(nw[[1]], "role.class")
  dat$attr$role.class <- role.class

  # Ins.quot
  ins.quot <- rep(NA, num)
  ins.quot[role.class == "I"]  <- 1
  ins.quot[role.class == "R"]  <- 0
  ins.quot[role.class == "V"]  <- runif(sum(role.class == "V"))
  dat$attr$ins.quot <- ins.quot

  # HIV-related attributes
  dat <- init_status_msm_whamp(dat)

  # ## GC/CT status
  # idsUreth <- which(role.class %in% c("I", "V"))
  # idsRect <- which(role.class %in% c("R", "V"))
  # 
  # uGC <- rGC <- rep(0, num)
  # uCT <- rCT <- rep(0, num)
  # 
  # # Initialize GC infection at both sites
  # idsUGC <- sample(idsUreth, size = round(init$prev.ugc * num), FALSE)
  # uGC[idsUGC] <- 1
  # 
  # idsRGC <- sample(setdiff(idsRect, idsUGC), size = round(init$prev.rgc * num), FALSE)
  # rGC[idsRGC] <- 1
  # 
  # dat$attr$rGC <- rGC
  # dat$attr$uGC <- uGC
  # 
  # dat$attr$rGC.sympt <- dat$attr$uGC.sympt <- rep(NA, num)
  # dat$attr$rGC.sympt[rGC == 1] <- rbinom(sum(rGC == 1), 1, dat$param$rgc.sympt.prob)
  # dat$attr$uGC.sympt[uGC == 1] <- rbinom(sum(uGC == 1), 1, dat$param$ugc.sympt.prob)
  # 
  # dat$attr$rGC.infTime <- dat$attr$uGC.infTime <- rep(NA, length(dat$attr$active))
  # dat$attr$rGC.infTime[rGC == 1] <- 1
  # dat$attr$uGC.infTime[uGC == 1] <- 1
  # 
  # dat$attr$rGC.timesInf <- rep(0, num)
  # dat$attr$rGC.timesInf[rGC == 1] <- 1
  # dat$attr$uGC.timesInf <- rep(0, num)
  # dat$attr$uGC.timesInf[uGC == 1] <- 1
  # 
  # dat$attr$rGC.tx <- dat$attr$uGC.tx <- rep(NA, num)
  # dat$attr$rGC.tx.prep <- dat$attr$uGC.tx.prep <- rep(NA, num)
  # dat$attr$GC.cease <- rep(NA, num)
  # 
  # # Initialize CT infection at both sites
  # idsUCT <- sample(idsUreth, size = round(init$prev.uct * num), FALSE)
  # uCT[idsUCT] <- 1
  # 
  # idsRCT <- sample(setdiff(idsRect, idsUCT), size = round(init$prev.rct * num), FALSE)
  # rCT[idsRCT] <- 1
  # 
  # dat$attr$rCT <- rCT
  # dat$attr$uCT <- uCT
  # 
  # dat$attr$rCT.sympt <- dat$attr$uCT.sympt <- rep(NA, num)
  # dat$attr$rCT.sympt[rCT == 1] <- rbinom(sum(rCT == 1), 1, dat$param$rct.sympt.prob)
  # dat$attr$uCT.sympt[uCT == 1] <- rbinom(sum(uCT == 1), 1, dat$param$uct.sympt.prob)
  # 
  # dat$attr$rCT.infTime <- dat$attr$uCT.infTime <- rep(NA, num)
  # dat$attr$rCT.infTime[dat$attr$rCT == 1] <- 1
  # dat$attr$uCT.infTime[dat$attr$uCT == 1] <- 1
  # 
  # dat$attr$rCT.timesInf <- rep(0, num)
  # dat$attr$rCT.timesInf[rCT == 1] <- 1
  # dat$attr$uCT.timesInf <- rep(0, num)
  # dat$attr$uCT.timesInf[uCT == 1] <- 1
  # 
  # dat$attr$rCT.tx <- dat$attr$uCT.tx <- rep(NA, num)
  # dat$attr$rCT.tx.prep <- dat$attr$uCT.tx.prep <- rep(NA, num)
  # dat$attr$CT.cease <- rep(NA, num)
  # 

  # CCR5
  dat <- init_ccr5_msm_whamp(dat)


  # Network statistics
  dat$stats$nwstats <- list()


  # Prevalence Tracking
  dat$temp$deg.dists <- list()
  dat$temp$discl.list <- matrix(NA, nrow = 0, ncol = 3)
  colnames(dat$temp$discl.list) <- c("pos", "neg", "discl.time")

  dat <- prevalence_msm_whamp(dat, at = 1)

  class(dat) <- "dat"
  return(dat)
}


#' @title Removes any sexual partnerships prohibited by sexual role mixing
#'
#' @description Due to occassional issues in ERGM fitting, it is possible to
#'              have initial simulations in which there are pairings between
#'              exclusively insertive/insertive or receptive/receptive.
#'
#' @param nw An object of class \code{network}.
#'
#' @export
#' @keywords initiation utility msm
#'
remove_bad_roles_msm <- function(nw) {

  el <- as.edgelist(nw)

  rc <- get.vertex.attribute(nw, "role.class")
  rc.el <- matrix(rc[el], ncol = 2)

  rc.el.bad <- which((rc.el[, 1] == "R" & rc.el[, 2] == "R") |
                     (rc.el[, 1] == "I" & rc.el[, 2] == "I"))

  if (length(rc.el.bad) > 0) {
    el.bad <- el[rc.el.bad, , drop = FALSE]

    eid <- rep(NA, nrow(el.bad))
    for (i in 1:nrow(el.bad)) {
      eid[i] <- get.edgeIDs(nw, v = el.bad[i, 1], alter = el.bad[i, 2])
    }
    nw <- delete.edges(nw, eid)
  }

  return(nw)
}


#' @title Initialize the HIV status of persons in the network for the WHAMP model
#'
#' @description Sets the initial individual-level disease status of persons
#'              in the network, as well as disease-related attributes for
#'              infected persons. This function works by 1) assigning time since 
#'              infection as a draw from a uniform distribution with the max being
#'              the time since sexual debut, 2) assigning diagnosis based on time 
#'              since last test, which is sampled from a uniform distribution, or, 
#'              for non-screener types, time to symptomatic infection, 3) assigning treatment
#'              status, 4) assigning stage of infection based on time since infection and
#'              expected cumulative time on/off treatment, 5) assigning viral load  based on 
#'              stage of infection and treatment status
#'
#' @param dat Data object created in initialization module.
#'
#' @export
#' @keywords initiation utility msm
#'
init_status_msm_whamp <- function(dat) {
  
  age <- dat$attr$age
  race..wa <- dat$attr$race..wa
  region <- dat$attr$region
  
  ids.B..wa <- which(dat$attr$race..wa == "B")
  ids.H..wa <- which(dat$attr$race..wa == "H")
  ids.O..wa <- which(dat$attr$race..wa == "O")
  ids.KC.B <- which(dat$attr$race..wa == "B" & dat$attr$region == "KC")
  ids.KC.H <- which(dat$attr$race..wa == "H" & dat$attr$region == "KC")
  ids.KC.O <- which(dat$attr$race..wa == "O" & dat$attr$region == "KC")
  ids.OW.B <- which(dat$attr$race..wa == "B" & dat$attr$region == "OW")
  ids.OW.H <- which(dat$attr$race..wa == "H" & dat$attr$region == "OW")
  ids.OW.O <- which(dat$attr$race..wa == "O" & dat$attr$region == "OW")
  ids.EW.B <- which(dat$attr$race..wa == "B" & dat$attr$region == "EW")
  ids.EW.H <- which(dat$attr$race..wa == "H" & dat$attr$region == "EW")
  ids.EW.O <- which(dat$attr$race..wa == "O" & dat$attr$region == "EW")
  
  num.B..wa <- dat$init$num.B..wa
  num.H..wa <- dat$init$num.H..wa
  num.O..wa <- dat$init$num.O..wa
  num.KC.B <- length(ids.KC.B)
  num.KC.H <- length(ids.KC.H)
  num.KC.O <- length(ids.KC.O)
  num.OW.B <- length(ids.OW.B)
  num.OW.H <- length(ids.OW.H)
  num.OW.O <- length(ids.OW.O)
  num.EW.B <- length(ids.EW.B)
  num.EW.H <- length(ids.EW.H)
  num.EW.O <- length(ids.EW.O)
  
  num <- num.H..wa + num.B..wa + num.O..wa

  # Number initially infected by race/ethnicity
  nInfH..wa <- round(dat$init$prev.H..wa * num.H..wa)
  nInfB..wa <- round(dat$init$prev.B..wa * num.B..wa)
  nInfO..wa <- round(dat$init$prev.O..wa * num.O..wa)

  # Age-based infection probability
  probInfCrH..wa <- age[ids.H..wa] * dat$init$init.prev.age.slope.H..wa 
  probInfH..wa <- probInfCrH..wa + (nInfH..wa - sum(probInfCrH..wa)) / num.H..wa
  
  probInfCrB..wa <- age[ids.B..wa] * dat$init$init.prev.age.slope.B..wa 
  probInfB..wa <- probInfCrB..wa + (nInfB..wa - sum(probInfCrB..wa)) / num.B..wa

  probInfCrO..wa <- age[ids.O..wa] * dat$init$init.prev.age.slope.O..wa 
  probInfO..wa <- probInfCrO..wa + (nInfO..wa - sum(probInfCrO..wa)) / num.O..wa
  
  if (any(probInfH..wa <= 0) | any(probInfB..wa <= 0) | any(probInfO..wa <= 0)) { 
    stop("Slope of initial prevalence by age must be sufficiently low to ",
         "avoid non-positive probabilities.", call. = FALSE)
  }
    
    
  # Assign infection status
  status <- rep(0, num)
  while (sum(status[ids.H..wa]) != nInfH..wa) {
    status[ids.H..wa] <- rbinom(num.H..wa, 1, probInfH..wa)
  }
  while (sum(status[ids.B..wa]) != nInfB..wa) {
    status[ids.B..wa] <- rbinom(num.B..wa, 1, probInfB..wa)
  }
  while (sum(status[ids.O..wa]) != nInfO..wa) {
    status[ids.O..wa] <- rbinom(num.O..wa, 1, probInfO..wa)
  }
  dat$attr$status <- status

  # Calculate intertest interval as a function of age
  centered.age <- (age - dat$param$mean.age.iti)
  test.int <- (dat$param$iti.coefs[1] + centered.age * dat$param$iti.coefs[2] + centered.age^2 * dat$param$iti.coefs[3]) / dat$param$time.unit

  # Treatment trajectory
  tt.traj <- rep(NA, num)

  tt.traj[ids.KC.B] <- sample(apportion_lr(num.KC.B, c(1, 2, 3, 4),
                                        dat$param$tt.traj.KC.B.prob))
  tt.traj[ids.KC.H] <- sample(apportion_lr(num.KC.H, c(1, 2, 3, 4),
                                           dat$param$tt.traj.KC.H.prob))
  tt.traj[ids.KC.O] <- sample(apportion_lr(num.KC.O, c(1, 2, 3, 4),
                                           dat$param$tt.traj.KC.O.prob))
  tt.traj[ids.OW.B] <- sample(apportion_lr(num.OW.B, c(1, 2, 3, 4),
                                           dat$param$tt.traj.OW.B.prob))
  tt.traj[ids.OW.H] <- sample(apportion_lr(num.OW.H, c(1, 2, 3, 4),
                                           dat$param$tt.traj.OW.H.prob))
  tt.traj[ids.OW.O] <- sample(apportion_lr(num.OW.O, c(1, 2, 3, 4),
                                           dat$param$tt.traj.OW.O.prob))
  tt.traj[ids.EW.B] <- sample(apportion_lr(num.EW.B, c(1, 2, 3, 4),
                                           dat$param$tt.traj.EW.B.prob))
  tt.traj[ids.EW.H] <- sample(apportion_lr(num.EW.H, c(1, 2, 3, 4),
                                           dat$param$tt.traj.EW.H.prob))
  tt.traj[ids.EW.O] <- sample(apportion_lr(num.EW.O, c(1, 2, 3, 4),
                                           dat$param$tt.traj.EW.O.prob))
  dat$attr$tt.traj <- tt.traj


  ## Infection-related attributes
  stage <- rep(NA, num)
  stage.time <- rep(NA, num)
  inf.time <- rep(NA, num)
  vl <- rep(NA, num)
  diag.status <- rep(NA, num)
  diag.time <- rep(NA, num)
  last.neg.test <- rep(NA, num)
  tx.status <- rep(NA, num)
  tx.init.time <- rep(NA, num)
  cum.time.on.tx <- rep(NA, num)
  cum.time.off.tx <- rep(NA, num)
  infector <- rep(NA, num)
  inf.role <- rep(NA, num)
  inf.type <- rep(NA, num)
  inf.diag <- rep(NA, num)
  inf.tx <- rep(NA, num)
  inf.stage <- rep(NA, num)
  
  ## Define vectors and values used to calculate and assign the above attributes
  time.since.inf <- rep(NA, num)
  time.to.dx <- rep(NA, num)
  time.to.tx <- rep(NA, num)
  prop.time.on.tx <- rep(NA, num)
 
  vlar.int <- dat$param$vl.acute.rise.int
  vlap <- dat$param$vl.acute.peak
  vlaf.int <- dat$param$vl.acute.fall.int
  vlsp <- dat$param$vl.set.point
  vldo.int <- dat$param$vl.aids.onset.int
  vl.aids.int <- dat$param$vl.aids.int
  vlaidsp  <- dat$param$vl.aids.peak
  sympt.int <- dat$param$sympt.onset.int
  vlds <- (vlaidsp - vlsp) / vl.aids.int
  vl.acute.int <- vlar.int + vlaf.int
  tx.halt.full <- dat$param$tx.halt.full
  tx.halt.part.rr <- dat$param$tx.halt.part.rr
  tx.reinit.part.rr <- dat$param$tx.reinit.part.rr
  
  tx.init.int <- rep(NA, num)
  tx.init.int[race..wa == "B" & region == "KC"] <- dat$param$tx.init.int.KC.B
  tx.init.int[race..wa == "H" & region == "KC"] <- dat$param$tx.init.int.KC.H
  tx.init.int[race..wa == "O" & region == "KC"] <- dat$param$tx.init.int.KC.O
  tx.init.int[race..wa == "B" & region == "OW"] <- dat$param$tx.init.int.OW.B
  tx.init.int[race..wa == "H" & region == "OW"] <- dat$param$tx.init.int.OW.H
  tx.init.int[race..wa == "O" & region == "OW"] <- dat$param$tx.init.int.OW.O
  tx.init.int[race..wa == "B" & region == "EW"] <- dat$param$tx.init.int.EW.B
  tx.init.int[race..wa == "H" & region == "EW"] <- dat$param$tx.init.int.EW.H
  tx.init.int[race..wa == "O" & region == "EW"] <- dat$param$tx.init.int.EW.O
  
  tx.reinit.full <- rep(NA, num)
  tx.reinit.full[race..wa == "B" & region == "KC"] <- dat$param$tx.reinit.full.KC.B
  tx.reinit.full[race..wa == "H" & region == "KC"] <- dat$param$tx.reinit.full.KC.H
  tx.reinit.full[race..wa == "O" & region == "KC"] <- dat$param$tx.reinit.full.KC.O
  tx.reinit.full[race..wa == "B" & region == "OW"] <- dat$param$tx.reinit.full.OW.B
  tx.reinit.full[race..wa == "H" & region == "OW"] <- dat$param$tx.reinit.full.OW.H
  tx.reinit.full[race..wa == "O" & region == "OW"] <- dat$param$tx.reinit.full.OW.O
  tx.reinit.full[race..wa == "B" & region == "EW"] <- dat$param$tx.reinit.full.EW.B
  tx.reinit.full[race..wa == "H" & region == "EW"] <- dat$param$tx.reinit.full.EW.H
  tx.reinit.full[race..wa == "O" & region == "EW"] <- dat$param$tx.reinit.full.EW.O

  time.sex.active <- pmax(1,
                          round((365 / dat$param$time.unit) * age - (365 / dat$param$time.unit) *
                                    min(dat$init$ages), 0))
  
  ### 1. Non-screener type - partial suppressors
  selected <- which(status == 1 & tt.traj == 1)
 
  ## Assign infection time
  time.since.inf[selected] <- ceiling(runif(length(selected), max = time.sex.active[selected])) # set max time since infection to time since sexual debut (assumed to be age 18)
  inf.time[selected] <- 1 - time.since.inf[selected]
  
  ## Assign diagnosis status
  diag.status[selected] <- 0
  diag.status[selected][time.since.inf[selected] >= sympt.int] <- 1
  time.to.dx[selected][diag.status[selected] == 1] <- sympt.int
  diag.time[selected][diag.status[selected] == 1] <- 1 - time.since.inf[selected][diag.status[selected] == 1] + sympt.int
  
  ## Assign treatment status
  time.to.tx[selected] <- sympt.int + tx.init.int[selected]
  prop.time.on.tx[selected] <- (tx.reinit.full[selected] * tx.reinit.part.rr) /
    ((tx.halt.full * tx.halt.part.rr) + (tx.reinit.full[selected] * tx.reinit.part.rr))
  tx.status[selected] <- 0
  tx.status[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])] <-
    rbinom(sum(time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])),
           1, prop.time.on.tx[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])])
  tx.init.time[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])] <- 
    1 - time.since.inf[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])] + 
    time.to.tx[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])]

  ### 2. Non-screener type - full suppressors
  selected <- which(status == 1 & tt.traj == 2)
  
  ## Assign infection time
  time.since.inf[selected] <- ceiling(runif(length(selected), max = time.sex.active[selected])) # set max time since infection to time since sexual debut (assumed to be age 18)
  inf.time[selected] <- 1 - time.since.inf[selected]
  
  ## Assign diagnosis status
  diag.status[selected] <- 0
  diag.status[selected][time.since.inf[selected] >= sympt.int] <- 1
  time.to.dx[selected][diag.status[selected] == 1] <- sympt.int
  diag.time[selected][diag.status[selected] == 1] <- 1 - time.since.inf[selected][diag.status[selected] == 1] + sympt.int
  
  ## Assign treatment status
  time.to.tx[selected] <- sympt.int + tx.init.int[selected]
  prop.time.on.tx[selected] <- (tx.reinit.full[selected]) /
    ((tx.halt.full) + (tx.reinit.full[selected]))
  tx.status[selected] <- 0
  tx.status[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])] <-
    rbinom(sum(time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])),
           1, prop.time.on.tx[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])])
  tx.init.time[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])] <- 
    1 - time.since.inf[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])] + 
    time.to.tx[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])]
  
  ### 3. Regular screeners - partial suppressors
  selected <- which(status == 1 & tt.traj == 3)
  
  ## Assign infection time
  time.since.inf[selected] <- ceiling(runif(length(selected), max = time.sex.active[selected])) # set max time since infection to time since sexual debut (assumed to be age 18)
  inf.time[selected] <- 1 - time.since.inf[selected]
  
  ## Assign diagnosis status
  if (dat$param$testing.pattern == "interval") {
      tslt <- ceiling(runif(length(selected),
                            min = 0,
                            max = test.int[selected]))
  }
  if (dat$param$testing.pattern == "memoryless") {
      stop("Intertest interval parameter calculated assuming interval method. Revise parameter estimation procedure for memoryless process.")
  }
  
  diag.status[selected][tslt > time.since.inf[selected]] <- 0 # Last test was before infection
  diag.status[selected][time.since.inf[selected] < dat$param$test.window.int] <- 0 # Infection < test window period
  diag.status[selected][tslt <= time.since.inf[selected] & tslt > (time.since.inf[selected] - dat$param$test.window.int)] <- 0 # Last test < test window period
  diag.status[selected][tslt <= (time.since.inf[selected] - dat$param$test.window.int)] <- 1 # Last test occurred since infection but after the window period

  ##' Time to dx - calculate time to first test based on test interval (assigned last test used for determining if diagnosis would have occurred may not have been first pos test 
  ##' if infection occurred long ago. So calculate time to first post test as: time since infection - (time since last test + floor((time since infection - time since last test - test.window) / test.int) * test.int)
  time.to.dx[selected][diag.status[selected] == 1] <- time.since.inf[selected][diag.status[selected] == 1] - 
    (tslt[diag.status[selected] == 1] + floor((time.since.inf[selected][diag.status[selected] == 1] - tslt[diag.status[selected] == 1] - dat$param$test.window.int)/test.int[selected][diag.status[selected] == 1]) * test.int[selected][diag.status[selected] == 1])  
  
  diag.time[selected][diag.status[selected] == 1] <- 1 - time.since.inf[selected][diag.status[selected] == 1] + time.to.dx[selected][diag.status[selected] == 1]
  
  last.neg.test[selected][diag.status[selected] == 0] <- -tslt[diag.status[selected] == 0]
  last.neg.test[selected][diag.status[selected] == 1] <- NA
  
  # Assign treatment status
  time.to.tx[selected][diag.status[selected] == 1] <- time.to.dx[selected][diag.status[selected] == 1] + tx.init.int[selected][diag.status[selected] == 1]
  prop.time.on.tx[selected] <- (tx.reinit.full[selected] * tx.reinit.part.rr) /
    ((tx.halt.full * tx.halt.part.rr) + (tx.reinit.full[selected] * tx.reinit.part.rr))
  tx.status[selected] <- 0
  tx.status[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])] <-
    rbinom(sum(time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])),
           1, prop.time.on.tx[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])])
  tx.init.time[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])] <- 
    1 - time.since.inf[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])] + 
    time.to.tx[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])]
  
  ### 4. Regular screeners - full suppressors
  selected <- which(status == 1 & tt.traj == 4)
  
  ## Assign infection time
  time.since.inf[selected] <- ceiling(runif(length(selected), max = time.sex.active[selected])) # set max time since infection to time since sexual debut (assumed to be age 18)
  inf.time[selected] <- 1 - time.since.inf[selected]
  
  ## Assign diagnosis status
  if (dat$param$testing.pattern == "interval") {
      tslt <- ceiling(runif(length(selected),
                            min = 0,
                            max = test.int[selected]))
  }
  if (dat$param$testing.pattern == "memoryless") {
      stop("Intertest interval parameter calculated assuming interval method. Revise parameter estimation procedure for memoryless process.")
  }
  
  diag.status[selected][tslt > time.since.inf[selected]] <- 0 # Last test before infection
  diag.status[selected][time.since.inf[selected] < dat$param$test.window.int] <- 0 # Infection < test window period
  diag.status[selected][tslt <= time.since.inf[selected] & tslt > (time.since.inf[selected] - dat$param$test.window.int)] <- 0 # Last test < test window period
  diag.status[selected][tslt <= (time.since.inf[selected] - dat$param$test.window.int)] <- 1 # Last test occurred since infection but after the window period
  
  ##' Time to dx - calculate time to first test based on test interval (assigned last test used for determining if diagnosis would have occurred may not have been first pos test 
  ##' if infection occurred long ago. So calculate time to first post test as: time since infection - (time since last test + floor((time since infection - time since last test - test.window) / test.int) * test.int)
  time.to.dx[selected][diag.status[selected] == 1] <- time.since.inf[selected][diag.status[selected] == 1] - 
    (tslt[diag.status[selected] == 1] + floor((time.since.inf[selected][diag.status[selected] == 1] - tslt[diag.status[selected] == 1] - dat$param$test.window.int)/test.int[selected][diag.status[selected] == 1]) * test.int[selected][diag.status[selected] == 1])  
  
  diag.time[selected][diag.status[selected] == 1] <- 1 - time.since.inf[selected][diag.status[selected] == 1] + time.to.dx[selected][diag.status[selected] == 1]
  
  last.neg.test[selected][diag.status[selected] == 0] <- -tslt[diag.status[selected] == 0]
  last.neg.test[selected][diag.status[selected] == 1] <- NA
  
  # Assign treatment status
  time.to.tx[selected][diag.status[selected] == 1] <- time.to.dx[selected][diag.status[selected] == 1] + tx.init.int[selected][diag.status[selected] == 1]
  prop.time.on.tx[selected] <- (tx.reinit.full[selected]) /
      ((tx.halt.full) + (tx.reinit.full[selected]))
  tx.status[selected] <- 0
  tx.status[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])] <-
    rbinom(sum(time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])),
           1, prop.time.on.tx[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])])
  tx.init.time[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])] <- 
    1 - time.since.inf[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])] + 
    time.to.tx[selected][time.since.inf[selected] >= time.to.tx[selected] & !is.na(time.to.tx[selected])]
  
  ### All tt.traj groups: Calculate cumulative time on and off tx to determine progression to AIDS
  selected <- which(status == 1)
  
  # If not yet diagnosed or initiated treatment, cum.time.off = time since infection and cum.time.on = 0
  cum.time.off.tx[selected][diag.status[selected] == 0] <- time.since.inf[selected][diag.status[selected] == 0]
  cum.time.off.tx[selected][time.since.inf[selected] < time.to.tx[selected] & !is.na(time.to.tx[selected])] <- 
    time.since.inf[selected][time.since.inf[selected] < time.to.tx[selected] & !is.na(time.to.tx[selected])]
  cum.time.on.tx[selected][diag.status[selected] == 0] <- 0
  cum.time.on.tx[selected][time.since.inf[selected] < time.to.tx[selected] & !is.na(time.to.tx[selected])] <- 0
  
  # If initiate treatment in the current time step, cum.time.off = 0 time since infection and cum.time.on = 0
  cum.time.off.tx[selected][time.since.inf[selected] == time.to.tx[selected] & !is.na(time.to.tx[selected])] <- 
    time.since.inf[selected][time.since.inf[selected] == time.to.tx[selected] & !is.na(time.to.tx[selected])]
  cum.time.on.tx[selected][time.since.inf[selected] == time.to.tx[selected] & !is.na(time.to.tx[selected])] <- 0
  
  # If intitated tx, cum.time.off = time to tx initiation + expected cuml time off since starting tx, and cum.time.on = expected cuml time on since starting tx
  #' Expected cumulative time off since starting is estimated with offon: column 1 indicates the expected cumulative time off treatment and column 2 indicates the
  #' expected cumulative time on treatment at each time point since treatment was initiated. We calculate each of these vectors assuming a Markov process based on 
  #' the calculated proportion of time on treatment, and we define the "offon" matrix separately for each group defined by race, region, and full vs. partial 
  #' suppressoin status, since the proportion of time on treatment varies for each of these groups.
  
  # First define "offon" matrices
  max.steps <- max((dat$param$max.time.off.tx.int - tx.init.int[selected]) / (1 - max(prop.time.on.tx[selected]))) # determine the max number of steps to extend out the vector to make sure it covers everyone's possible time from tx to AIDS
  offon.KC.B.part <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "B" & region[selected] == "KC" & tt.traj[selected] %in% c(1,3)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "B" & region[selected] == "KC" & tt.traj[selected] %in% c(1,3)][1] * 1:max.steps)
  offon.KC.B.full <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "B" & region[selected] == "KC" & tt.traj[selected] %in% c(2,4)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "B" & region[selected] == "KC" & tt.traj[selected] %in% c(2,4)][1] * 1:max.steps)
  offon.KC.H.part <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "H" & region[selected] == "KC" & tt.traj[selected] %in% c(1,3)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "H" & region[selected] == "KC" & tt.traj[selected] %in% c(1,3)][1] * 1:max.steps)
  offon.KC.H.full <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "H" & region[selected] == "KC" & tt.traj[selected] %in% c(2,4)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "H" & region[selected] == "KC" & tt.traj[selected] %in% c(2,4)][1] * 1:max.steps)
  offon.KC.O.part <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "O" & region[selected] == "KC" & tt.traj[selected] %in% c(1,3)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "O" & region[selected] == "KC" & tt.traj[selected] %in% c(1,3)][1] * 1:max.steps)
  offon.KC.O.full <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "O" & region[selected] == "KC" & tt.traj[selected] %in% c(2,4)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "O" & region[selected] == "KC" & tt.traj[selected] %in% c(2,4)][1] * 1:max.steps)
  offon.OW.B.part <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "B" & region[selected] == "OW" & tt.traj[selected] %in% c(1,3)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "B" & region[selected] == "OW" & tt.traj[selected] %in% c(1,3)][1] * 1:max.steps)
  offon.OW.B.full <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "B" & region[selected] == "OW" & tt.traj[selected] %in% c(2,4)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "B" & region[selected] == "OW" & tt.traj[selected] %in% c(2,4)][1] * 1:max.steps)
  offon.OW.H.part <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "H" & region[selected] == "OW" & tt.traj[selected] %in% c(1,3)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "H" & region[selected] == "OW" & tt.traj[selected] %in% c(1,3)][1] * 1:max.steps)
  offon.OW.H.full <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "H" & region[selected] == "OW" & tt.traj[selected] %in% c(2,4)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "H" & region[selected] == "OW" & tt.traj[selected] %in% c(2,4)][1] * 1:max.steps)
  offon.OW.O.part <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "O" & region[selected] == "OW" & tt.traj[selected] %in% c(1,3)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "O" & region[selected] == "OW" & tt.traj[selected] %in% c(1,3)][1] * 1:max.steps)
  offon.OW.O.full <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "O" & region[selected] == "OW" & tt.traj[selected] %in% c(2,4)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "O" & region[selected] == "OW" & tt.traj[selected] %in% c(2,4)][1] * 1:max.steps)
  offon.EW.B.part <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "B" & region[selected] == "EW" & tt.traj[selected] %in% c(1,3)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "B" & region[selected] == "EW" & tt.traj[selected] %in% c(1,3)][1] * 1:max.steps)
  offon.EW.B.full <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "B" & region[selected] == "EW" & tt.traj[selected] %in% c(2,4)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "B" & region[selected] == "EW" & tt.traj[selected] %in% c(2,4)][1] * 1:max.steps)
  offon.EW.H.part <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "H" & region[selected] == "EW" & tt.traj[selected] %in% c(1,3)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "H" & region[selected] == "EW" & tt.traj[selected] %in% c(1,3)][1] * 1:max.steps)
  offon.EW.H.full <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "H" & region[selected] == "EW" & tt.traj[selected] %in% c(2,4)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "H" & region[selected] == "EW" & tt.traj[selected] %in% c(2,4)][1] * 1:max.steps)
  offon.EW.O.part <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "O" & region[selected] == "EW" & tt.traj[selected] %in% c(1,3)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "O" & region[selected] == "EW" & tt.traj[selected] %in% c(1,3)][1] * 1:max.steps)
  offon.EW.O.full <- cbind((1 - prop.time.on.tx[selected][race..wa[selected] == "O" & region[selected] == "EW" & tt.traj[selected] %in% c(2,4)][1]) * 1:max.steps, 
                           prop.time.on.tx[selected][race..wa[selected] == "O" & region[selected] == "EW" & tt.traj[selected] %in% c(2,4)][1] * 1:max.steps)
  
  # Define new indices to subset based on having started treatment 1+ time steps ago, being in the specified group defined by race/ethnicity, region, and partial vs. full adherence
  selected.on.tx.KC.B.part <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "B" & 
                                         region[selected] == "KC" & tt.traj[selected] %in% c(1,3)]
  selected.on.tx.KC.B.full <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "B" & 
                                         region[selected] == "KC" & tt.traj[selected] %in% c(2,4)]
  selected.on.tx.KC.H.part <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "H" & 
                                         region[selected] == "KC" & tt.traj[selected] %in% c(1,3)]
  selected.on.tx.KC.H.full <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "H" & 
                                         region[selected] == "KC" & tt.traj[selected] %in% c(2,4)]
  selected.on.tx.KC.O.part <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "O" & 
                                         region[selected] == "KC" & tt.traj[selected] %in% c(1,3)]
  selected.on.tx.KC.O.full <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "O" & 
                                         region[selected] == "KC" & tt.traj[selected] %in% c(2,4)]
  selected.on.tx.OW.B.part <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "B" & 
                                         region[selected] == "OW" & tt.traj[selected] %in% c(1,3)]
  selected.on.tx.OW.B.full <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "B" & 
                                         region[selected] == "OW" & tt.traj[selected] %in% c(2,4)]
  selected.on.tx.OW.H.part <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "H" & 
                                         region[selected] == "OW" & tt.traj[selected] %in% c(1,3)]
  selected.on.tx.OW.H.full <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "H" & 
                                         region[selected] == "OW" & tt.traj[selected] %in% c(2,4)]
  selected.on.tx.OW.O.part <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "O" & 
                                         region[selected] == "OW" & tt.traj[selected] %in% c(1,3)]
  selected.on.tx.OW.O.full <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "O" & 
                                         region[selected] == "OW" & tt.traj[selected] %in% c(2,4)]
  selected.on.tx.EW.B.part <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "B" & 
                                         region[selected] == "EW" & tt.traj[selected] %in% c(1,3)]
  selected.on.tx.EW.B.full <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "B" & 
                                         region[selected] == "EW" & tt.traj[selected] %in% c(2,4)]
  selected.on.tx.EW.H.part <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "H" & 
                                         region[selected] == "EW" & tt.traj[selected] %in% c(1,3)]
  selected.on.tx.EW.H.full <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "H" & 
                                         region[selected] == "EW" & tt.traj[selected] %in% c(2,4)]
  selected.on.tx.EW.O.part <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "O" & 
                                         region[selected] == "EW" & tt.traj[selected] %in% c(1,3)]
  selected.on.tx.EW.O.full <- selected[time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]) & race..wa[selected] == "O" & 
                                         region[selected] == "EW" & tt.traj[selected] %in% c(2,4)]
  
  #' For each group defined by race/ethnicity, region, and partial vs. full adherence, calculate the expected cumulative time off/on treatment using the
  #' offon matrices and time since treatment initiation + time from infection to treatment
  cum.time.off.tx[selected.on.tx.KC.B.part] <- round(time.to.tx[selected.on.tx.KC.B.part] + offon.KC.B.part[(-tx.init.time[selected.on.tx.KC.B.part] + 1), 1])
  cum.time.on.tx[selected.on.tx.KC.B.part] <- round(offon.KC.B.part[(-tx.init.time[selected.on.tx.KC.B.part] + 1), 2])
  cum.time.off.tx[selected.on.tx.KC.B.full] <- round(time.to.tx[selected.on.tx.KC.B.full] + offon.KC.B.full[(-tx.init.time[selected.on.tx.KC.B.full] + 1), 1])
  cum.time.on.tx[selected.on.tx.KC.B.full] <- round(offon.KC.B.full[(-tx.init.time[selected.on.tx.KC.B.full] + 1), 2])
  cum.time.off.tx[selected.on.tx.KC.H.part] <- round(time.to.tx[selected.on.tx.KC.H.part] + offon.KC.H.part[(-tx.init.time[selected.on.tx.KC.H.part] + 1), 1])
  cum.time.on.tx[selected.on.tx.KC.H.part] <- round(offon.KC.H.part[(-tx.init.time[selected.on.tx.KC.H.part] + 1), 2])
  cum.time.off.tx[selected.on.tx.KC.H.full] <- round(time.to.tx[selected.on.tx.KC.H.full] + offon.KC.H.full[(-tx.init.time[selected.on.tx.KC.H.full] + 1), 1])
  cum.time.on.tx[selected.on.tx.KC.H.full] <- round(offon.KC.H.full[(-tx.init.time[selected.on.tx.KC.H.full] + 1), 2])
  cum.time.off.tx[selected.on.tx.KC.O.part] <- round(time.to.tx[selected.on.tx.KC.O.part] + offon.KC.O.part[(-tx.init.time[selected.on.tx.KC.O.part] + 1), 1])
  cum.time.on.tx[selected.on.tx.KC.O.part] <- round(offon.KC.O.part[(-tx.init.time[selected.on.tx.KC.O.part] + 1), 2])
  cum.time.off.tx[selected.on.tx.KC.O.full] <- round(time.to.tx[selected.on.tx.KC.O.full] + offon.KC.O.full[(-tx.init.time[selected.on.tx.KC.O.full] + 1), 1])
  cum.time.on.tx[selected.on.tx.KC.O.full] <- round(offon.KC.O.full[(-tx.init.time[selected.on.tx.KC.O.full] + 1), 2])
  
  cum.time.off.tx[selected.on.tx.OW.B.part] <- round(time.to.tx[selected.on.tx.OW.B.part] + offon.OW.B.part[(-tx.init.time[selected.on.tx.OW.B.part] + 1), 1])
  cum.time.on.tx[selected.on.tx.OW.B.part] <- round(offon.OW.B.part[(-tx.init.time[selected.on.tx.OW.B.part] + 1), 2])
  cum.time.off.tx[selected.on.tx.OW.B.full] <- round(time.to.tx[selected.on.tx.OW.B.full] + offon.OW.B.full[(-tx.init.time[selected.on.tx.OW.B.full] + 1), 1])
  cum.time.on.tx[selected.on.tx.OW.B.full] <- round(offon.OW.B.full[(-tx.init.time[selected.on.tx.OW.B.full] + 1), 2])
  cum.time.off.tx[selected.on.tx.OW.H.part] <- round(time.to.tx[selected.on.tx.OW.H.part] + offon.OW.H.part[(-tx.init.time[selected.on.tx.OW.H.part] + 1), 1])
  cum.time.on.tx[selected.on.tx.OW.H.part] <- round(offon.OW.H.part[(-tx.init.time[selected.on.tx.OW.H.part] + 1), 2])
  cum.time.off.tx[selected.on.tx.OW.H.full] <- round(time.to.tx[selected.on.tx.OW.H.full] + offon.OW.H.full[(-tx.init.time[selected.on.tx.OW.H.full] + 1), 1])
  cum.time.on.tx[selected.on.tx.OW.H.full] <- round(offon.OW.H.full[(-tx.init.time[selected.on.tx.OW.H.full] + 1), 2])
  cum.time.off.tx[selected.on.tx.OW.O.part] <- round(time.to.tx[selected.on.tx.OW.O.part] + offon.OW.O.part[(-tx.init.time[selected.on.tx.OW.O.part] + 1), 1])
  cum.time.on.tx[selected.on.tx.OW.O.part] <- round(offon.OW.O.part[(-tx.init.time[selected.on.tx.OW.O.part] + 1), 2])
  cum.time.off.tx[selected.on.tx.OW.O.full] <- round(time.to.tx[selected.on.tx.OW.O.full] + offon.OW.O.full[(-tx.init.time[selected.on.tx.OW.O.full] + 1), 1])
  cum.time.on.tx[selected.on.tx.OW.O.full] <- round(offon.OW.O.full[(-tx.init.time[selected.on.tx.OW.O.full] + 1), 2])
  
  cum.time.off.tx[selected.on.tx.EW.B.part] <- round(time.to.tx[selected.on.tx.EW.B.part] + offon.EW.B.part[(-tx.init.time[selected.on.tx.EW.B.part] + 1), 1])
  cum.time.on.tx[selected.on.tx.EW.B.part] <- round(offon.EW.B.part[(-tx.init.time[selected.on.tx.EW.B.part] + 1), 2])
  cum.time.off.tx[selected.on.tx.EW.B.full] <- round(time.to.tx[selected.on.tx.EW.B.full] + offon.EW.B.full[(-tx.init.time[selected.on.tx.EW.B.full] + 1), 1])
  cum.time.on.tx[selected.on.tx.EW.B.full] <- round(offon.EW.B.full[(-tx.init.time[selected.on.tx.EW.B.full] + 1), 2])
  cum.time.off.tx[selected.on.tx.EW.H.part] <- round(time.to.tx[selected.on.tx.EW.H.part] + offon.EW.H.part[(-tx.init.time[selected.on.tx.EW.H.part] + 1), 1])
  cum.time.on.tx[selected.on.tx.EW.H.part] <- round(offon.EW.H.part[(-tx.init.time[selected.on.tx.EW.H.part] + 1), 2])
  cum.time.off.tx[selected.on.tx.EW.H.full] <- round(time.to.tx[selected.on.tx.EW.H.full] + offon.EW.H.full[(-tx.init.time[selected.on.tx.EW.H.full] + 1), 1])
  cum.time.on.tx[selected.on.tx.EW.H.full] <- round(offon.EW.H.full[(-tx.init.time[selected.on.tx.EW.H.full] + 1), 2])
  cum.time.off.tx[selected.on.tx.EW.O.part] <- round(time.to.tx[selected.on.tx.EW.O.part] + offon.EW.O.part[(-tx.init.time[selected.on.tx.EW.O.part] + 1), 1])
  cum.time.on.tx[selected.on.tx.EW.O.part] <- round(offon.EW.O.part[(-tx.init.time[selected.on.tx.EW.O.part] + 1), 2])
  cum.time.off.tx[selected.on.tx.EW.O.full] <- round(time.to.tx[selected.on.tx.EW.O.full] + offon.EW.O.full[(-tx.init.time[selected.on.tx.EW.O.full] + 1), 1])
  cum.time.on.tx[selected.on.tx.EW.O.full] <- round(offon.EW.O.full[(-tx.init.time[selected.on.tx.EW.O.full] + 1), 2])
  
  #' For each group, calculate the expected time step from infection at which AIDS will occur based on the offon matrix and time from infection to tx initiation
  exp.aids.time <- rep(NA, length(selected))
  
  calc_aids_timestep <- function(group) {
    selected <- get(paste0("selected.on.tx.", group))
    offon <- get(paste0("offon.", group))
    if (length(selected) > 0) {
      for (i in 1:length(selected)) {
        exp.aids.time[selected][i] <- time.to.tx[selected][i] + 
          min(which(offon[ ,1] > (dat$param$max.time.off.tx.int - time.to.tx[selected][i])))
      }
    }
    return(exp.aids.time)
  }
  
  exp.aids.time <- calc_aids_timestep("KC.B.part")
  exp.aids.time <- calc_aids_timestep("KC.B.full")
  exp.aids.time <- calc_aids_timestep("KC.H.part")
  exp.aids.time <- calc_aids_timestep("KC.H.full")
  exp.aids.time <- calc_aids_timestep("KC.O.part")
  exp.aids.time <- calc_aids_timestep("KC.O.full")
  
  exp.aids.time <- calc_aids_timestep("OW.B.part")
  exp.aids.time <- calc_aids_timestep("OW.B.full")
  exp.aids.time <- calc_aids_timestep("OW.H.part")
  exp.aids.time <- calc_aids_timestep("OW.H.full")
  exp.aids.time <- calc_aids_timestep("OW.O.part")
  exp.aids.time <- calc_aids_timestep("OW.O.full")
  
  exp.aids.time <- calc_aids_timestep("EW.B.part")
  exp.aids.time <- calc_aids_timestep("EW.B.full")
  exp.aids.time <- calc_aids_timestep("EW.H.part")
  exp.aids.time <- calc_aids_timestep("EW.H.full")
  exp.aids.time <- calc_aids_timestep("EW.O.part")
  exp.aids.time <- calc_aids_timestep("EW.O.full")
  
  ### All tt.traj groups: Assign stage of infection 
  stage[selected] <- (time.since.inf[selected] <= vlar.int) * 1 +
    (time.since.inf[selected] > vlar.int) * (time.since.inf[selected] <= vl.acute.int) * 2 +
    (time.since.inf[selected] > vl.acute.int) * (cum.time.off.tx[selected] < dat$param$max.time.off.tx.int) * 3 +
    (time.since.inf[selected] > vl.acute.int) * (cum.time.off.tx[selected] >= dat$param$max.time.off.tx.int) * 4
  stage.time[selected][stage[selected] == 1] <- time.since.inf[selected][stage[selected] == 1]
  stage.time[selected][stage[selected] == 2] <- time.since.inf[selected][stage[selected] == 2] - vlar.int
  stage.time[selected][stage[selected] == 3] <- time.since.inf[selected][stage[selected] == 3] - vl.acute.int
  stage.time[selected][stage[selected] == 4 & (time.since.inf[selected] <= time.to.tx[selected] | is.na(time.to.tx[selected]))] <- 
    time.since.inf[selected][stage[selected] == 4 & (time.since.inf[selected] <= time.to.tx[selected] | is.na(time.to.tx[selected]))] - dat$param$max.time.off.tx.int # If not yet initiated tx (or just initiated) and in stage 4, stage.time = time since infection - number of days off treatment for a before onset of AIDS
  stage.time[selected][stage[selected] == 4 & (time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]))] <- 
    time.since.inf[selected][stage[selected] == 4 & (time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]))] - 
    exp.aids.time[selected][stage[selected] == 4 & (time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected]))]
  
  
  ### All tt.traj groups: Assign VL (assuming a linear rate of change in VL up to peak viremia in acute phase and from peak down to set point)
  
  ##' To set VL in the AIDS phase, define a variable average consecutive time off tx
  ##' - For those who have not yet inititated tx, set cons.time.off.tx to time since infected.
  ##' - For those who disontinued, set to avg time off before reinitiating
  cons.time.off.tx <- rep(NA, length(selected))
  cons.time.off.tx[tx.status[selected] == 1] <- 0
  cons.time.off.tx[diag.status[selected] == 0] <- time.since.inf[selected][diag.status[selected] == 0]
  cons.time.off.tx[time.since.inf[selected] <= time.to.tx[selected] & !is.na(time.to.tx[selected])] <- 
    time.since.inf[selected][time.since.inf[selected] <= time.to.tx[selected] & !is.na(time.to.tx[selected])]
  cons.time.off.tx[tx.status[selected] == 0 & (time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected])) & tt.traj[selected] %in% c(1,3)] <- 
    (1 / (tx.reinit.full[selected][tx.status[selected] == 0 & (time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected])) & tt.traj[selected] %in% c(1,3)] * tx.reinit.part.rr))
  cons.time.off.tx[tx.status[selected] == 0 & (time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected])) & tt.traj[selected] %in% c(2,4)] <- 
    (1 / (tx.reinit.full[selected][tx.status[selected] == 0 & (time.since.inf[selected] > time.to.tx[selected] & !is.na(time.to.tx[selected])) & tt.traj[selected] %in% c(2,4)]))
  
  ##' Assume that, upon stopping ART in the AIDS phase, VL rises at the same rate as it would in the chronic phase up to set point and then
  ##' rises at the slope that VL rises in untreated AIDS to go from set point to the max level in 2 years. To implement this, first define 
  ##' the time it takes for VL to go back to set point. Because HIV-related mortality is modeled as a constant increased hazard of death rather
  ##' than upon progression to the end of AIDS, VL is constrained not to exceed 7 log10. 
  vl.supp <- rep(NA, length(selected))
  vl.up.slope <- rep(NA, length(selected))
  vl.supp[tt.traj[selected] %in% c(1,3)] <- dat$param$vl.part.supp
  vl.supp[tt.traj[selected] %in% c(2,4)] <- dat$param$vl.full.supp
  vl.up.slope[tt.traj[selected] %in% c(1,3)] <- dat$param$part.supp.up.slope
  vl.up.slope[tt.traj[selected] %in% c(2,4)] <- dat$param$full.supp.up.slope
  
  vl.rise.int <- rep(NA, length(selected))
  vl.rise.int <- (vlsp - vl.supp)/vl.up.slope
  
  ## if stage 4, set vl according to slope of change in VL and time since tx. If in stage 3, assume all off tx are at set point
  vl[selected] <- (time.since.inf[selected] <= vlar.int) * (vlap * time.since.inf[selected] / vlar.int) +
    (time.since.inf[selected] > vlar.int) * (time.since.inf[selected] <= vl.acute.int) *
    ((vlsp - vlap) * (time.since.inf[selected] - vlar.int) / vlaf.int + vlap) +
    (time.since.inf[selected] > vl.acute.int) * (cum.time.off.tx[selected] < dat$param$max.time.off.tx.int) * vlsp +
    (cum.time.off.tx[selected] >= dat$param$max.time.off.tx.int) * (cons.time.off.tx < vl.rise.int) * 
    (vl.supp + (cons.time.off.tx * vl.up.slope)) + 
    (cum.time.off.tx[selected] >= dat$param$max.time.off.tx.int) * (cons.time.off.tx >= vl.rise.int) * 
    pmin((vlsp + (cons.time.off.tx - vl.rise.int) * vlds), vlaidsp)
  ##' We assume that all who are on tx are at their suppressed levels. This doesn't account for the fact that some people may have 
  ##' recently reinitiated and not yet achieved on-treatment levels. But we assume all men in the chronic phase who are off tx are 
  ##' at set point, so it balances out for them.
  vl[selected][tx.status[selected] == 1] <- vl.supp[tx.status[selected] == 1]
  
  
  ### Last neg test before present for negatives
  selected <- which(status == 0 & tt.traj %in% c(3, 4))
  
  if (dat$param$testing.pattern == "interval") {
    tslt <- ceiling(runif(length(selected),
                          min = 0,
                          max = test.int[selected]))
  }
  if (dat$param$testing.pattern == "memoryless") {
    stop("Intertest interval parameter calculated assuming interval method. Revise parameter estimation procedure for memoryless process.")
  }
  last.neg.test[selected] <- -tslt
  

  ### Set all onto dat$attr
  dat$attr$stage <- stage
  dat$attr$stage.time <- stage.time
  dat$attr$inf.time <- inf.time
  dat$attr$vl <- vl
  dat$attr$diag.status <- diag.status
  dat$attr$diag.time <- diag.time
  dat$attr$last.neg.test <- last.neg.test
  dat$attr$tx.status <- tx.status
  dat$attr$tx.init.time <- tx.init.time
  dat$attr$cum.time.on.tx <- cum.time.on.tx
  dat$attr$cum.time.off.tx <- cum.time.off.tx
  dat$attr$infector <- infector
  dat$attr$inf.role <- inf.role
  dat$attr$inf.type <- inf.type
  dat$attr$inf.diag <- inf.diag
  dat$attr$inf.tx <- inf.tx
  dat$attr$inf.stage <- inf.stage

  return(dat)

}


#' @title Sets the CCR5 genetic status of persons for the WHAMP model
#'
#' @description Initializes the CCR5-delta-32 genetic allele of the men in the
#'              population, based on parameters defining the probability
#'              distribution.
#'
#' @param dat Data object created in initialization module.
#'
#' @export
#' @keywords initiation utility msm
#'
init_ccr5_msm_whamp <- function(dat) {

  num.H..wa <- dat$init$num.H..wa
  num.B..wa <- dat$init$num.B..wa
  num.O..wa <- dat$init$num.O..wa
  num <- num.H..wa + num.B..wa + num.O..wa

  ids.H..wa <- which(dat$attr$race..wa == "H")
  ids.B..wa <- which(dat$attr$race..wa == "B")
  ids.O..wa <- which(dat$attr$race..wa == "O")
  
  race..wa <- dat$attr$race..wa
  status <- dat$attr$status

  nInfB <- sum(race..wa == "B" & status == 1)
  nInfH <- sum(race..wa == "H" & status == 1)
  nInfO <- sum(race..wa == "O" & status == 1)
  
  ##  CCR5 genotype
  ccr5.heteroz.rr <- dat$param$ccr5.heteroz.rr
  ccr5 <- rep("WW", num)

  ### Black MSM
  # homozygotes for deletion
  num.ccr5.DD.B <- dat$param$ccr5.B.prob[1] * num.B..wa
  # heterozygotes
  num.ccr5.DW.B <- dat$param$ccr5.B.prob[2] * num.B..wa
  # homozygotes for deletion
  num.ccr5.WW.B <- num.B..wa - num.ccr5.DD.B - num.ccr5.DW.B
  # DD's can't be infected
  num.uninf.ccr5.DD.B <- round(num.ccr5.DD.B)
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.B <- round(num.ccr5.DW.B * nInfB * ccr5.heteroz.rr /
                             (num.ccr5.WW.B + num.ccr5.DW.B * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.B <- round(num.ccr5.DW.B - num.inf.ccr5.DW.B)
  inf.B <- which(status == 1 & race..wa == "B")
  inf.ccr5.DW.B <- sample(inf.B, num.inf.ccr5.DW.B, replace = FALSE)
  ccr5[inf.ccr5.DW.B] <- "DW"
  uninf.B <- which(status == 0 & race..wa == "B")
  uninf.ccr5.DWDD.B <- sample(uninf.B, num.uninf.ccr5.DW.B + num.uninf.ccr5.DD.B, replace = FALSE)
  uninf.ccr5.DW.B <- sample(uninf.ccr5.DWDD.B, num.uninf.ccr5.DW.B)
  uninf.ccr5.DD.B <- setdiff(uninf.ccr5.DWDD.B, uninf.ccr5.DW.B)
  ccr5[uninf.ccr5.DW.B] <- "DW"
  ccr5[uninf.ccr5.DD.B] <- "DD"
  
  ### Hispanic MSM
  # homozygotes for deletion
  num.ccr5.DD.H <- dat$param$ccr5.H.prob[1] * num.H..wa
  # heterozygotes
  num.ccr5.DW.H <- dat$param$ccr5.H.prob[2] * num.H..wa
  # homozygotes for deletion
  num.ccr5.WW.H <- num.H..wa - num.ccr5.DD.H - num.ccr5.DW.H
  # DD's can't be infected
  num.uninf.ccr5.DD.H <- round(num.ccr5.DD.H)
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.H <- round(num.ccr5.DW.H * nInfH * ccr5.heteroz.rr /
                               (num.ccr5.WW.H + num.ccr5.DW.H * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.H <- round(num.ccr5.DW.H - num.inf.ccr5.DW.H)
  inf.H <- which(status == 1 & race..wa == "H")
  inf.ccr5.DW.H <- sample(inf.H, num.inf.ccr5.DW.H, replace = FALSE)
  ccr5[inf.ccr5.DW.H] <- "DW"
  uninf.H <- which(status == 0 & race..wa == "H")
  uninf.ccr5.DWDD.H <- sample(uninf.H, num.uninf.ccr5.DW.H + num.uninf.ccr5.DD.H, replace = FALSE)
  uninf.ccr5.DW.H <- sample(uninf.ccr5.DWDD.H, num.uninf.ccr5.DW.H)
  uninf.ccr5.DD.H <- setdiff(uninf.ccr5.DWDD.H, uninf.ccr5.DW.H)
  ccr5[uninf.ccr5.DW.H] <- "DW"
  ccr5[uninf.ccr5.DD.H] <- "DD"
  
  ### Other MSM
  # homozygotes for deletion
  num.ccr5.DD.O <- dat$param$ccr5.O.prob[1] * num.O..wa
  # heterozygotes
  num.ccr5.DW.O <- dat$param$ccr5.O.prob[2] * num.O..wa
  # homozygotes for deletion
  num.ccr5.WW.O <- num.O..wa - num.ccr5.DD.O - num.ccr5.DW.O
  # DD's can't be infected
  num.uninf.ccr5.DD.O <- round(num.ccr5.DD.O)
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.O <- round(num.ccr5.DW.O * nInfO * ccr5.heteroz.rr /
                               (num.ccr5.WW.O + num.ccr5.DW.O * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.O <- round(num.ccr5.DW.O - num.inf.ccr5.DW.O)
  inf.O <- which(status == 1 & race..wa == "O")
  inf.ccr5.DW.O <- sample(inf.O, num.inf.ccr5.DW.O, replace = FALSE)
  ccr5[inf.ccr5.DW.O] <- "DW"
  uninf.O <- which(status == 0 & race..wa == "O")
  uninf.ccr5.DWDD.O <- sample(uninf.O, num.uninf.ccr5.DW.O + num.uninf.ccr5.DD.O, replace = FALSE)
  uninf.ccr5.DW.O <- sample(uninf.ccr5.DWDD.O, num.uninf.ccr5.DW.O)
  uninf.ccr5.DD.O <- setdiff(uninf.ccr5.DWDD.O, uninf.ccr5.DW.O)
  ccr5[uninf.ccr5.DW.O] <- "DW"
  ccr5[uninf.ccr5.DD.O] <- "DD"

  dat$attr$ccr5 <- ccr5

  return(dat)
}


#' @title Re-Initialization Module
#'
#' @description This function reinitializes an epidemic model to restart at a
#'              specified time step given an input \code{netsim} object.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netsim}}.
#' @inheritParams initialize_msm
#'
#' @details
#' Currently, the necessary components that must be on \code{x} for a simulation
#' to be restarted must be: param, control, nwparam, epi, attr, temp, el, p.
#' TODO: describe this more.
#'
#' @return
#' This function resets the data elements on the \code{dat} master data object
#' in the needed ways for the time loop to function.
#'
#' @export
#' @keywords module msm
#'
reinit_msm <- function(x, param, init, control, s) {

  need.for.reinit <- c("param", "control", "nwparam", "epi", "attr", "temp", "el", "p")
  if (!all(need.for.reinit %in% names(x))) {
    stop("x must contain the following elements for restarting: ",
         "param, control, nwparam, epi, attr, temp, el, p",
         call. = FALSE)
  }

  if (length(x$el) == 1) {
    s <- 1
  }

  dat <- list()

  dat$param <- param
  dat$param$modes <- 1
  dat$control <- control
  dat$nwparam <- x$nwparam

  dat$epi <- sapply(x$epi, function(var) var[s])
  names(dat$epi) <- names(x$epi)

  dat$el <- x$el[[s]]
  dat$p <- x$p[[s]]

  dat$attr <- x$attr[[s]]

  if (!is.null(x$stats)) {
    dat$stats <- list()
    if (!is.null(x$stats$nwstats)) {
      dat$stats$nwstats <- x$stats$nwstats[[s]]
    }
  }

  dat$temp <- x$temp[[s]]

  class(dat) <- "dat"
  return(dat)
}
