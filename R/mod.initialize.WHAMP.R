
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

    ##--ORIGINAL - DELETE AFTER FINISH DEBUGGING
    dat$attr$race <- get.vertex.attribute(nw[[1]], "race")
    num.B <- dat$init$num.B
    num.W <- dat$init$num.W
    ids.B <- which(dat$attr$race == "B")
    ids.W <- which(dat$attr$race == "W")
    
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
  p1 <- dat$param$cond.pers.always.prob
  p2 <- dat$param$cond.inst.always.prob
  rho <- dat$param$cond.always.prob.corr
  uai.always <- bindata::rmvbin(num, c(p1, p2), bincorr = (1 - rho) * diag(2) + rho)
  dat$attr$cond.always.pers <- uai.always[, 1]
  dat$attr$cond.always.inst <- uai.always[, 2]

  # Arrival and departure
  dat$attr$arrival.time <- rep(1, num)

  # Circumcision
  circ <- rep(NA, num)
  circ[ids.H..wa] <- sample(apportion_lr(num.H..wa, 0:1, 1 - param$circ.H.prob))
  circ[ids.B..wa] <- sample(apportion_lr(ids.B..wa, 0:1, 1 - param$circ.B.prob))
  circ[ids.O..wa] <- sample(apportion_lr(ids.O..wa, 0:1, 1 - param$circ.O.prob))
  dat$attr$circ <- circ

  # PrEP Attributes
  dat$attr$prepClass <- rep(NA, num)
  dat$attr$prepElig <- rep(NA, num)
  dat$attr$prepStat <- rep(0, num)
  dat$attr$prepStartTime <- rep(NA, num)
  dat$attr$prepLastRisk <- rep(NA, num)
  dat$attr$prepLastStiScreen <- rep(NA, num)

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
#'              infected persons.
#'
#' @param dat Data object created in initialization module.
#'
#' @export
#' @keywords initiation utility msm
#'
init_status_msm_whamp <- function(dat) {

  age <- dat$attr$age
  race <- dat$attr$race #-- Delete this old code when finish debugging
  race..wa <- dat$attr$race..wa
  region <- dat$attr$region
  
  ids.B <- which(dat$attr$race == "B")
  ids.W <- which(dat$attr$race == "W")
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
  
  num.B <- dat$init$num.B #-- Delete this old code eventually 
  num.W <- dat$init$num.W #-- Delete this old code eventually 
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
    nInfB <- round(dat$init$prev.B * num.B) #-- Delete this old code when finish debugging
    nInfW <- round(dat$init$prev.W * num.W) #-- Delete this old code when finish debugging
  
  nInfH..wa <- round(dat$init$prev.H..wa * num.H..wa)
  nInfB..wa <- round(dat$init$prev.B..wa * num.B..wa)
  nInfO..wa <- round(dat$init$prev.O..wa * num.O..wa)

  # Age-based infection probability
    probInfCrB <- age[ids.B] * dat$init$init.prev.age.slope.B #-- Delete this old code when finish debugging
    probInfB <- probInfCrB + (nInfB - sum(probInfCrB)) / num.B #-- Delete this old code when finish debugging

    probInfCrW <- age[ids.W] * dat$init$init.prev.age.slope.W #-- Delete this old code when finish debugging
    probInfW <- probInfCrW + (nInfW - sum(probInfCrW)) / num.W #-- Delete this old code when finish debugging

    if (any(probInfB <= 0) | any(probInfW <= 0)) {  #-- Delete this old code when finish debugging
      stop("Slope of initial prevalence by age must be sufficiently low to ",
           "avoid non-positive probabilities.", call. = FALSE)
    }

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
  avg.test.int <- dat$param$iti.coefs[1] + centered.age * dat$param$iti.coefs[2] + centered.age^2 * dat$param$iti.coefs[3]

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

  time.sex.active <- pmax(1,
                          round((365 / dat$param$time.unit) * age - (365 / dat$param$time.unit) *
                                  min(dat$init$ages), 0))

  vlar.int <- dat$param$vl.acute.rise.int
  vlap <- dat$param$vl.acute.peak
  vlaf.int <- dat$param$vl.acute.fall.int
  vlsp <- dat$param$vl.set.point
  vldo.int <- dat$param$vl.aids.onset.int
  vl.aids.int <- dat$param$vl.aids.int
  vlf  <- dat$param$vl.fatal
  vlds <- (vlf - vlsp) / vl.aids.int
  vl.acute.int <- vlar.int + vlaf.int


  ### Non-treater type: tester and non-tester
  selected <- which(status == 1 & tt.traj %in% c(1, 2))
  max.inf.time <- pmin(time.sex.active[selected], vldo.int + vl.aids.int)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  tx.status[selected] <- 0
  cum.time.on.tx[selected] <- 0
  cum.time.off.tx[selected] <- time.since.inf

  stage[selected[time.since.inf <= vlar.int]] <- 1
  stage[selected[time.since.inf > vlar.int & time.since.inf <= vl.acute.int]] <- 2
  stage[selected[time.since.inf > vl.acute.int & time.since.inf <= vldo.int]] <- 3
  stage[selected[time.since.inf > vldo.int]] <- 4

  stage.time[selected][stage[selected] == 1] <- time.since.inf[stage[selected] == 1]
  stage.time[selected][stage[selected] == 2] <- time.since.inf[stage[selected] == 2] -
                                                   vlar.int
  stage.time[selected][stage[selected] == 3] <- time.since.inf[stage[selected] == 3] -
                                                  vl.acute.int
  stage.time[selected][stage[selected] == 4] <- time.since.inf[stage[selected] == 4] -
                                                  vldo.int

  
  # Assumes a linear rate of change in VL up to peak viremia in acute phase and from peak down to set point
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) * (time.since.inf <= vldo.int) * (vlsp) +
                  (time.since.inf > vldo.int) * (vlsp + (time.since.inf - vldo.int) * vlds)

  selected <- which(status == 1 & tt.traj == 1)
  diag.status[selected] <- 0

  selected <- which(status == 1 & tt.traj == 2)

  # Time to next test
  if (dat$param$testing.pattern == "interval") {
    ttntest <- ceiling(runif(length(selected),
                             min = 0,
                             max = avg.test.int))
  }
  if (dat$param$testing.pattern == "memoryless") {
    stop("Intertest interval parameter calculated assuming interval method. Revise parameter estimation procedure for memoryless process.")
  }

  twind.int <- dat$param$test.window.int
  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
                           -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]

  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1


  ### Full adherent type

  # Create set of expected values for (cum.time.off.tx, cum.time.on.tx)
  tx.init.time.B <- twind.int + avg.test.int + 1 / dat$param$tx.init.B.prob

  # Stage for Blacks
  prop.time.on.tx.B <- dat$param$tx.reinit.B.prob /
                       (dat$param$tx.halt.B.prob + dat$param$tx.reinit.B.prob)
  offon.B <- matrix(c(1:tx.init.time.B, rep(0, tx.init.time.B)),
                    nrow = tx.init.time.B)
  numsteps.B <- (dat$param$max.time.off.tx.full.int - tx.init.time.B) /
                (1 - prop.time.on.tx.B)
  offon.B <- rbind(offon.B,
                   cbind(tx.init.time.B + (1 - prop.time.on.tx.B) * 1:numsteps.B,
                         prop.time.on.tx.B * 1:numsteps.B))
  offon.B <- round(offon.B)
  exp.dur.chronic.B <- nrow(offon.B) - vl.acute.int
  exp.onset.aids.B <- nrow(offon.B)
  offon.last.B <- offon.B[nrow(offon.B), ]
  offon.B <- rbind(offon.B,
                   matrix(c(offon.last.B[1] + (1:vl.aids.int),
                            rep(offon.last.B[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.B <- nrow(offon.B)
  offon.B[, 2] <- (1:max.possible.inf.time.B) - offon.B[, 1]
  stage.B <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B, vl.aids.int))
  stage.time.B <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B, 1:vl.aids.int)

  # Stage for Whites
  prop.time.on.tx.W <- dat$param$tx.reinit.W.prob /
    (dat$param$tx.halt.W.prob + dat$param$tx.reinit.W.prob)
  offon.W <- matrix(c(1:tx.init.time.W, rep(0, tx.init.time.W)),
                    nrow = tx.init.time.W)
  numsteps.W <- (dat$param$max.time.off.tx.full.int - tx.init.time.W) /
    (1 - prop.time.on.tx.W)
  offon.W <- rbind(offon.W,
                   cbind(tx.init.time.W + (1 - prop.time.on.tx.W) * 1:numsteps.W,
                         prop.time.on.tx.W * 1:numsteps.W))
  offon.W <- round(offon.W)
  exp.dur.chronic.W <- nrow(offon.W) - vl.acute.int
  exp.onset.aids.W <- nrow(offon.W)
  offon.last.W <- offon.W[nrow(offon.W), ]
  offon.W <- rbind(offon.W,
                   matrix(c(offon.last.W[1] + (1:vl.aids.int),
                            rep(offon.last.W[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.W <- nrow(offon.W)
  offon.W[, 2] <- (1:max.possible.inf.time.W) - offon.W[, 1]
  stage.W <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W, vl.aids.int))
  stage.time.W <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W, 1:vl.aids.int)

  # Vl for Blacks
  selected <- which(status == 1 & tt.traj == 4 & race == "B")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B[time.since.inf, 1]
  stage[selected] <- stage.B[time.since.inf]
  stage.time[selected] <- stage.time.B[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.B) * (vlsp) +
                  (time.since.inf > exp.onset.aids.B) *
                  (vlsp + (time.since.inf - exp.onset.aids.B) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp

  # VL for Whites
  selected <- which(status == 1 & tt.traj == 4 & race == "W")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W[time.since.inf, 1]
  stage[selected] <- stage.W[time.since.inf]
  stage.time[selected] <- stage.time.W[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.W) * (vlsp) +
                  (time.since.inf > exp.onset.aids.W) *
                  (vlsp + (time.since.inf - exp.onset.aids.W) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp

  # Diagnosis
  selected <- which(status == 1 & tt.traj == 4)
  
  if (dat$param$testing.pattern == "interval") {
    ttntest <- ceiling(runif(length(selected),
                             min = 0,
                             max = avg.test.int))
  }
  if (dat$param$testing.pattern == "memoryless") {
    stop("Intertest interval parameter calculated assuming interval method. Revise parameter estimation procedure for memoryless process.")
  }

  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
                           -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]
  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1
  diag.status[selected][cum.time.on.tx[selected] > 0] <- 1
  last.neg.test[selected][cum.time.on.tx[selected] > 0] <- NA


  ### Part adherent type

  # Create set of expected values for (cum.time.off.tx,cum.time.on.tx)

  prop.time.on.tx.B <- dat$param$tx.reinit.B.prob /
                       (dat$param$tx.halt.B.prob + dat$param$tx.reinit.B.prob)
  offon.B <- matrix(c(1:tx.init.time.B, rep(0, tx.init.time.B)),
                    nrow = tx.init.time.B)
  while (offon.B[nrow(offon.B), 1] / dat$param$max.time.off.tx.part.int +
         offon.B[nrow(offon.B), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.B <- rbind(offon.B,
                     offon.B[nrow(offon.B), ] + c(1 - prop.time.on.tx.B,
                                                      prop.time.on.tx.B))
  }
  offon.B <- round(offon.B)
  exp.dur.chronic.B <- nrow(offon.B) - vl.acute.int
  exp.onset.aids.B <- nrow(offon.B)
  offon.last.B <- offon.B[nrow(offon.B), ]
  offon.B <- rbind(offon.B,
                   matrix(c(offon.last.B[1] + (1:vl.aids.int),
                            rep(offon.last.B[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.B <- nrow(offon.B)
  offon.B[, 2] <- (1:max.possible.inf.time.B) - offon.B[, 1]
  stage.B <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B, vl.aids.int))
  stage.time.B <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B, 1:vl.aids.int)

  prop.time.on.tx.W <- dat$param$tx.reinit.W.prob /
                       (dat$param$tx.halt.W.prob + dat$param$tx.reinit.W.prob)
  offon.W <- matrix(c(1:tx.init.time.W, rep(0, tx.init.time.W)),
                    nrow = tx.init.time.W)

  while (offon.W[nrow(offon.W), 1] / dat$param$max.time.off.tx.part.int +
         offon.W[nrow(offon.W), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.W <- rbind(offon.W,
                     offon.W[nrow(offon.W), ] + c(1 - prop.time.on.tx.W,
                                                  prop.time.on.tx.W))
  }
  offon.W <- round(offon.W)
  exp.dur.chronic.W <- nrow(offon.W) - vl.acute.int
  exp.onset.aids.W <- nrow(offon.W)
  offon.last.W <- offon.W[nrow(offon.W), ]
  offon.W <- rbind(offon.W,
                   matrix(c(offon.last.W[1] + (1:vl.aids.int),
                            rep(offon.last.W[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.W <- nrow(offon.W)
  offon.W[, 2] <- (1:max.possible.inf.time.W) - offon.W[, 1]
  stage.W <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W, vl.aids.int))
  stage.time.W <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W, 1:vl.aids.int)

  # VL for Blacks
  selected <- which(status == 1 & tt.traj == 3 & race == "B")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B[time.since.inf, 1]
  stage[selected] <- stage.B[time.since.inf]
  stage.time[selected] <- stage.time.B[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.B) * (vlsp) +
                  (time.since.inf > exp.onset.aids.B) *
                  (vlsp + (time.since.inf - exp.onset.aids.B) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp

  # VL for Whites
  selected <- which(status == 1 & tt.traj == 3 & race == "W")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W[time.since.inf, 1]
  stage[selected] <- stage.W[time.since.inf]
  stage.time[selected] <- stage.time.W[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.W) * (vlsp) +
                  (time.since.inf > exp.onset.aids.W) *
                  (vlsp + (time.since.inf - exp.onset.aids.W) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp

  # Implement diagnosis for both
  selected <- which(status == 1 & tt.traj == 3)
  if (dat$param$testing.pattern == "interval") {
    ttntest <- ceiling(runif(length(selected),
                             min = 0,
                             max = avg.test.int))
  }

  if (dat$param$testing.pattern == "memoryless") {
    stop("Intertest interval parameter calculated assuming interval method. Revise parameter estimation procedure for memoryless process.")
  }


  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
    -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]

  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1
  diag.status[selected][cum.time.on.tx[selected] > 0] <- 1
  last.neg.test[selected][cum.time.on.tx[selected] > 0] <- NA


  # Last neg test before present for negatives
  selected <- which(status == 0 & tt.traj %in% c(2, 3, 4))

  if (dat$param$testing.pattern == "interval") {
    tslt <- ceiling(runif(length(selected),
                          min = 0,
                          max = avg.test.int))
  }
  if (dat$param$testing.pattern == "memoryless") {
    stop("Intertest interval parameter calculated assuming interval method. Revise parameter estimation procedure for memoryless process.")
  }
  last.neg.test[selected] <- -tslt


  ## Set all onto dat$attr
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

  num.B <- dat$init$num.B
  num.W <- dat$init$num.W
  num.H..wa <- dat$init$num.H..wa
  num.B..wa <- dat$init$num.B..wa
  num.O..wa <- dat$init$num.O..wa
  num <- num.H..wa + num.B..wa + num.O..wa
  
  ids.B <- which(dat$attr$race == "B")
  ids.W <- which(dat$attr$race == "W")
  ids.H..wa <- which(dat$attr$race..wa == "H")
  ids.B..wa <- which(dat$attr$race..wa == "B")
  ids.O..wa <- which(dat$attr$race..wa == "O")
  
  race <- dat$attr$race #-- Delete this code eventually
  race..wa <- dat$attr$race..wa
  region <- dat$attr$region
  status <- dat$attr$status

  nInfB <- sum(race == "B" & status == 1)
  nInfW <- sum(race == "W" & status == 1)

  ##  CCR5 genotype
  ccr5.heteroz.rr <- dat$param$ccr5.heteroz.rr
  ccr5 <- rep("WW", num)

  # homozygotes for deletion
  num.ccr5.DD.B <- dat$param$ccr5.B.prob[1] * num.B
  # heterozygotes
  num.ccr5.DW.B <- dat$param$ccr5.B.prob[2] * num.B
  # homozygotes for deletion
  num.ccr5.WW.B <- num.B - num.ccr5.DD.B - num.ccr5.DW.B
  # DD's can't be infected
  num.uninf.ccr5.DD.B <- round(num.ccr5.DD.B)
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.B <- round(num.ccr5.DW.B * nInfB * ccr5.heteroz.rr /
                             (num.ccr5.WW.B + num.ccr5.DW.B * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.B <- round(num.ccr5.DW.B - num.inf.ccr5.DW.B)
  inf.B <- which(status == 1 & race == "B")
  inf.ccr5.DW.B <- sample(inf.B, num.inf.ccr5.DW.B, replace = FALSE)
  ccr5[inf.ccr5.DW.B] <- "DW"
  uninf.B <- which(status == 0 & race == "B")
  uninf.ccr5.DWDD.B <- sample(uninf.B, num.uninf.ccr5.DW.B + num.uninf.ccr5.DD.B)
  uninf.ccr5.DW.B <- sample(uninf.ccr5.DWDD.B, num.uninf.ccr5.DW.B)
  uninf.ccr5.DD.B <- setdiff(uninf.ccr5.DWDD.B, uninf.ccr5.DW.B)
  ccr5[uninf.ccr5.DW.B] <- "DW"
  ccr5[uninf.ccr5.DD.B] <- "DD"

  num.ccr5.DD.W <- dat$param$ccr5.W.prob[1] * num.W
  num.ccr5.DW.W <- dat$param$ccr5.W.prob[2] * num.W
  num.ccr5.WW.W <- num.W - num.ccr5.DD.W - num.ccr5.DW.W
  num.uninf.ccr5.DD.W <- round(num.ccr5.DD.W)
  num.inf.ccr5.DW.W <- round(num.ccr5.DW.W * nInfW * ccr5.heteroz.rr /
                             (num.ccr5.WW.W + num.ccr5.DW.W * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.W <- round(num.ccr5.DW.W - num.inf.ccr5.DW.W)
  inf.W <- which(status == 1 & race == "W")
  inf.ccr5.DW.W <- sample(inf.W, num.inf.ccr5.DW.W)
  ccr5[inf.ccr5.DW.W] <- "DW"
  uninf.W <- which(status == 0 & race == "W")
  uninf.ccr5.DWDD.W <- sample(uninf.W, num.uninf.ccr5.DW.W + num.uninf.ccr5.DD.W)
  uninf.ccr5.DW.W <- sample(uninf.ccr5.DWDD.W, num.uninf.ccr5.DW.W)
  uninf.ccr5.DD.W <- setdiff(uninf.ccr5.DWDD.W, uninf.ccr5.DW.W)
  ccr5[uninf.ccr5.DW.W] <- "DW"
  ccr5[uninf.ccr5.DD.W] <- "DD"

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
