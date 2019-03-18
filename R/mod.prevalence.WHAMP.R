
#' @title Prevalence Calculations within Time Steps
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation for the WHAMP model
#'
#' @inheritParams aging_msm
#'
#' @details
#' Summary statistic calculations are of two broad forms: prevalence and
#' incidence. This function establishes the summary statistic vectors for both
#' prevalence and incidence at time 1, and then calculates the prevalence
#' statistics for times 2 onward. Incidence statistics (e.g., number of new
#' infections or deaths) are calculated within the modules as they depend on
#' vectors that are not stored external to the module.
#'
#' @return
#' This function returns the \code{dat} object with an updated summary of current
#' attributes stored in \code{dat$epi}.
#'
#' @keywords module msm
#'
#' @export
#'
prevalence_msm_whamp <- function(dat, at) {

  # Attributes

  active <- dat$attr$active
  race..wa <- dat$attr$race..wa
  region <- dat$attr$region
  tt.traj <- dat$attr$tt.traj
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  tx.status <- dat$attr$tx.status
  vl <- dat$attr$vl
  prepStat <- dat$attr$prepStat
  prepElig <- dat$attr$prepElig
  # rGC <- dat$attr$rGC
  # uGC <- dat$attr$uGC
  # rCT <- dat$attr$rCT
  # uCT <- dat$attr$uCT
  # rGC.sympt <- dat$attr$rGC.sympt
  # uGC.sympt <- dat$attr$uGC.sympt
  # rCT.sympt <- dat$attr$rCT.sympt
  # uCT.sympt <- dat$attr$uCT.sympt

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  if (at == 1) {
    dat$epi$num <- rNA
    dat$epi$num.H..wa <- rNA
    dat$epi$num.B..wa <- rNA
    dat$epi$num.O..wa <- rNA
    dat$epi$num.KC <- rNA
    dat$epi$num.OW <- rNA
    dat$epi$num.EW <- rNA
    
    dat$epi$dth.neg..wa <- rNA
    dat$epi$dth.pos..wa <- rNA
    dat$epi$dth.age <- rNA
    dat$epi$dth.pos.bkg <- rNA
    dat$epi$nBirths.B..wa <- rNA
    dat$epi$nBirths.H..wa <- rNA
    dat$epi$nBirths.O..wa <- rNA
    
    dat$epi$s.num <- rNA
    dat$epi$i.num <- rNA
    dat$epi$i.num.H..wa <- rNA
    dat$epi$i.num.B..wa <- rNA
    dat$epi$i.num.O..wa <- rNA
    dat$epi$i.num.KC <- rNA
    dat$epi$i.num.OW <- rNA
    dat$epi$i.num.EW <- rNA
    dat$epi$i.prev <- rNA
    dat$epi$i.prev.H..wa <- rNA
    dat$epi$i.prev.B..wa <- rNA
    dat$epi$i.prev.O..wa <- rNA
    dat$epi$i.prev.KC <- rNA
    dat$epi$i.prev.OW <- rNA
    dat$epi$i.prev.EW <- rNA
    dat$epi$incid <- rNA
    dat$epi$ir100 <- rNA
    dat$epi$prev.dx <- rNA
    dat$epi$prev.dx.H <- rNA
    dat$epi$prev.dx.B <- rNA
    dat$epi$prev.dx.O <- rNA
    dat$epi$prev.dx.KC <- rNA
    dat$epi$prev.dx.OW <- rNA
    dat$epi$prev.dx.EW <- rNA
    
    dat$epi$prepCurr <- rNA
    dat$epi$prepCov <- rNA
    dat$epi$prepElig <- rNA
    dat$epi$prepStart <- rNA
    dat$epi$i.num.prep0 <- rNA
    dat$epi$i.num.prep1 <- rNA
    dat$epi$i.prev.prep0 <- rNA
    dat$epi$i.prev.prep1 <- rNA
    
    dat$epi$tt.traj.1 <- rNA
    dat$epi$tt.traj.2 <- rNA
    dat$epi$tt.traj.3 <- rNA
    dat$epi$tt.traj.4 <- rNA
    dat$epi$undx <- rNA
    dat$epi$ontx <- rNA
    dat$epi$vs <- rNA

    # dat$epi$prev.rgc <- rNA
    # dat$epi$prev.ugc <- rNA
    # dat$epi$prev.gc <- rNA
    # dat$epi$prev.gc.sympt <- rNA
    # dat$epi$prev.gc.dual <- rNA
    # 
    # dat$epi$prev.rct <- rNA
    # dat$epi$prev.uct <- rNA
    # dat$epi$prev.ct <- rNA
    # dat$epi$prev.ct.sympt <- rNA
    # dat$epi$prev.ct.dual <- rNA
    # 
    # dat$epi$prev.rgcct <- rNA
    # dat$epi$prev.ugcct <- rNA
    # 
    # dat$epi$incid.rgc <- rNA
    # dat$epi$incid.ugc <- rNA
    # dat$epi$incid.gc <- rNA
    # dat$epi$incid.rct <- rNA
    # dat$epi$incid.uct <- rNA
    # dat$epi$incid.ct <- rNA
    # 
    # dat$epi$ir100.rgc <- rNA
    # dat$epi$ir100.ugc <- rNA
    # dat$epi$ir100.gc <- rNA
    # dat$epi$ir100.rct <- rNA
    # dat$epi$ir100.uct <- rNA
    # dat$epi$ir100.ct <- rNA
    # 
    # dat$epi$ir100.sti <- rNA
    # dat$epi$incid.gcct.prep <- rNA
    # 
    # dat$epi$recov.rgc <- rNA
    # dat$epi$recov.ugc <- rNA
    # dat$epi$recov.rct <- rNA
    # dat$epi$recov.uct <- rNA

    dat$epi$trans.main <- rNA
    dat$epi$trans.casl <- rNA
    dat$epi$trans.inst <- rNA

    # dat$epi$txGC <- rNA
    # dat$epi$txCT <- rNA
  }

  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)
  dat$epi$num.H..wa[at] <- sum(race..wa == "H", na.rm = TRUE)
  dat$epi$num.B..wa[at] <- sum(race..wa == "B", na.rm = TRUE)
  dat$epi$num.O..wa[at] <- sum(race..wa == "O", na.rm = TRUE)
  dat$epi$num.KC[at] <- sum(region == "KC", na.rm = TRUE)
  dat$epi$num.OW[at] <- sum(region == "OW", na.rm = TRUE)
  dat$epi$num.EW[at] <- sum(region == "EW", na.rm = TRUE)
  
  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)
  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$i.num.H..wa[at] <- sum(status == 1 & race..wa == "H", na.rm = TRUE)
  dat$epi$i.num.B..wa[at] <- sum(status == 1 & race..wa == "B", na.rm = TRUE)
  dat$epi$i.num.O..wa[at] <- sum(status == 1 & race..wa == "O", na.rm = TRUE)
  dat$epi$i.num.KC[at] <- sum(status == 1 & region == "KC", na.rm = TRUE)
  dat$epi$i.num.OW[at] <- sum(status == 1 & region == "OW", na.rm = TRUE)
  dat$epi$i.num.EW[at] <- sum(status == 1 & region == "EW", na.rm = TRUE)
  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  dat$epi$i.prev.H..wa[at] <- dat$epi$i.num.H..wa[at] / dat$epi$num.H..wa[at]
  dat$epi$i.prev.B..wa[at] <- dat$epi$i.num.B..wa[at] / dat$epi$num.B..wa[at]
  dat$epi$i.prev.O..wa[at] <- dat$epi$i.num.O..wa[at] / dat$epi$num.O..wa[at]
  dat$epi$i.prev.KC[at] <- dat$epi$i.num.KC[at] / dat$epi$num.KC[at]
  dat$epi$i.prev.OW[at] <- dat$epi$i.num.OW[at] / dat$epi$num.OW[at]
  dat$epi$i.prev.EW[at] <- dat$epi$i.num.EW[at] / dat$epi$num.EW[at]
  dat$epi$prev.dx[at] <- sum(diag.status == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.dx.H[at] <- sum(diag.status == 1 & race..wa == "H", na.rm = TRUE) / dat$epi$num.H..wa[at]
  dat$epi$prev.dx.B[at] <- sum(diag.status == 1 & race..wa == "B", na.rm = TRUE) / dat$epi$num.B..wa[at]
  dat$epi$prev.dx.O[at] <- sum(diag.status == 1 & race..wa == "O", na.rm = TRUE) / dat$epi$num.O..wa[at]
  dat$epi$prev.dx.KC[at] <- sum(diag.status == 1 & region == "KC", na.rm = TRUE) / dat$epi$num.KC[at]
  dat$epi$prev.dx.OW[at] <- sum(diag.status == 1 & region == "OW", na.rm = TRUE) / dat$epi$num.OW[at]
  dat$epi$prev.dx.EW[at] <- sum(diag.status == 1 & region == "EW", na.rm = TRUE) / dat$epi$num.EW[at]
  
  # Incidence rates per 100 per year overall and by attribute
  dat$epi$ir100[at] <- (dat$epi$incid[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid[at]) * 5200 
  dat$epi$ir100.B[at] <- (dat$epi$incid.B[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.B[at]) * 5200 
  dat$epi$ir100.H[at] <- (dat$epi$incid.H[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.H[at]) * 5200
  dat$epi$ir100.O[at] <- (dat$epi$incid.O[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.O[at]) * 5200
  dat$epi$ir100.KC[at] <- (dat$epi$incid.KC[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.KC[at]) * 5200
  dat$epi$ir100.OW[at] <- (dat$epi$incid.OW[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.OW[at]) * 5200
  dat$epi$ir100.EW[at] <- (dat$epi$incid.EW[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.EW[at]) * 5200
  dat$epi$ir100.18to24[at] <- (dat$epi$incid.18to24[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.18to24[at]) * 5200
  dat$epi$ir100.25to29[at] <- (dat$epi$incid.25to29[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.25to29[at]) * 5200
  dat$epi$ir100.30to34[at] <- (dat$epi$incid.30to34[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.30to34[at]) * 5200
  dat$epi$ir100.35to39[at] <- (dat$epi$incid.35to39[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.35to39[at]) * 5200
  dat$epi$ir100.40to44[at] <- (dat$epi$incid.40to44[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.40to44[at]) * 5200
  dat$epi$ir100.45to49[at] <- (dat$epi$incid.45to49[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.45to49[at]) * 5200
  dat$epi$ir100.50to54[at] <- (dat$epi$incid.50to54[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.50to54[at]) * 5200
  dat$epi$ir100.55to59[at] <- (dat$epi$incid.55to59[at] / (sum(status == 0, na.rm = TRUE)) + dat$epi$incid.55to59[at]) * 5200
  
  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  dat$epi$prepElig[at] <- sum(prepElig == 1, na.rm = TRUE)
  dat$epi$i.num.prep0[at] <- sum((is.na(prepStat) | prepStat == 0) & status == 1, na.rm = TRUE)
  dat$epi$i.num.prep1[at] <- sum(prepStat == 1 & status == 1, na.rm = TRUE)
  dat$epi$i.prev.prep0[at] <- dat$epi$i.num.prep0[at] /
    sum(prepStat == 0, na.rm = TRUE)
  if (at == 1) {
    dat$epi$i.prev.prep1[1] <- 0
  } else {
    dat$epi$i.prev.prep1[at] <- dat$epi$i.num.prep1[at] / sum(prepStat == 1, na.rm = TRUE)
  }
  
  dat$epi$tt.traj.1[at] <- sum(tt.traj == 1 & status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$tt.traj.2[at] <- sum(tt.traj == 2 & status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$tt.traj.3[at] <- sum(tt.traj == 3 & status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$tt.traj.4[at] <- sum(tt.traj == 4 & status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$undx[at] <- sum(status == 1 & !(diag.status == 1), na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$ontx[at] <- sum(status == 1 & tx.status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$vs[at] <- sum(status ==1 & vl <= 1.5, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)

  # dat$epi$prev.rgc[at] <- sum(rGC == 1, na.rm = TRUE) / dat$epi$num[at]
  # dat$epi$prev.ugc[at] <- sum(uGC == 1, na.rm = TRUE) / dat$epi$num[at]
  # dat$epi$prev.gc[at] <- sum((rGC == 1 | uGC == 1), na.rm = TRUE) / dat$epi$num[at]
  # dat$epi$prev.gc.sympt[at] <- sum((rGC.sympt == 1 | uGC.sympt == 1)) / dat$epi$num[at]
  # dat$epi$prev.gc.dual[at] <- sum((rGC == 1 & uGC == 1), na.rm = TRUE) / dat$epi$num[at]
  # 
  # dat$epi$prev.rct[at] <- sum(rCT == 1, na.rm = TRUE) / dat$epi$num[at]
  # dat$epi$prev.uct[at] <- sum(uCT == 1, na.rm = TRUE) / dat$epi$num[at]
  # dat$epi$prev.ct[at] <- sum((rCT == 1 | uCT == 1), na.rm = TRUE) / dat$epi$num[at]
  # dat$epi$prev.ct.sympt[at] <- sum((rCT.sympt == 1 | uCT.sympt == 1)) / dat$epi$num[at]
  # dat$epi$prev.ct.dual[at] <- sum((rCT == 1 & uCT == 1), na.rm = TRUE) / dat$epi$num[at]
  # 
  # dat$epi$prev.rgcct[at] <- sum(rGC == 1 | rCT == 1, na.rm = TRUE) / dat$epi$num[at]
  # dat$epi$prev.ugcct[at] <- sum(uGC == 1 | uCT == 1, na.rm = TRUE) / dat$epi$num[at]
  # 
  # dat$epi$ir100.rgc[at] <- (dat$epi$incid.rgc[at] / sum(rGC == 0, na.rm = TRUE)) * 5200
  # dat$epi$ir100.ugc[at] <- (dat$epi$incid.ugc[at] / sum(uGC == 0, na.rm = TRUE)) * 5200
  # dat$epi$ir100.gc[at] <- (dat$epi$incid.gc[at] /
  #                            (sum(rGC == 0, na.rm = TRUE) +
  #                               sum(uGC == 0, na.rm = TRUE))) * 5200
  # 
  # dat$epi$ir100.rct[at] <- (dat$epi$incid.rct[at] / sum(rCT == 0, na.rm = TRUE)) * 5200
  # dat$epi$ir100.uct[at] <- (dat$epi$incid.uct[at] / sum(uCT == 0, na.rm = TRUE)) * 5200
  # dat$epi$ir100.ct[at] <- (dat$epi$incid.ct[at] /
  #                            (sum(rCT == 0, na.rm = TRUE) +
  #                               sum(uCT == 0, na.rm = TRUE))) * 5200
  # 
  # dat$epi$prev.sti[at] <- sum(rGC == 1 | uGC == 1 |
  #                               rCT ==1 | uCT == 1, na.rm = TRUE) / dat$epi$num[at]
  # dat$epi$ir100.sti[at] <- ((dat$epi$incid.ct[at] + dat$epi$incid.gc[at]) /
  #                             (sum(rGC == 0, na.rm = TRUE) +
  #                                sum(uGC == 0, na.rm = TRUE) +
  #                                sum(rCT == 0, na.rm = TRUE) +
  #                                sum(uCT == 0, na.rm = TRUE))) * 5200
  # 
  # dat$epi$ir100.sti.prep[at] <- (dat$epi$incid.gcct.prep[at] /
  #                                 (sum(rGC == 0 & prepStat == 1, na.rm = TRUE) +
  #                                  sum(uGC == 0 & prepStat == 1, na.rm = TRUE) +
  #                                  sum(rCT == 0 & prepStat == 1, na.rm = TRUE) +
  #                                  sum(uCT == 0 & prepStat == 1, na.rm = TRUE))) * 5200

  return(dat)
}


#' @export
#' @rdname prevalence_msm
prevalence_het <- function(dat, at) {

  status <- dat$attr$status
  male <- dat$attr$male
  age <- dat$attr$age

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  # Initialize vectors
  if (at == 1) {
    dat$epi$i.num <- rNA
    dat$epi$num <- rNA

    dat$epi$i.num.male <- rNA
    dat$epi$i.num.feml <- rNA
    dat$epi$i.prev.male <- rNA
    dat$epi$i.prev.feml <- rNA

    dat$epi$num.male <- rNA
    dat$epi$num.feml <- rNA
    dat$epi$meanAge <- rNA
    dat$epi$propMale <- rNA

    dat$epi$si.flow <- rNA
    dat$epi$si.flow.male <- rNA
    dat$epi$si.flow.feml <- rNA

    dat$epi$b.flow <- rNA
    dat$epi$ds.flow <- dat$epi$di.flow <- rNA
  }

  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$num[at] <- length(status)

  dat$epi$i.num.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE)
  dat$epi$i.num.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE)
  dat$epi$i.prev.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE) /
    sum(male == 1, na.rm = TRUE)
  dat$epi$i.prev.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE) /
    sum(male == 0, na.rm = TRUE)

  dat$epi$num.male[at] <- sum(male == 1, na.rm = TRUE)
  dat$epi$num.feml[at] <- sum(male == 0, na.rm = TRUE)
  dat$epi$meanAge[at] <- mean(age, na.rm = TRUE)
  dat$epi$propMale[at] <- mean(male, na.rm = TRUE)

  return(dat)
}


whichVlSupp <- function(attr, param) {
  which(attr$status == 1 &
        attr$vlLevel <= log10(50) &
        (attr$age - attr$ageInf) * (365 / param$time.unit) >
        (param$vl.acute.topeak + param$vl.acute.toset))
}
