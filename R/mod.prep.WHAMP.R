
#' @title PrEP Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
prep_msm_whamp <- function(dat, at) {

  # Function Selection ------------------------------------------------------
  
  if (at >= dat$param$riskh.start) {
    dat <- riskhist_msm(dat, at)
  } else {
    return(dat)
  }
  
  if (at < dat$param$prep.start) {
    return(dat)
  }

  # Set attributes ----------------------------------------------------------
  
  # Attributes
  uid <- dat$attr$uid
  active <- dat$attr$active
  region <- dat$attr$region
  age <- dat$attr$age
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  lnt <- dat$attr$last.neg.test
  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  prepLastRisk <- dat$attr$prepLastRisk
  prepStartTime <- dat$attr$prepStartTime
  prepDiscont <- dat$attr$prepDiscont
  spontDisc <- dat$attr$spontDisc
  # prepLastStiScreen <- dat$attr$prepLastStiScreen

  # Parameters
  prep.risk.reassess.method <- dat$param$prep.risk.reassess.method
  prep.start.step <- dat$param$prep.start
  prep.init.rate <- dat$param$prep.init.rate
  prep.coverage.init.KC <- dat$param$prep.coverage.init.KC
  prep.coverage.init.oth <- dat$param$prep.coverage.init.oth
  prep.scaleup.KC <- dat$param$prep.scaleup.KC
  prep.scaleup.oth <- dat$param$prep.scaleup.oth
  prep.cov.max.KC <- dat$param$prep.cov.max.KC
  prep.cov.max.oth <- dat$param$prep.cov.max.oth
  prep.class.prob <- dat$param$prep.class.prob
  prep.discontinue <- dat$param$prep.discont
  prep.discont.prob <- dat$param$prep.discont.prob
  prep.start.step <- dat$param$prep.start
  
  # Set coverage quota for the current time step
  if (at == prep.start.step){
    prep.coverage.KC.18to24 <- prep.coverage.init.KC[1]
    prep.coverage.KC.25to29 <- prep.coverage.init.KC[2]
    prep.coverage.KC.30to49 <- prep.coverage.init.KC[3]
    prep.coverage.KC.50plus <- prep.coverage.init.KC[4]
    prep.coverage.oth.18to24 <- prep.coverage.init.oth[1]
    prep.coverage.oth.25to29 <- prep.coverage.init.oth[2]
    prep.coverage.oth.30to49 <- prep.coverage.init.oth[3]
    prep.coverage.oth.50plus <- prep.coverage.init.oth[4]
  }
  if (at > dat$param$prep.start.step & is.finite(prep.start.step)){
    prep.coverage.KC.18to24 <- min(prep.coverage.init.KC[1] + (at - prep.start.step)*prep.scaleup.KC[1], prep.cov.max.KC[1])
    prep.coverage.KC.25to29 <- min(prep.coverage.init.KC[2] + (at - prep.start.step)*prep.scaleup.KC[2], prep.cov.max.KC[2])
    prep.coverage.KC.30to49 <- min(prep.coverage.init.KC[3] + (at - prep.start.step)*prep.scaleup.KC[3], prep.cov.max.KC[3])
    prep.coverage.KC.50plus <- min(prep.coverage.init.KC[4] + (at - prep.start.step)*prep.scaleup.KC[4], prep.cov.max.KC[4])
    prep.coverage.oth.18to24 <- min(prep.coverage.init.oth[1] + (at - prep.start.step)*prep.scaleup.oth[1], prep.cov.max.oth[1])
    prep.coverage.oth.25to29 <- min(prep.coverage.init.oth[2] + (at - prep.start.step)*prep.scaleup.oth[2], prep.cov.max.oth[2])
    prep.coverage.oth.30to49 <- min(prep.coverage.init.oth[3] + (at - prep.start.step)*prep.scaleup.oth[3], prep.cov.max.oth[3])
    prep.coverage.oth.50plus <- min(prep.coverage.init.oth[4] + (at - prep.start.step)*prep.scaleup.oth[4], prep.cov.max.oth[4])
  }


  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(active == 1 & status == 0 & prepStat == 0 & lnt == at)

  # Core eligiblity
  ind1 <- dat$attr$prep.ind.discord.ongoing
  ind2 <- dat$attr$prep.ind.uai.risk

  twind <- at - dat$param$prep.risk.int
  idsEligStart <- intersect(which(ind1 == at | ind2 >= twind),
                            idsEligStart)

  prepElig[idsEligStart] <- 1

  
  # No longer indicated for PrEP
  idsNoIndic <- which((ind1 < at | is.na(ind1)) & 
                        (ind2 < twind | is.na(ind2)))
  
  prepElig[idsNoIndic] <- 0

  ## Discontinuation ------------------------------------------------------------------

  # With changes in risk behavior 
  if (prep.risk.reassess.method == "none") {
    idsEligStop <- NULL
  } else if (prep.risk.reassess.method == "inst") {
    idsRiskAssess <- which(active == 1 & prepStat == 1)
    prepLastRisk[idsRiskAssess] <- at
    idsEligStop <- intersect(idsNoIndic, idsRiskAssess)
  } else if (prep.risk.reassess.method == "year") {
    idsRiskAssess <- which(active == 1 & prepStat == 1 & lnt == at & (at - prepLastRisk) >= 52)
    prepLastRisk[idsRiskAssess] <- at
    idsEligStop <- intersect(idsNoIndic, idsRiskAssess)
  }
  
  
  # Spontaneous (memoryless) discontinuation
  discont.elig <- which(active == 1 & prepStat == 1 & prepDiscont == 1 & diag.status != 1)
  idsDiscont <- discont.elig[rbinom(length(discont.elig), 1,
                                    prep.discont.prob) == 1]
  
  spontDisc[idsDiscont] <- 1
  
  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & diag.status == 1)

  # Death
  idsStpDth <- which(active == 0 & prepStat == 1)

  # All discontinuation
  idsStp <- c(idsStpDx, idsStpDth, idsEligStop, idsDiscont)
  
  # Calculate mean time to discontinuation among those who stopped
  time.to.disc <- rep(NA, length(idsStp))
  time.to.disc <- at - prepStartTime[idsStp]
  
  # Reset PrEP status
  prepStat[idsStp] <- 0
  prepLastRisk[idsStp] <- NA
  prepStartTime[idsStp] <- NA
  # prepLastStiScreen[idsStp] <- NA
  

  ## Initiation ----------------------------------------------------------------

  # For each group, calculate the current coverage among eligible men
  prep.cov.curr.KC.18to24 <- sum(prepStat == 1 & region %in% "KC" & age <25, na.rm = TRUE)/
    sum(prepElig == 1 & region %in% "KC" & age <25, na.rm = TRUE)
  prep.cov.curr.KC.18to24 <- ifelse(is.nan(prep.cov.curr.KC.18to24), 0, prep.cov.curr.KC.18to24)
  prep.cov.curr.KC.25to29 <- sum(prepStat == 1 & region %in% "KC" & age >=25 & age <30, na.rm = TRUE)/
    sum(prepElig == 1 & region %in% "KC" & age >=25 & age <30, na.rm = TRUE)
  prep.cov.curr.KC.25to29 <- ifelse(is.nan(prep.cov.curr.KC.25to29), 0, prep.cov.curr.KC.25to29)
  prep.cov.curr.KC.30to49 <- sum(prepStat == 1 & region %in% "KC" & age >=30 & age <50, na.rm = TRUE)/
    sum(prepElig == 1 & region %in% "KC" & age >=30 & age <50, na.rm = TRUE)
  prep.cov.curr.KC.30to49 <- ifelse(is.nan(prep.cov.curr.KC.30to49), 0, prep.cov.curr.KC.30to49)
  prep.cov.curr.KC.50plus <- sum(prepStat == 1 & region %in% "KC" & age >=50, na.rm = TRUE)/
    sum(prepElig == 1 & region %in% "KC" & age >=50, na.rm = TRUE)
  prep.cov.curr.KC.50plus <- ifelse(is.nan(prep.cov.curr.KC.50plus), 0, prep.cov.curr.KC.50plus)
  
  prep.cov.curr.oth.18to24 <- sum(prepStat == 1 & region %in% c("OW", "EW") & age <25, na.rm = TRUE)/
    sum(prepElig == 1 & region %in% c("OW", "EW") & age <25, na.rm = TRUE)
  prep.cov.curr.oth.18to24 <- ifelse(is.nan(prep.cov.curr.oth.18to24), 0, prep.cov.curr.oth.18to24)
  prep.cov.curr.oth.25to29 <- sum(prepStat == 1 & region %in% c("OW", "EW") & age >=25 & age <30, na.rm = TRUE)/
    sum(prepElig == 1 & region %in% c("OW", "EW") & age >=25 & age <30, na.rm = TRUE)
  prep.cov.curr.oth.25to29 <- ifelse(is.nan(prep.cov.curr.oth.25to29), 0, prep.cov.curr.oth.25to29)
  prep.cov.curr.oth.30to49 <- sum(prepStat == 1 & region %in% c("OW", "EW") & age >=30 & age <50, na.rm = TRUE)/
    sum(prepElig == 1 & region %in% c("OW", "EW") & age >=30 & age <50, na.rm = TRUE)
  prep.cov.curr.oth.30to49 <- ifelse(is.nan(prep.cov.curr.oth.30to49), 0, prep.cov.curr.oth.30to49)
  prep.cov.curr.oth.50plus <- sum(prepStat == 1 & region %in% c("OW", "EW") & age >=50, na.rm = TRUE)/
    sum(prepElig == 1 & region %in% c("OW", "EW") & age >=50, na.rm = TRUE)
  prep.cov.curr.oth.50plus <- ifelse(is.nan(prep.cov.curr.oth.50plus), 0, prep.cov.curr.oth.50plus)
  
  # For each group, sample ids to start PrEP if current coverage < coverage threshold
  idsEligSt.KC.18to24 <- idsEligStart[which(region %in% c("KC") & age <25)]
  idsEligSt.KC.25to29 <- idsEligStart[which(region %in% c("KC") & age >=25 & age <30)]
  idsEligSt.KC.30to49 <- idsEligStart[which(region %in% c("KC") & age >=30 & age <50)]
  idsEligSt.KC.50plus <- idsEligStart[which(region %in% c("KC") & age >=50)]
  
  idsEligSt.oth.18to24 <- idsEligStart[which(region %in% c("OW", "EW") & age <25)]
  idsEligSt.oth.25to29 <- idsEligStart[which(region %in% c("OW", "EW") & age >=25 & age <30)]
  idsEligSt.oth.30to49 <- idsEligStart[which(region %in% c("OW", "EW") & age >=30 & age <50)]
  idsEligSt.oth.50plus <- idsEligStart[which(region %in% c("OW", "EW") & age >=50)]
  
  nStart.KC.18to24 <- pmax(0, pmin(length(idsEligSt.KC.18to24), round((prep.coverage.KC.18to24 - prep.cov.curr.KC.18to24) *
                                        sum(prepElig == 1 & region %in% c("KC") & age <25, na.rm = TRUE))))
  nStart.KC.25to29 <- pmax(0, pmin(length(idsEligSt.KC.25to29), round((prep.coverage.KC.25to29 - prep.cov.curr.KC.25to29) *
                                        sum(prepElig == 1 & region %in% c("KC") & age >=25 & age <30, na.rm = TRUE))))
  nStart.KC.30to49 <- pmax(0, pmin(length(idsEligSt.KC.30to49), round((prep.coverage.KC.30to49 - prep.cov.curr.KC.30to49) *
                                        sum(prepElig == 1 & region %in% c("KC") & age >=30 & age <50, na.rm = TRUE))))
  nStart.KC.50plus <- pmax(0, pmin(length(idsEligSt.KC.50plus), round((prep.coverage.KC.50plus - prep.cov.curr.KC.50plus) *
                                        sum(prepElig == 1 & region %in% c("KC") & age >=50, na.rm = TRUE))))
  nStart.oth.18to24 <- pmax(0, pmin(length(idsEligSt.oth.18to24), round((prep.coverage.oth.18to24 - prep.cov.curr.oth.18to24) *
                                        sum(prepElig == 1 & region %in% c("OW", "EW") & age <25, na.rm = TRUE))))
  nStart.oth.25to29 <- pmax(0, pmin(length(idsEligSt.oth.25to29), round((prep.coverage.oth.25to29 - prep.cov.curr.oth.25to29) *
                                        sum(prepElig == 1 & region %in% c("OW", "EW") & age >=25 & age <30, na.rm = TRUE))))
  nStart.oth.30to49 <- pmax(0, pmin(length(idsEligSt.oth.30to49), round((prep.coverage.oth.30to49 - prep.cov.curr.oth.30to49) *
                                        sum(prepElig == 1 & region %in% c("OW", "EW") & age >=30 & age <50, na.rm = TRUE))))
  nStart.oth.50plus <- pmax(0, pmin(length(idsEligSt.oth.50plus), round((prep.coverage.oth.50plus - prep.cov.curr.oth.50plus) *
                                        sum(prepElig == 1 & region %in% c("OW", "EW") & age >=50, na.rm = TRUE))))
  
  idsStart.KC.18to24 <- NULL
  if (nStart.KC.18to24 > 0) {
    if (prep.init.rate >= 1) {
      idsStart.KC.18to24 <- ssample(idsEligSt.KC.18to24, nStart.KC.18to24)
    } else {
      idsStart.KC.18to24 <- idsEligSt.KC.18to24[rbinom(nStart.KC.18to24, 1, prep.init.rate) == 1]
    }
  }
  idsStart.KC.25to29 <- NULL
  if (nStart.KC.25to29 > 0) {
    if (prep.init.rate >= 1) {
      idsStart.KC.25to29 <- ssample(idsEligSt.KC.25to29, nStart.KC.25to29)
    } else {
      idsStart.KC.25to29 <- idsEligSt.KC.25to29[rbinom(nStart.KC.25to29, 1, prep.init.rate) == 1]
    }
  }
  idsStart.KC.30to49 <- NULL
  if (nStart.KC.30to49 > 0) {
    if (prep.init.rate >= 1) {
      idsStart.KC.30to49 <- ssample(idsEligSt.KC.30to49, nStart.KC.30to49)
    } else {
      idsStart.KC.30to49 <- idsEligSt.KC.30to49[rbinom(nStart.KC.30to49, 1, prep.init.rate) == 1]
    }
  }
  idsStart.KC.50plus <- NULL
  if (nStart.KC.50plus > 0) {
    if (prep.init.rate >= 1) {
      idsStart.KC.50plus <- ssample(idsEligSt.KC.50plus, nStart.KC.50plus)
    } else {
      idsStart.KC.50plus <- idsEligSt.KC.50plus[rbinom(nStart.KC.50plus, 1, prep.init.rate) == 1]
    }
  }
  idsStart.oth.18to24 <- NULL
  if (nStart.oth.18to24 > 0) {
    if (prep.init.rate >= 1) {
      idsStart.oth.18to24 <- ssample(idsEligSt.oth.18to24, nStart.oth.18to24)
    } else {
      idsStart.oth.18to24 <- idsEligSt.oth.18to24[rbinom(nStart.oth.18to24, 1, prep.init.rate) == 1]
    }
  }
  idsStart.oth.25to29 <- NULL
  if (nStart.oth.25to29 > 0) {
    if (prep.init.rate >= 1) {
      idsStart.oth.25to29 <- ssample(idsEligSt.oth.25to29, nStart.oth.25to29)
    } else {
      idsStart.oth.25to29 <- idsEligSt.oth.25to29[rbinom(nStart.oth.25to29, 1, prep.init.rate) == 1]
    }
  }
  idsStart.oth.30to49 <- NULL
  if (nStart.oth.30to49 > 0) {
    if (prep.init.rate >= 1) {
      idsStart.oth.30to49 <- ssample(idsEligSt.oth.30to49, nStart.oth.30to49)
    } else {
      idsStart.oth.30to49 <- idsEligSt.oth.30to49[rbinom(nStart.oth.30to49, 1, prep.init.rate) == 1]
    }
  }
  idsStart.oth.50plus <- NULL
  if (nStart.oth.50plus > 0) {
    if (prep.init.rate >= 1) {
      idsStart.oth.50plus <- ssample(idsEligSt.oth.50plus, nStart.oth.50plus)
    } else {
      idsStart.oth.50plus <- idsEligSt.oth.50plus[rbinom(nStart.oth.50plus, 1, prep.init.rate) == 1]
    }
  }
  
  idsStart <- c(idsStart.KC.18to24, idsStart.KC.25to29, idsStart.KC.30to49, idsStart.KC.50plus, 
                idsStart.oth.18to24, idsStart.oth.25to29, idsStart.oth.30to49, idsStart.oth.50plus)
  
  # Attributes
  if (length(idsStart) > 0) {
    prepStat[idsStart] <- 1
    prepStartTime[idsStart] <- at
    prepLastRisk[idsStart] <- at
    spontDisc[idsStart] <- 0

    # PrEP class
    needPC <- which(is.na(prepClass[idsStart]))
    prepClass[idsStart[needPC]] <- sample(x = 1:3, size = length(needPC),
                                          replace = TRUE, prob = prep.class.prob)
    
    # Discontinuation group
    needgp <- which(is.na(prepDiscont[idsStart]))
    prepDiscont[idsStart[needgp]] <- rbinom(length(needgp, 1, prep.discontinue)) 
    
  }
  
  ## Discontinue if current coverage exceeds threshold due to men aging into groups with lower coverage-------
  ndrop.KC.18to24 <- pmax(0, round((prep.cov.curr.KC.18to24 - prep.coverage.KC.18to24) * 
                                     sum(prepElig == 1 & region %in% c("KC") & age <25, na.rm = TRUE)))
  ndrop.KC.25to29 <- pmax(0, round((prep.cov.curr.KC.25to29 - prep.coverage.KC.25to29) * 
                                     sum(prepElig == 1 & region %in% c("KC") & age >=25 & age <30, na.rm = TRUE)))
  ndrop.KC.30to49 <- pmax(0, round((prep.cov.curr.KC.30to49 - prep.coverage.KC.30to49) * 
                                     sum(prepElig == 1 & region %in% c("KC") & age >=30 & age <50, na.rm = TRUE)))
  ndrop.KC.50plus <- pmax(0, round((prep.cov.curr.KC.50plus - prep.coverage.KC.50plus) * 
                                     sum(prepElig == 1 & region %in% c("KC") & age >=50, na.rm = TRUE)))
  ndrop.oth.18to24 <- pmax(0, round((prep.cov.curr.oth.18to24 - prep.coverage.oth.18to24) * 
                                     sum(prepElig == 1 & region %in% c("OW", "EW") & age <25, na.rm = TRUE)))
  ndrop.oth.25to29 <- pmax(0, round((prep.cov.curr.oth.25to29 - prep.coverage.oth.25to29) * 
                                     sum(prepElig == 1 & region %in% c("OW", "EW") & age >=25 & age <30, na.rm = TRUE)))
  ndrop.oth.30to49 <- pmax(0, round((prep.cov.curr.oth.30to49 - prep.coverage.oth.30to49) * 
                                     sum(prepElig == 1 & region %in% c("OW", "EW") & age >=30 & age <50, na.rm = TRUE)))
  ndrop.oth.50plus <- pmax(0, round((prep.cov.curr.oth.50plus - prep.coverage.oth.50plus) * 
                                     sum(prepElig == 1 & region %in% c("OW", "EW") & age >=50, na.rm = TRUE)))
  
  # Select ndrop men to discontinue
  on.prep.KC.18to24 <- which(active == 1 & prepStat == 1 & region %in% c("KC") & age <25)
  idsStopQuota.KC.18to24 <- ssample(on.prep.KC.18to24, ndrop.KC.18to24)
  
  on.prep.KC.25to29 <- which(active == 1 & prepStat == 1 & region %in% c("KC") & age >=25 & age <30)
  idsStopQuota.KC.25to29 <- ssample(on.prep.KC.25to29, ndrop.KC.25to29)

  on.prep.KC.30to49 <- which(active == 1 & prepStat == 1 & region %in% c("KC") & age >=30 & age <50)
  idsStopQuota.KC.30to49 <- ssample(on.prep.KC.30to49, ndrop.KC.30to49)
  
  on.prep.KC.50plus <- which(active == 1 & prepStat == 1 & region %in% c("KC") & age >=50)
  idsStopQuota.KC.50plus <- ssample(on.prep.KC.50plus, ndrop.KC.50plus)
  
  on.prep.oth.18to24 <- which(active == 1 & prepStat == 1 & region %in% c("OW", "EW") & age <25)
  idsStopQuota.oth.18to24 <- ssample(on.prep.oth.18to24, ndrop.oth.18to24)
  
  on.prep.oth.25to29 <- which(active == 1 & prepStat == 1 & region %in% c("OW", "EW") & age >=25 & age <30)
  idsStopQuota.oth.25to29 <- ssample(on.prep.oth.25to29, ndrop.oth.25to29)
  
  on.prep.oth.30to49 <- which(active == 1 & prepStat == 1 & region %in% c("OW", "EW") & age >=30 & age <50)
  idsStopQuota.oth.30to49 <- ssample(on.prep.oth.30to49, ndrop.oth.30to49)
  
  on.prep.oth.50plus <- which(active == 1 & prepStat == 1 & region %in% c("OW", "EW") & age >=50)
  idsStopQuota.oth.50plus <- ssample(on.prep.oth.50plus, ndrop.oth.50plus)
  
  idsStopQuota <- c(idsStopQuota.KC.18to24, idsStopQuota.KC.25to29, idsStopQuota.KC.30to49, idsStopQuota.KC.50plus,
                    idsStopQuota.oth.18to24, idsStopQuota.oth.25to29, idsStopQuota.oth.30to49, idsStopQuota.oth.50plus)
  
  # Reset PrEP status
  prepStat[idsStopQuota] <- 0
  prepLastRisk[idsStopQuota] <- NA
  prepStartTime[idsStopQuota] <- NA
  
  # Calculate elapsed time on PrEP among men discontinued due to coverage quotas by age
  time.to.disc <- c(time.to.disc, c(at - prepStartTime[idsStopQuota]))
  
  ## Calculate time on PrEP among current users -----------------------
  time.since.prep.start <- rep(NA, sum(prepStat == 1))
  time.since.prep.start[prepStat == 1] <- (at - prepStartTime[prepStat == 1]) 
  
  
  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepStartTime <- prepStartTime
  dat$attr$prepClass <- prepClass
  dat$attr$prepLastRisk <- prepLastRisk
  dat$attr$prepLastStiScreen <- prepLastStiScreen
  dat$attr$prepDiscont <- prepDiscont
  dat$attr$spontDisc <- spontDisc

  # Summary Statistics
  prep.cov <- sum(prepStat == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  prep.cov.KC <- sum(prepStat == 1 & region %in% "KC", na.rm = TRUE)/
    sum(prepElig == 1 & region %in% "KC", na.rm = TRUE)
  prep.cov.oth <- sum(prepStat == 1 & region %in% c("OW", "EW"), na.rm = TRUE)/
    sum(prepElig == 1 & region %in% c("OW", "EW"), na.rm = TRUE)
  
  dat$epi$prep.cov[at] <- prep.cov
  dat$epi$prep.cov.KC[at] <- prep.cov.KC
  dat$epi$prep.cov.oth[at] <- prep.cov.oth
  dat$epi$prepStart[at] <- length(idsStart)
  dat$epi$time.on.prep <- time.on.prep # Vector with the total time on PrEP among those who discontinued in this time step
  dat$epi$time.since.prep.start <- time.since.prep.start # Vector with the elapsed time on PrEP among current users

  return(dat)
}


#' @title Risk History Sub-Module
#'
#' @description Sub-Module function to track the risk history of uninfected persons
#'              for purpose of PrEP targeting in the WHAMP model.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'

riskhist_msm_whamp <- function(dat, at) {
  
  if (at < dat$param$riskh.start) {
    return(dat)
  }
  
  ## Attributes
  n <- length(dat$attr$active)
  uid <- dat$attr$uid
  dx <- dat$attr$diag.status
  since.test <- at - dat$attr$last.neg.test
  # rGC.tx <- dat$attr$rGC.tx
  # uGC.tx <- dat$attr$uGC.tx
  # rCT.tx <- dat$attr$rCT.tx
  # uCT.tx <- dat$attr$uCT.tx
  
  ## Parameters
  time.unit <- dat$param$time.unit
  
  ## Edgelist, adds uai summation per partnership from act list
  pid <- NULL # For R CMD Check
  al <- as.data.frame(dat$temp$al)
  by_pid <- group_by(al, pid)
  uai <- summarise(by_pid, uai = sum(uai))[, 2]
  el <- as.data.frame(cbind(dat$temp$el, uai))
  
  if (max(el[, 1:2]) > n) stop("riskhist max(el) > n")
  
  # Remove concordant positive edges
  el2 <- el[el$st2 == 0, ]
  
  # Initialize attributes
  if (is.null(dat$attr$prep.ind.discord.ongoing)) {
    dat$attr$prep.ind.discord.ongoing <- rep(NA, n)
    dat$attr$prep.ind.uai.risk <- rep(NA, n)
  }
  
  ## Degree ##
  main.deg <- get_degree(dat$el[[1]])
  casl.deg <- get_degree(dat$el[[2]])
  inst.deg <- get_degree(dat$el[[3]])
  
  
  ## Preconditions ##
  
  # Any UAI
  uai.any <- unique(c(el2$p1[el2$uai > 0],
                      el2$p2[el2$uai > 0]))
  
  # Monogamous partnerships: 2-sided
  tot.deg <- main.deg + casl.deg + inst.deg
  uai.mono <- intersect(which(tot.deg == 1), uai.any)
  
  # "Negative" partnerships
  tneg <- unique(c(el2$p1[el2$st1 == 0], el2$p2[el2$st1 == 0])) ##-- should'nt this be el2$p2[el2$st2] ==0] ?
  fneg <- unique(c(el2$p1[which(dx[el2$p1] == 0)], el2$p2[which(dx[el2$p1] == 0)])) ##-- should'nt this be el2$p2[which(dx[el2$p2] ==0] ?
  all.neg <- c(tneg, fneg)
  
  ## Condition 1: ongoing positive partner who has disclosed
  discord <- el2[el2$st1 == 1 & el2$ptype %in% 1:2, ]
  
  # Disclosure
  discl.list <- dat$temp$discl.list
  disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
  delt.cdl <- uid[discord[, 1]] * 1e7 + uid[discord[, 2]]
  discl <- (delt.cdl %in% disclose.cdl)
  ai.sd <- discord$p2[discl == TRUE]
  
  dat$attr$prep.ind.discord.ongoing[ai.sd] <- at
  
  ## Condition 2: UAI outside of a 2-sided "monogamous" partnership
  ##               with a partner tested negative in past 6 months
  mono.neg <- intersect(uai.mono, all.neg)
  part.id1 <- c(el2[el2$p1 %in% mono.neg, 2], el2[el2$p2 %in% mono.neg, 1])
  part.recently.tested <- since.test[part.id1] <= (180/time.unit)
  mono.neg.recently.tested <- mono.neg[which(part.recently.tested == TRUE)]
  
  uai.risk <- setdiff(uai.any, mono.neg.recently.tested)
  dat$attr$prep.ind.uai.risk[uai.risk] <- at
  
  return(dat)
}

