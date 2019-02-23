
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

  if (at < dat$param$prep.start) {
    return(dat)
  }

  ## Variables

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
  # prepLastStiScreen <- dat$attr$prepLastStiScreen

  # Parameters
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


  ## Discontinuation ------------------------------------------------------------------

  # No longer indicated for PrEP
  if (dat$param$prep.risk.reassess == TRUE) {  # If TRUE, reassess at every testing visit
    idsRiskAssess <- which(active == 1 & prepStat == 1 & lnt == at)

    idsEligStop <- intersect(which(ind1 < at & ind2 < twind),
                             idsRiskAssess)
  
    prepElig[idsEligStop] <- 0

  } 
  if (dat$param$prep.risk.reassess == FALSE) { #  If FALSE, reassess at every year
    idsRiskAssess <- which(active == 1 & prepStat == 1 & lnt == at & (at - prepLastRisk) >= 52)
    prepLastRisk[idsRiskAssess] <- at
    
    idsEligStop <- intersect(which(ind1 < at & ind2 < twind),
                             idsRiskAssess)
    
    prepElig[idsEligStop] <- 0
    
  } 
  
  # Spontaneous discontinuation
  discont.elig <- which(active == 1 & prepStat == 1 & prepDiscont == 1)
  idsDiscont <- discont.elig[rbinom(length(discont.elig), 1,
                                     prep.discont.prob) == 1]

  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & diag.status == 1)

  # Death
  idsStpDth <- which(active == 0 & prepStat == 1)

  # Reset PrEP status
  idsStp <- c(idsStpDx, idsStpDth, idsEligStop, idsDiscont)
  prepStat[idsStp] <- 0
  prepLastRisk[idsStp] <- NA
  prepStartTime[idsStp] <- NA

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
  
  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepStartTime <- prepStartTime
  dat$attr$prepClass <- prepClass
  dat$attr$prepLastRisk <- prepLastRisk
  dat$attr$prepLastStiScreen <- prepLastStiScreen
  dat$attr$prepDiscont <- prepDiscont

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

  return(dat)
}
