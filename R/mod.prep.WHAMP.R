
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
    dat <- riskhist_msm_whamp(dat, at)
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
  tt.traj <- dat$attr$tt.traj
  lnt <- dat$attr$last.neg.test
  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  prepLastRisk <- dat$attr$prepLastRisk
  prepStartTime <- dat$attr$prepStartTime
  prepDiscont <- dat$attr$prepDisc
  spontDisc <- dat$attr$spontDisc
  # prepLastStiScreen <- dat$attr$prepLastStiScreen

  # Parameters
  prep.risk.reassess.method <- dat$param$prep.risk.reassess.method
  prep.start.step <- dat$param$prep.start
  prep.init.rate <- dat$param$prep.init.rate
  prep.uptake.scen <- dat$param$prep.uptake.scen
  if(prep.uptake.scen == "lower"){
    prep.coverage.init <- dat$param$prep.coverage.init.low
  } else {
    prep.coverage.init <- dat$param$prep.coverage.init
  }
  prep.cov.max <- dat$param$prep.cov.max
  prep.discontinue <- dat$param$prep.discont
  prep.discont.prob <- dat$param$prep.discont.prob
  if(dat$param$prep.adh == "high"){
    prep.class.prob <- dat$param$prep.class.prob.high
  } else {
    prep.class.prob <- dat$param$prep.class.prob.low
  }

  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(active == 1 & (diag.status == 0 | is.na(diag.status)) & 
                          prepStat == 0 & lnt == at & tt.traj %in% c(3,4))

  # Core eligiblity
  ind1 <- dat$attr$prep.ind.discord.ongoing
  ind2 <- dat$attr$prep.ind.uai.risk

  twind <- at - dat$param$prep.risk.int
  
  prepElig[which((ind1 == at | ind2 >= twind) & (diag.status == 0 | is.na(diag.status)))] <- 1
  
  idsEligStart <- intersect(which(prepElig == 1),
                            idsEligStart)
  
  # No longer indicated for PrEP
  idsNoIndic <- which(((ind1 < at) | is.na(ind1)) & 
                        ((ind2 < twind) | is.na(ind2)))
  
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
  discont.elig <- which(active == 1 & prepStat == 1 & prepDiscont == 1 & (diag.status == 0 | is.na(diag.status)))
  idsDiscont <- discont.elig[rbinom(length(discont.elig), 1,
                                    prep.discont.prob) == 1]
  
  spontDisc[idsDiscont] <- 1
  
  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & diag.status == 1)
  prepElig[idsStpDx] <- 0

  # Death
  idsStpDth <- which(active == 0 & prepStat == 1)

  # All discontinuation
  idsStp <- c(idsStpDx, idsStpDth, idsEligStop, idsDiscont)
  
  # Time to discontinuation among those who stopped
    # time.to.disc <- rep(NA, length(idsStp))
    # time.to.disc <- at - prepStartTime[idsStp]
 
  
  # Reset PrEP status
  prepStat[idsStp] <- 0
  prepLastRisk[idsStp] <- NA
  prepStartTime[idsStp] <- NA
  # prepLastStiScreen[idsStp] <- NA
  

  ## Initiation ----------------------------------------------------------------
  
  # Calculate the current coverage among eligible men
  prep.cov.curr <- sum(prepStat == 1, na.rm = TRUE) / sum(prepElig ==1, na.rm = TRUE)
  
  # Determine coverage quota for the current time step
  if (at == prep.start.step){
    prep.coverage <- prep.coverage.init
  }
  
  if (at > prep.start.step & is.finite(prep.start.step)){
    # Stable PrEP use scenarios
    if (prep.uptake.scen %in% c("stable", "lower")){
      prep.coverage <- prep.coverage.init
    } 
    # 50% uptake by 2020 (over 3*52 time steps from PrEP start in 2017)
    if (prep.uptake.scen == "fifty"){
      prep.cov.max <- 0.5
      prep.scaleup.rate <- (prep.cov.max - prep.coverage.init) / (3*52)
      prep.coverage <- min((prep.coverage.init + prep.scaleup.rate*(at - prep.start.step)), prep.cov.max)
    }
    # Increase to maximum coverage in each region, after reaching 50% overall by 2020
    if (prep.uptake.scen == "max"){
      prep.cov.2020 <- 0.5
      prep.scaleup.rate <- (prep.cov.2020 - prep.coverage.init) / (3*52) # set scaleup rate based on rate needed to reach 50% by 2020
      prep.coverage <- min((prep.coverage.init + prep.scaleup.rate*(at - prep.start.step)), prep.cov.max)
    }
  }
  
  # Assign PrEP use
  
    ## If at = prep.start.step, assing PrEP to eligible individuals such that coverage = prep.coverage.init.
    ## For this step, do not require men to test in the current time step to be eligible, but only regular testers can start PrEP
    if(at == prep.start.step){
      idsSusc <- which(active == 1 & (diag.status == 0 | is.na(diag.status)) & prepStat == 0 & tt.traj %in% c(3,4))
      idsElig.init <- intersect(which(prepElig == 1), idsSusc)
      
      nStart <- pmax(0, pmin(length(idsElig.init), round(prep.coverage * sum(prepElig == 1, na.rm = TRUE))))
      
      idsStart <- NULL
      if (nStart > 0){
        if (prep.init.rate >= 1) {
          idsStart <- ssample(idsElig.init, nStart)
        } else {
          idsStart <- idsElig.init[rbinom(nStart, 1, prep.init.rate) == 1]
        }
      }
    }
  
  
    ## If at > prep.start.step, sample ids to start PrEP if current coverage < coverage threshold
    if (at > prep.start.step & is.finite(prep.start.step)){
      nStart <- pmax(0, pmin(length(idsEligStart), round((prep.coverage - prep.cov.curr)*sum(prepElig == 1, na.rm = TRUE))))
      
      idsStart <- NULL
      if (nStart > 0) {
        if (prep.init.rate >= 1) {
          idsStart <- ssample(idsEligStart, nStart)
        } else {
          idsStart <- idsEligStart[rbinom(nStart, 1, prep.init.rate) == 1]
        }
      }
    }
  
 
  # Assign PrEP-related attributes
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
    prepDiscont[idsStart[needgp]] <- rbinom(length(needgp), 1, prep.discontinue) 
    
  }
  
  ## Calculate time on PrEP among current users -----------------------
  # time.since.prep.start <- rep(NA, sum(prepStat == 1, na.rm = TRUE))
  # time.since.prep.start <- (at - prepStartTime[prepStat == 1]) 

  
  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepStartTime <- prepStartTime
  dat$attr$prepClass <- prepClass
  dat$attr$prepLastRisk <- prepLastRisk
  # dat$attr$prepLastStiScreen <- prepLastStiScreen
  dat$attr$prepDisc <- prepDiscont
  dat$attr$spontDisc <- spontDisc

  # Summary Statistics
  prep.cov.elig <- sum(prepStat == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  prep.cov.all <- sum(prepStat == 1, na.rm = TRUE)/sum((diag.status == 0 | is.na(diag.status)), na.rm = TRUE)

  dat$epi$prep.cov.elig[at] <- prep.cov.elig
  dat$epi$prep.cov.all[at] <- prep.cov.all
  
  dat$epi$prepStart[at] <- length(idsStart)
  # dat$epi$time.to.disc <- mean(time.to.disc)
  # dat$epi$time.since.prep.start <- mean(time.since.prep.start)

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
  prep.ind.plt.int <- dat$param$prep.ind.plt.int
  
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
  
  
  ## Condition 1: ongoing positive partner who has disclosed
  
  discord.ong <- el2[el2$st1 == 1 & el2$ptype %in% c(1:2), ]
  
  # Disclosure
  discl.list <- dat$temp$discl.list
  disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
  delt.cdl <- uid[discord.ong[, 1]] * 1e7 + uid[discord.ong[, 2]]
  discl <- (delt.cdl %in% disclose.cdl)
  ai.sd <- discord.ong$p2[discl == TRUE]
  
  dat$attr$prep.ind.discord.ongoing[ai.sd] <- at
  
  ## Condition 2: UAI outside of a 2-sided "monogamous" partnership
  ##               with a partner tested negative in past 6 months
  
  # Any UAI
  uai.any <- unique(c(el2$p1[el2$uai > 0],
                      el2$p2[el2$uai > 0]))
  
  # Monogamous partnerships (2-sided)
  tot.deg <- main.deg + casl.deg + inst.deg
  tot.deg.p1 <- tot.deg[el2$p1]
  tot.deg.p2 <- tot.deg[el2$p2]
  monog.2sided <- unique(c(el2$p1[tot.deg.p1 == 1 & tot.deg.p2 == 1],
                           el2$p2[tot.deg.p1 == 1 & tot.deg.p2 == 1]))
  
  # "Negative" partnerships
  tneg <- unique(c(el2$p1[el2$st1 == 0], el2$p2[el2$st1 == 0]))
  fneg <- unique(c(el2$p1[which(dx[el2$p1] == 0)], el2$p2[which(dx[el2$p1] == 0)]))
  all.neg <- c(tneg, fneg)
  

  uai.mono <- intersect(monog.2sided, uai.any)
  uai.mono.neg <- intersect(uai.mono, all.neg)
  part.id1 <- c(el2[el2$p1 %in% uai.mono.neg, 2], el2[el2$p2 %in% uai.mono.neg, 1])
  part.recently.tested <- since.test[part.id1] <= prep.ind.plt.int
  mono.neg.recently.tested <- uai.mono.neg[which(part.recently.tested == TRUE)]
  
  uai.risk <- setdiff(uai.any, mono.neg.recently.tested)
  dat$attr$prep.ind.uai.risk[uai.risk] <- at
  
  return(dat)
}

