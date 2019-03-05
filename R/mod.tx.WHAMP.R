
#' @title Treatment Module
#'
#' @description Module function for anti-retroviral treatment initiation and
#'              adherence over time for the WHAMP model.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Persons enter into the simulation with one of four ART "patterns": non-screeners
#' who are diagnosed upon progression to AIDS-associated symptoms and treat to partial 
#' HIV viral suppression, non-screeners who acheive full viral suppression, regular
#' screeners who treat to partial viral suppression, and regular screeners who treat
#' to full viral suppression. These types are stored as individual-level attributes 
#' in \code{tt.traj}). This module initiates ART for treatment-naive persons and 
#' cycles them on and off treatment conditional on overall discontinuation rates for
#' full and partial suppessors and reinitiation rates specific to each group defined 
#' by race/ethnicity, region, and full vs. partial adherence. ART initiation is deterministic,
#' whereas discontinuation and reinitiation are stochastically simulated based on binomial 
#' statistical models.
#'
#' @return
#' This function returns the \code{dat} object with updated \code{tx.status},
#' \code{tx.init.time}, \code{cum.time.on.tx}, \code{cum.time.off.tx} attributes.
#'
#' @keywords module msm
#'
#' @export
#'
tx_msm_whamp <- function(dat, at) {

  ## Variables

  # Attributes
  race..wa <- dat$attr$race..wa
  region <- dat$attr$region
  status <- dat$attr$status
  tx.status <- dat$attr$tx.status
  diag.status <- dat$attr$diag.status
  diag.time <- dat$attr$diag.time
  tt.traj <- dat$attr$tt.traj
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  stage <- dat$attr$stage

  # Parameters
  tx.init.int <- rep(NA, length(race..wa))
  tx.init.int[race..wa == "B" & region == "KC"] <- dat$param$tx.init.int.KC.B
  tx.init.int[race..wa == "H" & region == "KC"] <- dat$param$tx.init.int.KC.H
  tx.init.int[race..wa == "O" & region == "KC"] <- dat$param$tx.init.int.KC.O
  tx.init.int[race..wa == "B" & region == "OW"] <- dat$param$tx.init.int.OW.B
  tx.init.int[race..wa == "H" & region == "OW"] <- dat$param$tx.init.int.OW.H
  tx.init.int[race..wa == "O" & region == "OW"] <- dat$param$tx.init.int.OW.O
  tx.init.int[race..wa == "B" & region == "EW"] <- dat$param$tx.init.int.EW.B
  tx.init.int[race..wa == "H" & region == "EW"] <- dat$param$tx.init.int.EW.H
  tx.init.int[race..wa == "O" & region == "EW"] <- dat$param$tx.init.int.EW.O
  
  tx.halt.full.prob <- dat$param$tx.halt.full
  tx.halt.part.rr <- dat$param$tx.halt.part.rr
  
  tx.reinit.full.KC.B.prob <- dat$param$tx.reinit.full.KC.B
  tx.reinit.full.KC.H.prob <- dat$param$tx.reinit.full.KC.H
  tx.reinit.full.KC.O.prob <- dat$param$tx.reinit.full.KC.O
  tx.reinit.full.OW.B.prob <- dat$param$tx.reinit.full.OW.B
  tx.reinit.full.OW.H.prob <- dat$param$tx.reinit.full.OW.H
  tx.reinit.full.OW.O.prob <- dat$param$tx.reinit.full.OW.O
  tx.reinit.full.EW.B.prob <- dat$param$tx.reinit.full.EW.B
  tx.reinit.full.EW.H.prob <- dat$param$tx.reinit.full.EW.H
  tx.reinit.full.EW.O.prob <- dat$param$tx.reinit.full.EW.O
  
  tx.reinit.part.rr <- dat$param$tx.reinit.part.rr
  

  ## Initiation
  time.since.dx <- at - diag.time
  tx.init <- which(status == 1 & diag.status == 1 &
                        time.since.dx >= tx.init.int &
                        tx.status == 0 & cum.time.on.tx == 0)

  dat$attr$tx.status[tx.init] <- 1
  dat$attr$tx.init.time[tx.init] <- at


  ## Halting
  tx.halt.elig.full <- which(tx.status == 1 & tt.traj %in% c(2,4))
  tx.halt.full <- tx.halt.elig.full[rbinom(length(tx.halt.elig.full), 1,
                                     tx.halt.full.prob) == 1]
  
  tx.halt.elig.part <- which(tx.status == 1 & tt.traj %in% c(1,3))
  tx.halt.part <- tx.halt.elig.part[rbinom(length(tx.halt.elig.part), 1,
                                           (tx.halt.full.prob * tx.halt.part.rr)) == 1]

  tx.halt <- c(tx.halt.full, tx.halt.part)
  dat$attr$tx.status[tx.halt] <- 0


  ## Reinitating
  
  tx_reinit_fun <- function(reg, race.eth, tx.group) {
    tx.reinit.prob <- get(paste0("tx.reinit.full.", reg, ".", race.eth, ".prob"))
    if (tx.group == "full") {
      tx.reinit.elig <- which(race..wa == race.eth & region == reg & tt.traj %in% c(2, 4) &
                                tx.status ==0 & cum.time.on.tx > 0)
      tx.reinit <- tx.reinit.elig[rbinom(length(tx.reinit.elig), 1, tx.reinit.prob) == 1]
    }
    if (tx.group == "part") {
      tx.reinit.elig <- which(race..wa == race.eth & region == reg & tt.traj %in% c(1, 3) &
                                tx.status ==0 & cum.time.on.tx > 0)
      tx.reinit <- tx.reinit.elig[rbinom(length(tx.reinit.elig), 1, tx.reinit.prob * tx.reinit.part.rr) == 1]
    }
    
    return(tx.reinit)
  }
  
  tx.reinit.KC.B.full <- tx_reinit_fun("KC", "B", "full")
  tx.reinit.KC.B.part <- tx_reinit_fun("KC", "B", "part")
  tx.reinit.KC.H.full <- tx_reinit_fun("KC", "H", "full")
  tx.reinit.KC.H.part <- tx_reinit_fun("KC", "H", "part")
  tx.reinit.KC.O.full <- tx_reinit_fun("KC", "O", "full")
  tx.reinit.KC.O.part <- tx_reinit_fun("KC", "O", "part")
  
  tx.reinit.OW.B.full <- tx_reinit_fun("OW", "B", "full")
  tx.reinit.OW.B.part <- tx_reinit_fun("OW", "B", "part")
  tx.reinit.OW.H.full <- tx_reinit_fun("OW", "H", "full")
  tx.reinit.OW.H.part <- tx_reinit_fun("OW", "H", "part")
  tx.reinit.OW.O.full <- tx_reinit_fun("OW", "O", "full")
  tx.reinit.OW.O.part <- tx_reinit_fun("OW", "O", "part")
  
  tx.reinit.EW.B.full <- tx_reinit_fun("EW", "B", "full")
  tx.reinit.EW.B.part <- tx_reinit_fun("EW", "B", "part")
  tx.reinit.EW.H.full <- tx_reinit_fun("EW", "H", "full")
  tx.reinit.EW.H.part <- tx_reinit_fun("EW", "H", "part")
  tx.reinit.EW.O.full <- tx_reinit_fun("EW", "O", "full")
  tx.reinit.EW.O.part <- tx_reinit_fun("EW", "O", "part")
  
  tx.reinit <- c(tx.reinit.KC.B.full, tx.reinit.KC.B.part, tx.reinit.KC.H.full, tx.reinit.KC.H.part, tx.reinit.KC.O.full, tx.reinit.KC.O.part,
                 tx.reinit.OW.B.full, tx.reinit.OW.B.part, tx.reinit.OW.H.full, tx.reinit.OW.H.part, tx.reinit.OW.O.full, tx.reinit.OW.O.part,
                 tx.reinit.EW.B.full, tx.reinit.EW.B.part, tx.reinit.EW.H.full, tx.reinit.EW.H.part, tx.reinit.EW.O.full, tx.reinit.EW.O.part)
  
  dat$attr$tx.status[tx.reinit] <- 1


  ## Other output
  dat$attr$cum.time.on.tx <- dat$attr$cum.time.on.tx +
                             ((dat$attr$tx.status == 1) %in% TRUE)
  dat$attr$cum.time.off.tx <- dat$attr$cum.time.off.tx +
                              ((dat$attr$tx.status == 0) %in% TRUE)

  ## Summary statistics
  dat$epi$tx.init.inc[at] <- length(tx.init)
  dat$epi$tx.halt.inc[at] <- length(tx.halt)
  dat$epi$tx.resm.inc[at] <- length(tx.reinit)

  return(dat)
}


#' @export
#' @rdname tx_msm
tx_het <- function(dat, at) {

  # Variables ---------------------------------------------------------------
  dxStat <- dat$attr$dxStat
  txStat <- dat$attr$txStat
  txStartTime <- dat$attr$txStartTime
  txStops <- dat$attr$txStops
  txTimeOn <- dat$attr$txTimeOn
  txTimeOff <- dat$attr$txTimeOff
  txCD4start <- dat$attr$txCD4start

  cd4Count <- dat$attr$cd4Count
  tx.elig.cd4 <- dat$param$tx.elig.cd4
  tx.coverage <- dat$param$tx.coverage

  txType <- dat$attr$txType
  tx.adhere.full <- dat$param$tx.adhere.full
  tx.adhere.part <- dat$param$tx.adhere.part


  # Start tx for tx naive ---------------------------------------------------

  ## Calculate tx coverage
  allElig <- which((cd4Count < tx.elig.cd4 | !is.na(txStartTime)))
  txCov <- sum(!is.na(txStartTime[allElig]))/length(allElig)
  if (is.nan(txCov)) {
    txCov <- 0
  }

  idsElig <- which(dxStat == 1 & txStat == 0 &
                   is.na(txStartTime) & cd4Count < tx.elig.cd4)
  nElig <- length(idsElig)
  idsTx <- NULL


  ## Treatment coverage
  nStart <- max(0, min(nElig, round((tx.coverage - txCov) * length(allElig))))
  if (nStart > 0) {
    idsTx <- ssample(idsElig, nStart)
  }


  ## Treatment type assignment
  if (length(idsTx) > 0) {
    needtxType <- which(is.na(txType[idsTx]))
    if (length(needtxType) > 0) {
      txType[idsTx[needtxType]] <- rbinom(length(needtxType), 1, tx.adhere.full)
    }
    if (tx.adhere.part == 0) {
      idsTx <- intersect(idsTx, which(txType == 1))
    }
  }

  if (length(idsTx) > 0) {
    txStat[idsTx] <- 1
    txStartTime[idsTx] <- at
    txStops[idsTx] <- 0
    txTimeOn[idsTx] <- 0
    txTimeOff[idsTx] <- 0
    txCD4start[idsTx] <- cd4Count[idsTx]
  }


  # Stop tx -----------------------------------------------------------------
  idsStop <- NULL
  idsEligStop <- which(dat$attr$txStat == 1 & txType == 0)
  nEligStop <- length(idsEligStop)
  if (nEligStop > 0) {
    vecStop <- which(rbinom(nEligStop, 1, (1 - tx.adhere.part)) == 1)
    if (length(vecStop) > 0) {
      idsStop <- idsEligStop[vecStop]
      txStat[idsStop] <- 0
      txStops[idsStop] <- txStops[idsStop] + 1
    }
  }


  # Restart tx --------------------------------------------------------------
  idsRest <- NULL
  idsEligRest <- which(dat$attr$txStat == 0 & txStops > 0)
  nEligRest <- length(idsEligRest)
  if (nEligRest > 0) {
    vecRes <- which(rbinom(nEligRest, 1, tx.adhere.part) == 1)
    if (length(vecRes) > 0) {
      idsRest <- idsEligRest[vecRes]
      txStat[idsRest] <- 1
      dat$attr$vlSlope[idsRest] <- NA
    }
  }


  # Output ------------------------------------------------------------------
  idsOnTx <- which(txStat == 1)
  idsOffTx <- which(txStat == 0 & !is.na(txStartTime))
  txTimeOn[idsOnTx] <- txTimeOn[idsOnTx] + 1
  txTimeOff[idsOffTx] <- txTimeOff[idsOffTx] + 1

  dat$attr$txStat <- txStat
  dat$attr$txStartTime <- txStartTime
  dat$attr$txStops <- txStops
  dat$attr$txTimeOn <- txTimeOn
  dat$attr$txTimeOff <- txTimeOff
  dat$attr$txType <- txType
  dat$attr$txCD4start <- txCD4start

  return(dat)
}

