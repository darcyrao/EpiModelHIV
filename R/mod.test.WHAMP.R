
#' @title HIV Testing Module
#'
#' @description Module function for HIV diagnostic testing of infected persons for the WHAMP model.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This testing module supports testing as an interval process and diagnostic testing 
#' upon onset of AIDS-associated symptoms
#'
#' @return
#' This function returns the \code{dat} object with updated \code{last.neg.test},
#' \code{diag.status} and \code{diag.time} attributes.
#'
#' @keywords module msm
#'
#' @export
#'
test_msm_whamp <- function(dat, at) {

  ## Variables

  # Attributes
  age <- dat$attr$age
  diag.status <- dat$attr$diag.status
  diag.time <- dat$attr$diag.time
  race..wa <- dat$attr$race..wa
  tt.traj <- dat$attr$tt.traj
  status <- dat$attr$status
  inf.time <- dat$attr$inf.time

  prepStat <- dat$attr$prepStat
  prep.tst.int <- dat$param$prep.tst.int

  # Parameters
  testing.pattern <- dat$param$testing.pattern
  mean.age.iti <- dat$param$mean.age.iti
  iti.coefs <- dat$param$iti.coefs
  twind.int <- dat$param$test.window.int
  sympt.int <- dat$param$sympt.onset.int

  tsincelntst <- at - dat$attr$last.neg.test
  tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]

  # Calculate intertest interval as a function of age
  centered.age <- (age - mean.age.iti)
  test.int <- (iti.coefs[1] + centered.age * iti.coefs[2] + centered.age^2 * iti.coefs[3]) / dat$param$time.unit
  
  ## Process

  if (testing.pattern == "memoryless") {
    stop("Intertest interval parameter calculated assuming interval method. Revise parameter estimation procedure for memoryless process.")
  }

  # Regular screeners
  if (testing.pattern == "interval") {
    tst <- which(tt.traj %in% c(3, 4) &
                   (diag.status == 0 | is.na(diag.status)) &
                   tsincelntst >= test.int &
                   prepStat == 0)
    tst.nprep <- tst
  }
  
  # Diagnostic testing at symptom onset
  tst.sympt <- which(status == 1 &
                     (diag.status == 0 | is.na(diag.status)) &
                     (at - inf.time) >= sympt.int)
  
  # PrEP testing
  tst.prep <- which((diag.status == 0 | is.na(diag.status)) &
                    prepStat == 1 &
                    tsincelntst >= prep.tst.int)

  tst.all <- c(tst.nprep, tst.sympt, tst.prep)

  tst.pos <- tst.all[status[tst.all] == 1 & inf.time[tst.all] <= at - twind.int]
  tst.neg <- setdiff(tst.all, tst.pos)

  # Attributes
  dat$attr$last.neg.test[tst.neg] <- at
  dat$attr$diag.status[tst.pos] <- 1
  dat$attr$diag.time[tst.pos] <- at

  return(dat)
}


#' @export
#' @rdname test_msm
dx_het <- function(dat, at) {

  # Variables
  status <- dat$attr$status
  txCD4min <- dat$attr$txCD4min
  cd4Count <- dat$attr$cd4Count
  dxStat <- dat$attr$dxStat

  # Process
  tested <- which(status == 1 & dxStat == 0 & cd4Count <= txCD4min)


  # Results
  if (length(tested) > 0) {
    dat$attr$dxStat[tested] <- 1
    dat$attr$txStat[tested] <- 0
    dat$attr$dxTime[tested] <- at
  }

  return(dat)
}

