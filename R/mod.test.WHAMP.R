
#' @title HIV Testing Module
#'
#' @description Module function for HIV diagnostic testing of infected persons for the WHAMP model.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This testing module supports two testing parameterizations, input via the
#' \code{testing.pattern} parameter: memoryless for stochastic and
#' geometrically-distributed waiting times to test (constant hazard); and interval
#' for deterministic tested after defined waiting time intervals.
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
  race <- dat$attr$race #-- Delete when finish debugging
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

  tsincelntst <- at - dat$attr$last.neg.test
  tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]

  # Calculate intertest interval as a function of age
  centered.age <- (age - mean.age.iti)
  avg.test.int <- iti.coefs[1] + centered.age * iti.coefs[2] + centered.age^2 * iti.coefs[3]
  
  ## Process

  if (testing.pattern == "memoryless") {
    stop("Intertest interval parameter calculated assuming interval method. Revise parameter estimation procedure for memoryless process.")
  }

  if (testing.pattern == "interval") {
    tst <- which(tt.traj != 1 &
                   (diag.status == 0 | is.na(diag.status)) &
                   tsincelntst >= (avg.test.int) &
                   prepStat == 0)
    tst.nprep <- tst
  }

  # PrEP testing
  tst.prep <- which((diag.status == 0 | is.na(diag.status)) &
                    prepStat == 1 &
                    tsincelntst >= prep.tst.int)

  tst.all <- c(tst.nprep, tst.prep)

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

