
#' @title Death Module
#'
#' @description Module function for simulting both general and disease-related
#'              deaths, as well as aging out of the network, among population 
#'              members for the WHAMP model.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Deaths are divided into two categories: general deaths, for which demographic
#' data on age-specific mortality rates applies; and disease-related diseases,
#' for which the rate of death is a function of progression to end-stage AIDS.
#' Which nodes have died is determined stochastically for general deaths using
#' draws from a binomial distribution, and deterministically for disease-related
#' deaths after nodes have reach a maximum viral load value set in the
#' \code{vl.fatal} parameter. This module also tracks nodes that age out of the network.
#'
#' @return
#' This function returns the updated \code{dat} object accounting for deaths.
#' The deaths are deactivated from the main and casual networks, as those are in
#' \code{networkDynamic} class objects; dead nodes are not deleted from the
#' instant network until the \code{\link{simnet_msm}} module for bookkeeping
#' purposes.
#'
#' @keywords module msm
#' @export
#'
deaths_msm_whamp <- function(dat, at) {

  ## General deaths (includes deaths from aging out of the population b/c asmr = 1 at age 60)
  age <- floor(dat$attr$age)
  race <- dat$attr$race   #-- Delete this old code eventually
  race..wa <- dat$attr$race..wa
  
  #-- Delete this old code eventually
    alive.B <- which(race == "B")
    age.B <- age[alive.B]
    death.B.prob <- dat$param$asmr.B[age.B]
    deaths.B <- alive.B[rbinom(length(death.B.prob), 1, death.B.prob) == 1]
  
    alive.W <- which(race == "W")
    age.W <- age[alive.W]
    death.W.prob <- dat$param$asmr.W[age.W]
    deaths.W <- alive.W[rbinom(length(death.W.prob), 1, death.W.prob) == 1]
  
  alive.H..wa <- which(race..wa == "H")
  age.H..wa <- age[alive.H..wa]
  death.H.prob..wa <-  dat$param$asmr.H..wa[age.H..wa]
  deaths.H..wa <- alive.H..wa[rbinom(length(death.H.prob..wa), 1, death.H.prob..wa) == 1]
  
  alive.B..wa <- which(race..wa == "B")
  age.B..wa <- age[alive.B..wa]
  death.B.prob..wa <-  dat$param$asmr.B..wa[age.B..wa]
  deaths.B..wa <- alive.B..wa[rbinom(length(death.B.prob..wa), 1, death.B.prob..wa) == 1]
  
  alive.O..wa <- which(race..wa =="O")
  age.O..wa <- age[alive.O..wa]
  death.O.prob..wa <-  dat$param$asmr.O..wa[age.O..wa]
  deaths.O..wa <- alive.O..wa[rbinom(length(death.O.prob..wa), 1, death.O.prob..wa) == 1]

  dth.gen <- c(deaths.B, deaths.W) #-- Delete this when finish debugging
  dth.gen..wa <- c(deaths.H..wa, deaths.B..wa, deaths.O..wa)


  ## Disease deaths
  dth.dis <- which(dat$attr$stage == 4 &
                   dat$attr$vl >= dat$param$vl.fatal)

  dth.all <- NULL  #-- Delete this when finish debugging
  dth.all <- unique(c(dth.gen, dth.dis)) #-- Delete this when finish debugging
  dth.all..wa <- NULL
  dth.all..wa <- unique(c(dth.gen..wa, dth.dis))

  if (length(dth.all..wa) > 0) {
    dat$attr$active[dth.all..wa] <- 0
    for (i in 1:3) {
      dat$el[[i]] <- tergmLite::delete_vertices(dat$el[[i]], dth.all..wa)
    }
    dat$attr <- deleteAttr(dat$attr, dth.all..wa)
    if (unique(sapply(dat$attr, length)) != attributes(dat$el[[1]])$n) {
      stop("mismatch between el and attr length in death mod")
    }
  }

  ## Aging out of the network
  dth.age<-which(age >= dat$param$exit.age)
  
  ## Summary Output
  dat$epi$dth.gen..wa[at] <- max(0, length(dth.gen..wa)-length(dth.age)) #subtract age deaths to track them separately, b/c dth.gen includes # aging out
  dat$epi$dth.dis[at] <- max(0, length(dth.dis))
  dat$epi$dth.age[at] <- max(0, length(dth.age))

  return(dat)
}


#' @export
#' @rdname deaths_msm
deaths_het <- function(dat, at) {

  ### 1. Susceptible Deaths ###

  ## Variables
  male <- dat$attr$male
  age <- dat$attr$age
  cd4Count <- dat$attr$cd4Count

  di.cd4.aids <- dat$param$di.cd4.aids
  ds.exit.age <- dat$param$ds.exit.age

  ## Eligible are: active uninf, pre-death infected, unhealthy old
  idsEligSus <- which((is.na(cd4Count) |
                       cd4Count > di.cd4.aids |
                       (cd4Count <= di.cd4.aids & age > ds.exit.age)))
  nEligSus <- length(idsEligSus)

  # Set age-sex specific rates
  ds.rates <- dat$param$ds.rates
  if (nEligSus > 0) {
    rates <- ds.rates$mrate[100*male[idsEligSus] + age[idsEligSus]]
  }


  ## Process
  nDeathsSus <- 0; idsDeathsSus <- NULL
  if (nEligSus > 0) {
    vecDeathsSus <- which(rbinom(nEligSus, 1, rates) == 1)
    nDeathsSus <- length(vecDeathsSus)
  }


  ## Update Attributes
  if (nDeathsSus > 0) {
    idsDeathsSus <- idsEligSus[vecDeathsSus]
    dat$attr$active[idsDeathsSus] <- 0
  }


  ### 2. Infected Deaths ###

  ## Variables
  active <- dat$attr$active
  di.cd4.rate <- dat$param$di.cd4.rate

  ## Process
  nDeathsInf <- 0; idsDeathsInf <- NULL

  cd4Count <- dat$attr$cd4Count
  stopifnot(length(active) == length(cd4Count))

  idsEligInf <- which(active == 1 & cd4Count <= di.cd4.aids)
  nEligInf <- length(idsEligInf)

  if (nEligInf > 0) {
    vecDeathsInf <- which(rbinom(nEligInf, 1, di.cd4.rate) == 1)
    if (length(vecDeathsInf) > 0) {
      idsDeathsInf <- idsEligInf[vecDeathsInf]
      nDeathsInf <- length(idsDeathsInf)
    }
  }

  idsDeathsDet <- which(cd4Count <= 0)
  if (length(idsDeathsDet) > 0) {
    idsDeathsInf <- c(idsDeathsInf, idsDeathsDet)
    nDeathsInf <- nDeathsInf + length(idsDeathsDet)
  }


  ### 3. Update Attributes ###
  if (nDeathsInf > 0) {
    dat$attr$active[idsDeathsInf] <- 0
  }

  ## 4. Update Population Structure ##
  inactive <- which(dat$attr$active == 0)
  dat$el[[1]] <- tergmLite::delete_vertices(dat$el[[1]], inactive)
  dat$attr <- deleteAttr(dat$attr, inactive)

  if (unique(sapply(dat$attr, length)) != attributes(dat$el[[1]])$n) {
    stop("mismatch between el and attr length in death mod")
  }

  ### 5. Summary Statistics ###
  dat$epi$ds.flow[at] <- nDeathsSus
  dat$epi$di.flow[at] <- nDeathsInf

  return(dat)
}
