
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
#' data on age-specific mortality rates applies; and disease-related deaths,
#' which are calculated by multiplying general age-specific mortality rates by a 
#' scalar for persons living with HIV.
#' Which nodes have died is determined stochastically using draws from binomial 
#' distributions for both general and HIV-related deaths. This module also tracks 
#' nodes that age out of the network.
#'
#' @return
#' This function returns the updated \code{dat} object accounting for deaths.
#' The deaths are deactivated from the main and persistent networks, as those are in
#' \code{networkDynamic} class objects; dead nodes are not deleted from the
#' instant network until the \code{\link{simnet_msm}} module for bookkeeping
#' purposes.
#'
#' @keywords module msm
#' @export
#'
deaths_msm_whamp <- function(dat, at) {
  
  ## General deaths among HIV-neg (includes deaths from aging out of the population b/c asmr = 1 at age 60)
  age <- floor(dat$attr$age)
  race..wa <- dat$attr$race..wa
  status <- dat$attr$status
  
  alive.H..wa.neg <- which(race..wa == "H" & status == 0)
  age.H..wa.neg <- age[alive.H..wa.neg]
  death.H.prob..wa.neg <-  dat$param$asmr.H..wa[age.H..wa.neg]
  deaths.H..wa.neg <- alive.H..wa.neg[rbinom(length(death.H.prob..wa.neg), 1, death.H.prob..wa.neg) == 1]
  
  alive.B..wa.neg <- which(race..wa == "B" & status == 0)
  age.B..wa.neg <- age[alive.B..wa.neg]
  death.B.prob..wa.neg <-  dat$param$asmr.B..wa[age.B..wa.neg]
  deaths.B..wa.neg <- alive.B..wa.neg[rbinom(length(death.B.prob..wa.neg), 1, death.B.prob..wa.neg) == 1]
  
  alive.O..wa.neg <- which(race..wa == "O" & status == 0)
  age.O..wa.neg <- age[alive.O..wa.neg]
  death.O.prob..wa.neg <-  dat$param$asmr.O..wa[age.O..wa.neg]
  deaths.O..wa.neg <- alive.O..wa.neg[rbinom(length(death.O.prob..wa.neg), 1, death.O.prob..wa.neg) == 1]
  
  dth.neg..wa <- c(deaths.H..wa.neg, deaths.B..wa.neg, deaths.O..wa.neg)
  
  
  ## Deaths among HIV-pos
  asmr.H..wa.pos <-  c(dat$param$asmr.H..wa[1:17], 
                       dat$param$asmr.H..wa[18:44]*dat$param$asmr.rr.pos[1], 
                       dat$param$asmr.H..wa[45:54]*dat$param$asmr.rr.pos[2],
                       dat$param$asmr.H..wa[55:59]*dat$param$asmr.rr.pos[3],
                       dat$param$asmr.H..wa[60])
  
  asmr.B..wa.pos <-  c(dat$param$asmr.B..wa[1:17], 
                       dat$param$asmr.B..wa[18:44]*dat$param$asmr.rr.pos[1], 
                       dat$param$asmr.B..wa[45:54]*dat$param$asmr.rr.pos[2],
                       dat$param$asmr.B..wa[55:59]*dat$param$asmr.rr.pos[3],
                       dat$param$asmr.B..wa[60])
  
  asmr.O..wa.pos <-  c(dat$param$asmr.O..wa[1:17], 
                       dat$param$asmr.O..wa[18:44]*dat$param$asmr.rr.pos[1], 
                       dat$param$asmr.O..wa[45:54]*dat$param$asmr.rr.pos[2],
                       dat$param$asmr.O..wa[55:59]*dat$param$asmr.rr.pos[3],
                       dat$param$asmr.O..wa[60])
  
  
  alive.H..wa.pos <- which(race..wa == "H" & status == 1)
  age.H..wa.pos <- age[alive.H..wa.pos]
  death.H.prob..wa.pos <-  asmr.H..wa.pos[age.H..wa.pos]
  deaths.H..wa.pos <- alive.H..wa.pos[rbinom(length(death.H.prob..wa.pos), 1, death.H.prob..wa.pos) == 1]
  
  alive.B..wa.pos <- which(race..wa == "B" & status == 1)
  age.B..wa.pos <- age[alive.B..wa.pos]
  death.B.prob..wa.pos <-  asmr.B..wa.pos[age.B..wa.pos]
  deaths.B..wa.pos <- alive.B..wa.pos[rbinom(length(death.B.prob..wa.pos), 1, death.B.prob..wa.pos) == 1]
  
  alive.O..wa.pos <- which(race..wa == "O" & status == 1)
  age.O..wa.pos <- age[alive.O..wa.pos]
  death.O.prob..wa.pos <-  asmr.O..wa.pos[age.O..wa.pos]
  deaths.O..wa.pos <- alive.O..wa.pos[rbinom(length(death.O.prob..wa.pos), 1, death.O.prob..wa.pos) == 1]
  
  dth.pos..wa <- c(deaths.H..wa.pos, deaths.B..wa.pos, deaths.O..wa.pos)
  
  ### Calculate the excepted number of deaths among HIV-pos if they experienced the background general pop ASMRs 
  ### (for use in determining the number of births to balance out background mortality
  death.H.prob..wa.pos.bkg <-  dat$param$asmr.H..wa[age.H..wa.pos]
  deaths.H..wa.pos.bkg <- sum(rbinom(length(death.H.prob..wa.pos.bkg), 1, death.H.prob..wa.pos.bkg) == 1)
  
  death.B.prob..wa.pos.bkg <-  dat$param$asmr.B..wa[age.B..wa.pos]
  deaths.B..wa.pos.bkg <- sum(rbinom(length(death.B.prob..wa.pos.bkg), 1, death.B.prob..wa.pos.bkg) == 1)
  
  death.O.prob..wa.pos.bkg <-  dat$param$asmr.O..wa[age.O..wa.pos]
  deaths.O..wa.pos.bkg <- sum(rbinom(length(death.O.prob..wa.pos.bkg), 1, death.O.prob..wa.pos.bkg) == 1)
  
  dth.pos.bkg <- sum(deaths.H..wa.pos.bkg, deaths.B..wa.pos.bkg, deaths.O..wa.pos.bkg)
  
  ## Combine
  dth.all..wa <- NULL
  dth.all..wa <- unique(c(dth.neg..wa, dth.pos..wa))
  
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
  dth.age <- which(age >= dat$param$exit.age)
  dth.age.neg <- which(age >= dat$param$exit.age & status == 0) # define these separately so can subract them below
  dth.age.pos <- which(age >= dat$param$exit.age & status == 1)
  
  ## Summary Output
  dat$epi$dth.neg..wa[at] <- max(0, length(dth.neg..wa) - length(dth.age.neg)) #subtract age deaths to track them separately, b/c dth.neg and dth.pos include aging out
  dat$epi$dth.pos..wa[at] <- max(0, length(dth.pos..wa) - length(dth.age.pos))
  dat$epi$dth.age[at] <- max(0, length(dth.age))
  dat$epi$dth.pos.bkg[at] <- max(0, dth.pos.bkg - length(dth.age.pos))
  
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
