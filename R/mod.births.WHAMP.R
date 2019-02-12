
#' @title Births Module
#'
#' @description Module function for births or entries into the sexually active
#'              population for the WHAMP model
#'
#' @inheritParams aging_msm
#'
#' @details
#' New population members are added based on expected numbers of entries,
#' stochastically determined with draws from Poisson distributions. The 
#' expected number of entries is set to balance exits due to aging out and  
#' background mortality, defined as the number of deaths among HIV-negative men plus
#' the expected number of deaths among HIV+ if they were subject to general pop
#' age-specific mortality (i.e. subtracting excess mortality due to disease).
#' For each new entry, a set of attributes is added for that node,
#' and the nodes are added onto the network objects. Only attributes that are
#' a part of the network model formulae are updated as vertex attributes on the
#' network objects.
#'
#' @return
#' This function updates the \code{attr} list with new attributes for each new
#' population member, and the \code{nw} objects with new vertices.
#'
#' @keywords module msm
#' @export
#'
births_msm_whamp <- function(dat, at){

  ## Determine the number of births needed to balance general mortality and the number aging out
  
  nBirths.gen <- dat$epi$dth.neg..wa[at] + dat$epi$dth.pos.bkg[at] # Births to replace those who died from general mortality
  nBirths.age <- dat$epi$dth.age[at] # Births to replace those who aged out
  
  ## Determine the number of additional births to add to achieve the specified population growth
  nBirths.growth <- sum(dat$epi$num.H..wa[at - 1], dat$epi$num.B..wa[at - 1], dat$epi$num.O..wa[at - 1])*(dat$param$growth.rate - 1)

  nBirths <- nBirths.gen + nBirths.age + nBirths.growth
  
  ## Update Attr
  
  if (nBirths > 0) {
    dat <- setBirthAttr_msm_whamp(dat, at, nBirths)
  }


  # Update Networks
  if (nBirths > 0) {
    for (i in 1:3) {
      dat$el[[i]] <- tergmLite::add_vertices(dat$el[[i]], nBirths)
    }
  }


  ## Output
  dat$epi$nBirths[at] <- nBirths

  return(dat)
}


setBirthAttr_msm_whamp <- function(dat, at, nBirths) {

  # Set all attributes NA by default
  dat$attr <- lapply(dat$attr, {
    function(x)
      c(x, rep(NA, nBirths))
  })
  newIds <- which(is.na(dat$attr$active))

  ## Demographic
  dat$attr$active[newIds] <- rep(1, nBirths)
  dat$attr$uid[newIds] <- dat$temp$max.uid + (1:nBirths)
  dat$temp$max.uid <- dat$temp$max.uid + nBirths

  dat$attr$arrival.time[newIds] <- rep(at, nBirths)
  
  # Distribute new births by race/ethnicity and region according to proportion of the population in each group
    
    #-- Delete when finish debugging
    race <- sample(apportion_lr(nBirths, c("B", "W"), c(0.5, 0.5))) #-- Assign old race attribute until finish revising code
    newB <- which(race == "B") #-- Delete when finish de-bugging
    newW <- which(race == "W") #-- Delete when finish de-bugging
    nBirths.B <- length(newB) #-- Delete when finish de-bugging
    nBirths.W <- length(newW) #-- Delete when finish de-bugging
    dat$attr$race[newIds] <- race #-- Delete when finish de-bugging
  
  ## Vector with the proportion of Hispanic, black, and other race/ethnicity men in King County, other western WA, and eastern WA
  prop.race.region <- sumto1(c(0.0549, 0.0421, 0.4739, 0.0309, 0.0166, 0.2807, 0.0222, 0.0021, 0.0767))
  
  race.region <- sample(apportion_lr(nBirths, c("H.KC", "B.KC", "O.KC", "H.OW", "B.OW", "O.OW", "H.EW", "B.EW", "O.EW"), 
                                    prop.race.region))
  
  race..wa <- rep(NA, nBirths)
  race..wa[race.region %in% c("H.KC", "H.OW", "H.EW")] <- "H"
  race..wa[race.region %in% c("B.KC", "B.OW", "B.EW")] <- "B"
  race..wa[race.region %in% c("O.KC", "O.OW", "O.EW")] <- "O"
  newH..wa <- which(race..wa == "H")
  newB..wa <- which(race..wa == "B")
  newO..wa <- which(race..wa == "O")
  nBirths.H..wa <- length(newH..wa)
  nBirths.B..wa <- length(newB..wa)
  nBirths.O..wa <- length(newO..wa)
  
  region <- rep(NA, nBirths)
  region[race.region %in% c("H.KC", "B.KC", "O.KC")] <- "KC"
  region[race.region %in% c("H.OW", "B.OW", "O.OW")] <- "OW"
  region[race.region %in% c("H.EW", "B.EW", "O.EW")] <- "EW"
  
  newKC <- which(region == "KC") #-- May not need this section of code. if don't use it in other parts of this file, delete
  newOW <- which(region == "OW")
  newEW <- which(region == "EW")
  nBirths.KC <- length(newKC)
  nBirths.OW <- length(newOW)
  nBirths.EW <- length(newEW)
  
  dat$attr$race..wa[newIds] <- race..wa
  dat$attr$region[newIds] <- region
  

  # New entries will enter at age 18
  dat$attr$age[newIds] <- rep(dat$param$birth.age, nBirths)
  dat$attr$sqrt.age[newIds] <- sqrt(dat$attr$age[newIds])

  # Disease status and related
  dat$attr$status[newIds] <- rep(0, nBirths)

  dat$attr$tt.traj[newIds[newB]] <- sample(c(1, 2, 3, 4),
                                           nBirths.B, replace = TRUE,
                                           prob = dat$param$tt.traj.B.prob)
  dat$attr$tt.traj[newIds[newW]] <- sample(c(1, 2, 3, 4),
                                           nBirths.W, replace = TRUE,
                                           prob = dat$param$tt.traj.W.prob)

  # Circumcision
  dat$attr$circ[newIds[newB]] <- rbinom(nBirths.B, 1, dat$param$circ.B.prob)
  dat$attr$circ[newIds[newW]] <- rbinom(nBirths.W, 1, dat$param$circ.W.prob)

  # Role
  dat$attr$role.class[newIds[newB]] <- sample(c("I", "R", "V"),
                                              nBirths.B, replace = TRUE,
                                              prob = dat$param$role.B.prob)
  dat$attr$role.class[newIds[newW]] <- sample(c("I", "R", "V"),
                                              nBirths.W, replace = TRUE,
                                              prob = dat$param$role.W.prob)

  ins.quot <- rep(NA, nBirths)
  ins.quot[dat$attr$role.class[newIds] == "I"]  <- 1
  ins.quot[dat$attr$role.class[newIds] == "R"]  <- 0
  ins.quot[dat$attr$role.class[newIds] == "V"]  <-
                                  runif(sum(dat$attr$role.class[newIds] == "V"))
  dat$attr$ins.quot[newIds] <- ins.quot

  # CCR5
  ccr5.B.prob <- dat$param$ccr5.B.prob
  ccr5.H.prob <- dat$param$ccr5.H.prob
  ccr5.O.prob <- dat$param$ccr5.O.prob
  dat$attr$ccr5[newIds[newB..wa]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.B..wa, replace = TRUE,
                                        prob = c(1 - sum(ccr5.B.prob),
                                                 ccr5.B.prob[2], ccr5.B.prob[1]))
  dat$attr$ccr5[newIds[newH..wa]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.H..wa, replace = TRUE,
                                        prob = c(1 - sum(ccr5.H.prob),
                                                 ccr5.H.prob[2], ccr5.H.prob[1]))
  dat$attr$ccr5[newIds[newO..wa]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.O..wa, replace = TRUE,
                                        prob = c(1 - sum(ccr5.O.prob),
                                                  ccr5.O.prob[2], ccr5.O.prob[1]))

  # Degree
  dat$attr$deg.main[newIds] <- 0
  dat$attr$deg.pers[newIds] <- 0

  # One-off risk group (all enter at age 18 so are in the "young" group)
  dat$attr$riskg[newIds] <- sample(c("Y1", "Y2", "Y3", "Y4"), nBirths, TRUE)

  # UAI group
  p1 <- dat$param$cond.pers.always.prob
  p2 <- dat$param$cond.inst.always.prob
  rho <- dat$param$cond.always.prob.corr
  uai.always <- bindata::rmvbin(nBirths, c(p1, p2), bincorr = (1 - rho) * diag(2) + rho)
  dat$attr$cond.always.pers[newIds] <- uai.always[, 1]
  dat$attr$cond.always.inst[newIds] <- uai.always[, 2]

  # PrEP
  dat$attr$prepStat[newIds] <- 0

  return(dat)
}



#' @export
#' @rdname births_msm
births_het <- function(dat, at) {

  # Variables
  b.rate.method <- dat$param$b.rate.method
  b.rate <- dat$param$b.rate
  active <- dat$attr$active


  # Process
  nBirths <- 0
  if (b.rate.method == "stgrowth") {
    exptPopSize <- dat$epi$num[1] * (1 + b.rate*at)
    numNeeded <- exptPopSize - sum(active == 1)
    if (numNeeded > 0) {
      nBirths <- rpois(1, numNeeded)
    }
  }
  if (b.rate.method == "totpop") {
    nElig <- dat$epi$num[at - 1]
    if (nElig > 0) {
      nBirths <- rpois(1, nElig * b.rate)
    }
  }
  if (b.rate.method == "fpop") {
    nElig <- dat$epi$num.feml[at - 1]
    if (nElig > 0) {
      nBirths <- rpois(1, nElig * b.rate)
    }
  }


  # Update Population Structure
  if (nBirths > 0) {
    dat <- setBirthAttr_het(dat, at, nBirths)
    dat$el[[1]] <- tergmLite::add_vertices(dat$el[[1]], nBirths)
  }

  if (unique(sapply(dat$attr, length)) != attributes(dat$el[[1]])$n) {
    stop("mismatch between el and attr length in births mod")
  }

  # Output
  dat$epi$b.flow[at] <- nBirths

  return(dat)
}


setBirthAttr_het <- function(dat, at, nBirths) {

  # Set attributes for new births to NA
  dat$attr <- lapply(dat$attr, function(x) c(x, rep(NA, nBirths)))
  newIds <- which(is.na(dat$attr$active))


  # Network Status
  dat$attr$active[newIds] <- rep(1, nBirths)
  dat$attr$entTime[newIds] <- rep(at, nBirths)


  # Demography
  prop.male <- ifelse(is.null(dat$param$b.propmale),
                      dat$epi$propMale[1],
                      dat$param$b.propmale)
  dat$attr$male[newIds] <- rbinom(nBirths, 1, prop.male)

  dat$attr$age[newIds] <- rep(18, nBirths)

  # Circumcision
  entTime <- dat$attr$entTime

  idsNewMale <- which(dat$attr$male == 1 & entTime == at)

  if (length(idsNewMale) > 0) {
    age <- dat$attr$age[idsNewMale]
    newCirc <- rbinom(length(idsNewMale), 1, dat$param$circ.prob.birth)
    isCirc <- which(newCirc == 1)

    newCircTime <- rep(NA, length(idsNewMale))
    newCircTime[isCirc] <- round(-age[isCirc] * (365 / dat$param$time.unit))

    dat$attr$circStat[idsNewMale] <- newCirc
    dat$attr$circTime[idsNewMale] <- newCircTime
  }


  # Epi/Clinical
  dat$attr$status[newIds] <- rep(0, nBirths)

  if (length(unique(sapply(dat$attr, length))) != 1) {
    sapply(dat$attr, length)
    stop("Attribute dimensions not unique")
  }

  return(dat)
}
