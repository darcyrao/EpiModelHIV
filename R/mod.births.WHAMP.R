
#' @title Births Module
#'
#' @description Module function for births or entries into the sexually active
#'              population for the WHAMP model
#'
#' @inheritParams aging_msm
#'
#' @details
#' New population members are added based on expected numbers of entries,
#' stochastically determined with draws from Poisson distributions.
#'For each new entry, a set of attributes is added for that node,
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
  
  nBirths.gen <- dat$epi$dth.gen..wa[at] # Births to replace those who died from general mortality
  nBirths.age <- dat$epi$dth.age[at] # Births to replace those who aged out
  nBirths <- nBirths.gen + nBirths.age

  ## Update Attr
  
  if (nBirths > 0) {
    dat <- setBirthAttr_msm(dat, at, nBirths)
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
  prop.H.KC <- 0.0549
  prop.B.KC <- 0.0421
  prop.O.KC <- 0.4739
  prop.H.OW <- 0.0309
  prop.B.OW <- 0.0166
  prop.O.OW <- 0.2807
  prop.H.EW <- 0.0222
  prop.B.EW <- 0.0021
  prop.O.EW <- 1 - sum(prop.H.KC, prop.B.KC, prop.O.KC, prop.H.OW, prop.B.OW, prop.O.OW, prop.H.EW, prop.B.EW) #set value for last group to remainder to avoid rounding error
  
  race <- sample(rep(c("B", "W"), c(0.5*nBirths, 0.5*nBirths))) #-- Assign old race attribute until finish revising code
  new.B <- which(race == "B") #-- Delete when finish de-bugging
  new.W <- which(race == "W") #-- Delete when finish de-bugging
  dat$attr$race[newIds] <- race #-- Delete when finish de-bugging
  race.region <- sample(rep(c("H.KC", "B.KC", "O.KC", "H.OW", "B.OW", "O.OW", "H.EW", "B.EW", "O.EW"), 
                            c(prop.H.KC*nBirths, prop.B.KC*nBirths, prop.O.KC*nBirths, 
                              prop.H.OW*nBirths, prop.B.OW*nBirths, prop.O.OW*nBirths, 
                              prop.H.EW*nBirths, prop.B.EW*nBirths, prop.O.EW*nBirths)))
  race..wa <- ifelse(race.region %in% c("H.KC", "H.OW", "H.EW"), "H",
                     ifelse(race.region %in% c("B.KC", "B.OW", "B.EW"), "B",
                            ifelse(race.region %in% c("O.KC", "O.OW", "O.EW"), "O", NA))) 
  region <- ifelse(race.region %in% c("H.KC", "B.KC", "O.KC"), "KC",
                   ifelse(race.region %in% c("H.OW", "B.OW", "O.OW"), "OW",
                          ifelse(race.region %in% c("H.EW", "B.EW", "O.EW"), "EW", NA)))
  
  newH..wa <- which(race..wa == "H")
  newB..wa <- which(race..wa == "B")
  newO..wa <- which(race..wa == "O")
  new.KC <- which(region == "KC")
  new.OW <- which(region == "OW")
  new.EW <- which(region == "EW")
  
  dat$attr$race..wa[newIds] <- race..wa
  dat$attr$region[newIds] <- region
  
  # New entries will enter at age 18
  dat$attr$age[newIds] <- rep(dat$param$birth.age, nBirths)
  dat$attr$sqrt.age[newIds] <- sqrt(dat$attr$age[newIds])

  # Disease status and related
  dat$attr$status[newIds] <- rep(0, nBirths)

  dat$attr$inst.ai.class[newIds] <- sample(1:dat$param$num.inst.ai.classes,
                                           nBirths, replace = TRUE)

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
  ccr5.W.prob <- dat$param$ccr5.W.prob
  dat$attr$ccr5[newIds[newB]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.B, replace = TRUE,
                                        prob = c(1 - sum(ccr5.B.prob),
                                                 ccr5.B.prob[2], ccr5.B.prob[1]))
  dat$attr$ccr5[newIds[newW]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.W, replace = TRUE,
                                        prob = c(1 - sum(ccr5.W.prob),
                                                 ccr5.W.prob[2], ccr5.W.prob[1]))


  # Degree
  dat$attr$deg.main[newIds] <- 0
  dat$attr$deg.pers[newIds] <- 0

  # One-off risk group
  dat$attr$riskg[newIds] <- sample(1:5, nBirths, TRUE)

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
