
#' @title Condom Use Module
#'
#' @description Module function stochastically simulates potential condom use
#'              for each act on the discordant edgelist for the WHAMP model.
#'
#' @inheritParams aging_msm
#'
#' @details
#' For each act on the discordant edgelist, condom use is stochastically simulated
#' based on the partnership type and, for main partnerships, age combination. Other
#' modifiers for the probability of condom use in that pair are diagnosis of
#' disease, disclosure of status, full or partial HIV viral suppression
#' given HIV anti-retroviral therapy, and PrEP use. If \code{rcomp.discont} is TRUE,
#' condom use stays at on-PrEP risk compensation levels after spontaneous 
#' discontinuation.
#'
#' @return
#' Updates the discordant edgelist with a \code{uai} variable indicating whether
#' condoms were used in that act.
#'
#' @keywords module msm
#' 
#' @export

condoms_msm_whamp <- function(dat, at) {

  # Attributes
  uid <- dat$attr$uid
  age <- dat$attr$age
  diag.status <- dat$attr$diag.status
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  spontDisc <- dat$attr$spontDisc

  # Parameters
  rcomp.prob <- dat$param$rcomp.prob
  if(dat$param$rcomp.adh.groups %in% "all"){
    rcomp.adh.groups <- c(1:3)
  } else if (dat$param$rcomp.adh.groups %in% "top2"){
    rcomp.adh.groups <- c(2:3)
  } else if (dat$param$rcomp.adh.groups %in% "top1"){
    rcomp.adh.groups <- 3
  }
  rcomp.main.only <- dat$param$rcomp.main.only
  rcomp.discl.only <- dat$param$rcomp.discl.only
  rcomp.discont <- dat$param$rcomp.discont

  el <- dat$temp$el

  for (type in c("main", "pers", "inst")) {

    ## Variables ##

    # Parameters
    cond.rr <- dat$param$cond.rr

    if (type == "main") {
      el.main <- el[el[, "ptype"] == 1, ]
      age.p1 <- ifelse(age[el.main[, 1]] <35, "Y", "O")
      age.p2 <- ifelse(age[el.main[, 2]] <35, "Y", "O")
      num.YY <- (age.p1 == "Y") + (age.p2 == "Y")
      cond.prob <- (num.YY == 2) * (dat$param$cond.main.YY.prob) +
        (num.YY == 1) * (dat$param$cond.main.other.prob) +
        (num.YY == 0) * (dat$param$cond.main.other.prob)
      cond.always <- NULL
      diag.beta <- dat$param$cond.diag.main.beta
      discl.beta <- dat$param$cond.discl.main.beta
      ptype <- 1
    }
    if (type == "pers") {
      cond.prob <- dat$param$cond.pers.prob
      diag.beta <- dat$param$cond.diag.pers.beta
      discl.beta <- dat$param$cond.discl.pers.beta
      cond.always <- dat$attr$cond.always.pers
      ptype <- 2
    }
    if (type == "inst") {
      cond.prob <- dat$param$cond.inst.prob
      diag.beta <- dat$param$cond.diag.inst.beta
      discl.beta <- dat$param$cond.discl.inst.beta
      cond.always <- dat$attr$cond.always.inst
      ptype <- 3
    }

    elt <- el[el[, "ptype"] == ptype, ]

    ## Process ##

    # Base condom probs
    cond.prob.base <- cond.prob * cond.rr
    
    if (type %in% c("pers", "inst")) {
    cond.prob.adj <- rep(cond.prob.base, nrow(elt))
    } else {
      cond.prob.adj <- cond.prob.base
    }
    
    # Transform base condom probs to UAI logit
    uai.prob <- 1 - cond.prob.adj
    uai.logodds <- log(uai.prob / (1 - uai.prob))

    # Diagnosis modifier
    pos.diag <- diag.status[elt[, 1]]
    isDx <- which(pos.diag == 1)
    uai.logodds[isDx] <- uai.logodds[isDx] + diag.beta

    # Disclosure modifier
    isDiscord <- which((elt[, "st1"] - elt[, "st2"]) == 1)
    delt <- elt[isDiscord, ]
    discl.list <- dat$temp$discl.list
    disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
    delt.cdl <- uid[delt[, 1]] * 1e7 + uid[delt[, 2]]
    discl.disc <- (delt.cdl %in% disclose.cdl)

    discl <- rep(NA, nrow(elt))
    discl[isDiscord] <- discl.disc

    isDisc <- which(discl == 1)
    uai.logodds[isDisc] <- uai.logodds[isDisc] + discl.beta

    # Back transform to prob
    old.uai.prob <- uai.prob
    uai.prob <- exp(uai.logodds) / (1 + exp(uai.logodds))

    uai.prob[is.na(uai.prob) & old.uai.prob == 0] <- 0
    uai.prob[is.na(uai.prob) & old.uai.prob == 1] <- 1

    # UAI group
    if (type %in% c("pers", "inst")) {
      ca1 <- cond.always[elt[, 1]]
      ca2 <- cond.always[elt[, 2]]
      uai.prob <- ifelse(ca1 == 1 | ca2 == 1, 0, uai.prob) # If either partner is an always condom user, uai.prob = 0
      if (type == "pers") {
        dat$epi$cprob.always.pers <- NULL
        # dat$epi$cprob.always.pers[at] <- mean(cond.prob.adj == 1)
      } else {
        dat$epi$cprob.always.inst <- NULL
        # dat$epi$cprob.always.inst[at] <- mean(cond.prob.adj == 1)
      }
    }
    
    
    # PrEP Status (risk compensation)
    if (rcomp.prob > 0) {
      
      idsRC <- which((prepStat[elt[, 1]] == 1 & prepClass[elt[, 1]] %in% rcomp.adh.groups) |
                       (prepStat[elt[, 2]] == 1 & prepClass[elt[, 2]] %in% rcomp.adh.groups))

      if (rcomp.main.only == TRUE & ptype > 1) {
        idsRC <- NULL
      }
      if (rcomp.discl.only == TRUE) {
        idsRC <- intersect(idsRC, isDisc)
      }
      uai.prob[idsRC] <- 1 - (1 - uai.prob[idsRC]) * (1 - rcomp.prob)
    }
    
    # Continued risk compensation after PrEP discontinuation - only if not diagnosed with HIV and the partner is not an always condom user 
    if (rcomp.discont == TRUE){
      ids.spontDisc <- which((prepStat[elt[, 1]] == 0 & (diag.status[elt[, 1]] == 0 | is.na(diag.status[elt[,1]])) & 
                                prepClass[elt[, 1]] %in% rcomp.adh.groups & spontDisc[elt[, 1]] == 1 & cond.always[elt[, 2]] == 0) |
                               (prepStat[elt[, 2]] == 0 & (diag.status[elt[, 2]] == 0 | is.na(diag.status[elt[, 2]])) & 
                                  prepClass[elt[, 2]] %in% rcomp.adh.groups & spontDisc[elt[, 2]] == 1 & cond.always[elt[, 1]] == 0))
      idsRC.disc <- setdiff(ids.spontDisc, idsRC) # If partner is on PrEP and risk compensates, don't add effect of continued risk compensation from discontinued partner
      uai.prob[idsRC.disc] <- 1 - (1 - uai.prob[idsRC.disc]) * (1 - rcomp.prob)
      
    }
  
    anyprep <- (prepStat[elt[,1]] == 1 | prepStat[elt[,2]] == 1)
    
    ai.vec <- elt[, "ai"]
    p1 <- rep(elt[, "p1"], ai.vec)
    p2 <- rep(elt[, "p2"], ai.vec)
    ptype <- rep(elt[, "ptype"], ai.vec)
    prep <- rep(anyprep, ai.vec)

    uai.prob.peract <- rep(uai.prob, ai.vec)
    uai <- rbinom(length(p1), 1, uai.prob.peract)

    if (type == "main") {
      pid <- rep(1:length(ai.vec), ai.vec)
      al <- cbind(p1, p2, ptype, uai, prep, pid)
    } else {
      pid <- rep(max(al[, "pid"]) + (1:length(ai.vec)), ai.vec)
      tmp.al <- cbind(p1, p2, ptype, uai, prep, pid)
      al <- rbind(al, tmp.al)
    }

  } # end ptype loop

  dat$temp$al <- al

  if (at == 2) {
    dat$epi$ai.events <- rep(NA, 2)
    dat$epi$uai.events <- rep(NA, 2)
    dat$epi$ai.events.main <- rep(NA, 2)
    dat$epi$uai.events.main <- rep(NA, 2)
    dat$epi$ai.events.pers <- rep(NA, 2)
    dat$epi$uai.events.pers <- rep(NA, 2)
    dat$epi$ai.events.inst <- rep(NA, 2)
    dat$epi$uai.events.inst <- rep(NA, 2)
    dat$epi$uai.events.prep <- rep(NA, 2)
    dat$epi$uai.events.noprep <- rep(NA, 2)
  }
  
  dat$epi$ai.events[at] <- nrow(al)
  dat$epi$uai.events[at] <- sum(al[, "uai"])
  dat$epi$ai.events.main[at] <- nrow(al[al[,"ptype"] == 1,])
  dat$epi$uai.events.main[at] <- sum(al[, "uai"][al[,"ptype"] == 1])
  dat$epi$ai.events.pers[at] <- nrow(al[al[,"ptype"] == 2,])
  dat$epi$uai.events.pers[at] <- sum(al[, "uai"][al[,"ptype"] == 2])
  dat$epi$ai.events.inst[at] <- nrow(al[al[,"ptype"] == 3,])
  dat$epi$uai.events.inst[at] <- sum(al[, "uai"][al[,"ptype"] == 3])
  dat$epi$ai.events.prep[at] <- nrow(al[al[,"prep"] == TRUE,])
  dat$epi$uai.events.prep[at] <- sum(al[, "uai"][al[,"prep"] == TRUE])
  dat$epi$ai.events.noprep[at] <- nrow(al[al[,"prep"] == FALSE,])
  dat$epi$uai.events.noprep[at] <- sum(al[, "uai"][al[,"prep"] == FALSE])
  
  return(dat)
}
