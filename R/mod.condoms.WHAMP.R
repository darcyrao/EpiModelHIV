
#' @title Condom Use Module
#'
#' @description Module function stochastically simulates potential condom use
#'              for each act on the discordant edgelist for the WHAMP model
#'
#' @inheritParams aging_msm
#'
#' @details
#' For each act on the discordant edgelist, condom use is stochastically simulated
#' based on the partnership type and, for main partnerships, age combination. Other
#' modifiers for the probability of condom use in that pair are diagnosis of
#' disease, disclosure of status, full or partial HIV viral suppression
#' given HIV anti-retroviral therapy, and PrEP use.
#'
#' @return
#' Updates the discordant edgelist with a \code{uai} variable indicating whether
#' condoms were used in that act.
#'
#' @keywords module msm
#' @export
#'
condoms_msm_whamp <- function(dat, at) {

  # Attributes
  uid <- dat$attr$uid
  diag.status <- dat$attr$diag.status
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  # Parameters
  rcomp.prob <- dat$param$rcomp.prob
  rcomp.adh.groups <- dat$param$rcomp.adh.groups
  rcomp.main.only <- dat$param$rcomp.main.only
  rcomp.discl.only <- dat$param$rcomp.discl.only

  el <- dat$temp$el

  for (type in c("main", "pers", "inst")) {

    ## Variables ##

    # Parameters
    cond.rr <- dat$param$cond.rr

    if (type == "main") {
      el.main <- el[el[, "ptype"] == 1, ]
      age.p1 <- ifelse(age[el.main[, 1]] %in% c(18:34), "Y", "O")
      age.p2 <- ifelse(age[el.main[, 2]] %in% c(18:34), "Y", "O")
      num.YY <- (age.p1 == "Y") + (age.p1 == "Y")
      cond.prob <- (num.YY == 2) * (dat$param$cond.main.YY.prob) +
        (num.B == 1) * (dat$param$cond.main.other.prob) +
        (num.B == 0) * (dat$param$cond.main.other.prob)
      diag.beta <- dat$param$cond.diag.main.beta
      discl.beta <- dat$param$cond.discl.main.beta
      cond.always <- NULL
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
    cond.prob.adj <- cond.prob * cond.rr


    # Transform to UAI logit
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
      uai.prob <- ifelse(ca1 == 1 | ca2 == 1, 0, uai.prob) # If either partner is an always condom user, uai = 0
      if (type == "pers") {
        dat$epi$cprob.always.pers <- NULL
        # dat$epi$cprob.always.pers[at] <- mean(uai.prob == 0)
      } else {
        dat$epi$cprob.always.inst <- NULL
        # dat$epi$cprob.always.inst[at] <- mean(uai.prob == 0)
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

    ai.vec <- elt[, "ai"]
    p1 <- rep(elt[, "p1"], ai.vec)
    p2 <- rep(elt[, "p2"], ai.vec)
    ptype <- rep(elt[, "ptype"], ai.vec)

    uai.prob.peract <- rep(uai.prob, ai.vec)
    uai <- rbinom(length(p1), 1, uai.prob.peract)

    if (type == "main") {
      pid <- rep(1:length(ai.vec), ai.vec)
      al <- cbind(p1, p2, ptype, uai, pid)
    } else {
      pid <- rep(max(al[, "pid"]) + (1:length(ai.vec)), ai.vec)
      tmp.al <- cbind(p1, p2, ptype, uai, pid)
      al <- rbind(al, tmp.al)
    }

  } # end ptype loop

  dat$temp$al <- al

  if (at == 2) {
    dat$epi$ai.events <- rep(NA, 2)
    dat$epi$uai.events <- rep(NA, 2)
  }
  dat$epi$ai.events[at] <- nrow(al)
  dat$epi$uai.events[at] <- sum(al[, "uai"])

  return(dat)
}