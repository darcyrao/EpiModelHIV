
#' @title Sexual Acts Module
#'
#' @description Module function for setting the number of sexual acts on the
#'              discordant edgelist for the WHAMP model.
#'
#' @inheritParams aging_msm
#'
#' @details
#' The number of acts at each time step is specified as a function of the expected values 
#' specific to each partnership type. For one-off partnerships, this is deterministically
#' set at 1, whereas for main and persistent partnerships it is a stochastic draw
#' from a Poisson distribution. The number of total acts may further be modified
#' by the level of HIV viral suppression in an infected person.
#'
#' @return
#' This function returns the \code{dat} object with the updated discordant act
#' list (\code{dal}). Each element of \code{dal} is a data frame with the ids of the
#' discordant pair repeated the number of times they have AI.
#'
#' @keywords module msm
#' @export
#'
acts_msm_whamp <- function(dat, at) {

  for (type in c("main", "pers", "inst")) {

    ## Variables ##

    # Attributes
    status <- dat$attr$status

    # Parameters
    ai.scale <- dat$param$ai.scale
    if (type == "main") {
      base.ai.rate <- dat$param$base.ai.main.rate
      fixed <- FALSE
      ptype <- 1
      el <- dat$el[[1]]
    }
    if (type == "pers") {
      base.ai.rate <- dat$param$base.ai.pers.rate
      fixed <- FALSE
      ptype <- 2
      el <- dat$el[[2]]
    }
    if (type == "inst") {
      base.ai.rate <- 1
      fixed <- ifelse(ai.scale != 1, FALSE, TRUE)
      ptype <- 3
      el <- dat$el[[3]]
    }

    ## Processes ##

    # Construct edgelist

    st1 <- status[el[, 1]]
    st2 <- status[el[, 2]]
    disc <- abs(st1 - st2) == 1
    el[which(disc == 1 & st2 == 1), ] <- el[which(disc == 1 & st2 == 1), 2:1] # If discordant, pos partners in col 1
    el <- cbind(el, status[el[, 1]], status[el[, 2]])
    colnames(el) <- c("p1", "p2", "st1", "st2")

    if (nrow(el) > 0) {

      # Base AI rates
      ai.rate <- base.ai.rate * ai.scale

      # ## STI associated cessation of activity
      # idsCease <- which(dat$attr$GC.cease == 1 | dat$attr$CT.cease == 1)
      # noActs <- el[, "p1"] %in% idsCease | el[, "p2"] %in% idsCease
      # ai.rate[noActs] <- 0

      # Final act number
      if (fixed == FALSE) {
        ai <- rpois(nrow(el), ai.rate)
      } else {
        ai <- round(ai.rate)
      }

      # Full edge list
      el <- cbind(el, ptype, ai)
      colnames(el)[5:6] <- c("ptype", "ai")

      if (type == "main") {
        dat$temp$el <- el
      } else {
        dat$temp$el <- rbind(dat$temp$el, el)
      }
    }

  } # loop over type end

  # Remove inactive edges from el
  dat$temp$el <- dat$temp$el[-which(dat$temp$el[, "ai"] == 0), ]

  return(dat)
}
