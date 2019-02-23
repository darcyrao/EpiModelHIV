
#' @title Risk History Module
#'
#' @description Module function to track the risk history of uninfected persons
#'              for purpose of intervention targeting for the WHAMP model.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
riskhist_msm_whamp <- function(dat, at) {

  if (at < dat$param$riskh.start) {
    return(dat)
  }

  ## Attributes
  uid <- dat$attr$uid
  dx <- dat$attr$diag.status
  since.test <- at - dat$attr$last.neg.test
  # rGC.tx <- dat$attr$rGC.tx
  # uGC.tx <- dat$attr$uGC.tx
  # rCT.tx <- dat$attr$rCT.tx
  # uCT.tx <- dat$attr$uCT.tx

  ## Parameters
  time.unit <- dat$param$time.unit

  ## Edgelist, adds uai summation per partnership from act list
  pid <- NULL # For R CMD Check
  al <- as.data.frame(dat$temp$al)
  by_pid <- group_by(al, pid)
  uai <- summarise(by_pid, uai = sum(uai))[, 2]
  el <- as.data.frame(cbind(dat$temp$el, uai))

  # Remove concordant positive edges
  el2 <- el[el$st2 == 0, ]

  # Initialize attributes
  if (is.null(dat$attr$prep.ind.discord.ongoing)) {
    dat$attr$prep.ind.discord.ongoing <- rep(NA, length(uid))
    dat$attr$prep.ind.uai.risk <- rep(NA, length(uid))
  }

  ## Degree ##
  main.deg <- get_degree(dat$el[[1]])
  casl.deg <- get_degree(dat$el[[2]])
  inst.deg <- get_degree(dat$el[[3]])


  ## Preconditions ##

  # Any UAI
  uai.any <- unique(c(el2$p1[el2$uai > 0],
                      el2$p2[el2$uai > 0]))

  # Monogamous partnerships: 2-sided
  tot.deg <- main.deg + casl.deg + inst.deg
  uai.mono <- intersect(which(tot.deg == 1), uai.any)
  
  # "Negative" partnerships
  tneg <- unique(c(el2$p1[el2$st1 == 0], el2$p2[el2$st1 == 0])) ##-- should'nt this be el2$p2[el2$st2] ==0] ?
  fneg <- unique(c(el2$p1[which(dx[el2$p1] == 0)], el2$p2[which(dx[el2$p1] == 0)])) ##-- should'nt this be el2$p2[which(dx[el2$p2] ==0] ?
  all.neg <- c(tneg, fneg)

  ## Condition 1: ongoing positive partner who has disclosed
  discord <- el2[el2$st1 == 1 & el2$ptype %in% 1:2, ]
  
  # Disclosure
  discl.list <- dat$temp$discl.list
  disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
  delt.cdl <- uid[discord[, 1]] * 1e7 + uid[discord[, 2]]
  discl <- (delt.cdl %in% disclose.cdl)
  ai.sd <- discord$p2[discl == TRUE]
  
  dat$attr$prep.ind.discord.ongoing[ai.sd] <- at
  
  ## Condition 2: UAI outside of a 2-sided "monogamous" partnership,
  ##               with a partner tested negative in past 6 months
  mono.neg <- intersect(uai.mono, all.neg)
  part.id1 <- c(el2[el2$p1 %in% mono.neg, 2], el2[el2$p2 %in% mono.neg, 1])
  part.recently.tested <- since.test[part.id1] <= (180/time.unit)
  mono.neg.recently.tested <- mono.neg[which(part.recently.tested == TRUE)]
  
  uai.risk <- setdiff(uai.any, mono.neg.recently.tested)
  dat$attr$prep.ind.uai.risk[uai.risk] <- at

  return(dat)
}
