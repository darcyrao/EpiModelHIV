
#' @title Update Role Class in One-Off Partnerships
#'
#' @description Module function for updating act class in instantaneous partnerships
#'              at age 50 for the WHAMP model.
#'
#' @inheritParams aging_msm
#'
#' @return
#' This function updates the individual-level attribute \code{riskg} on
#' \code{dat$attr}.
#'
#' @keywords module msm
#' 
#' @export
#'
update_aiclass_msm_whamp <- function(dat, at) {

  age <- dat$attr$age
  old.riskg <- dat$attr$riskg

  dat$attr$riskg[age >= 50 & old.riskg %in% "Y1"] <- "O1"
  dat$attr$riskg[age >= 50 & old.riskg %in% "Y2"] <- "O2"
  dat$attr$riskg[age >= 50 & old.riskg %in% "Y3"] <- "O3"
  dat$attr$riskg[age >= 50 & old.riskg %in% "Y4"] <- "O4"
  
  return(dat)
}
