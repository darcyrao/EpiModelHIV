

# MSM -----------------------------------------------------------------

#' @title Calculate Target Statistics for Network Model Estimation for the WHAMP model
#'
#' @description Calculates the target statistics for the formation and dissolution
#'              components of the network model to be estimated with \code{netest}. 
#'              For the WHAMP model, the base EpiModelHIV code was modified to define 
#'              three racial/ethnic groups (Hispanic, black, and other), age structure,
#'              add region as an attribute, and change the way risk group quantiles are 
#'              defined.
#'
#' @param time.unit Time unit relative to 1 for daily.
#' @param num.B Population size of black MSM.
#' @param num.W Population size of white MSM.
#' @param num.H..wa Population size of Hispanic MSM.
#' @param num.B..wa Population size of non-Hispanic black MSM.
#' @param num.O..wa Population size of non-Hispanic other MSM.
#' @param num.KC Population size of MSM in other counties in western Washington.
#' @param num.OW Population size of MSM in other counties in western Washington.
#' @param num.EW Population size of MSM in other counties in western Washington.
#' @param agestr Vector with the proportion of MSM in each age group 18-24, 25-29...50-59
#' @param deg.mp Degree distribution matrix for main and casual partners overall,
#'        as a 2 by 3 matrix.
#' @param deg.mp.H Degree distribution matrix for main and causal partners for
#'        Hispanic MSM, as a 2 by 3 matrix.
#' @param deg.mp.B Degree distribution matrix for main and causal partners for
#'        Black MSM, as a 2 by 3 matrix. 
#' @param deg.mp.O Degree distribution matrix for main and causal partners for
#'        Other race/ethnicity MSM, as a 2 by 3 matrix.     
#' @param deg.m.region Main degree distribution by region (EW, KC, OW)
#' @param prop.deg.p.EW Proportion of men with persistent degree 0, 1, and 2+ who live in eastern WA
#' @param prop.deg.p.KC Proportion of men with persistent degree 0, 1, and 2+ who live in King County
#' @param prop.deg.p.OW Proportion of men with persistent degree 0, 1, and 2+ who live in western WA                
#' @param mdeg.inst Mean degree, or rate, of one-off partnerships per day.
#' @param qnts.18to49 Means rates withing quartiles of the distribution of instantaneous partnerships 
#'        for MSM ages 18-49. Use \code{NA} to ignore these quantiles in the target statistics.
#' @param qnts.50to59 Means rates withing quartiles of the distribution of instantaneous partnerships
#'        for MSM ages 50-59. Use \code{NA} to ignore these quantiles in the target statistics.
#' @param inst.bho Mean rate of instantaneous partnerships by racial/ethnic group
#' @param inst.region Distribution of instantaneous partnerships by region (EW, KC, OW)
#' @param prop.hom.mpi.H A vector of length 3 for the proportion of main, casual,
#'        and one-off partnerships in same race for Hispanic MSM.
#' @param prop.hom.mpi.B A vector of length 3 for the proportion of main, casual,
#'        and instantaneous partnerships in same race for black MSM.
#' @param prop.hom.mpi.O A vector of length 3 for the proportion of main, casual,
#'        and instantaneous partnerships in same race for other race/ethnicity MSM.
#' @param sqrt.adiff.mpi Vector of length 3 with the mean absolute differences
#'        in the square root of ages in main, casual, and instantaneous partnerships.
#' @param prop.hom.region.mpi A vector of length 3 for the proportion of main, casual,
#'        and instantaneous partnerships that are within-region.
#' @param diss.main Dissolution model formula for main partnerships.
#' @param diss.pers Dissolution model formula for casual partnerships.
#' @param durs.main Duration of main partnerships in days.
#' @param durs.pers Duration of persistent partnerships in days. 
#' @param ages Integer vector of ages in years that defines range of possible
#'        initial ages in the population.
#' @param asmr.B Vector of length 40 defining the age-specific
#'        mortality rate for persons within that age slot, for black MSM.
#' @param asmr.W Vector of length 40 defining the age-specific
#'        mortality rate for persons within that age slot, for white MSM.
#' @param asmr.H..wa Vector of length 60 defining the age-specific
#'        mortality rates for for Hispanic MSM.
#' @param asmr.B..wa Vector of length 60 defining the age-specific
#'        mortality rates for for black MSM.
#' @param asmr.O..wa Vector of length 60 defining the age-specific
#'        mortality rates for for other race/ethnicity MSM.
#' @param role.prob Vector of length 3 for the probability of sexual role as
#'        insertive, receptive, and versatile.
#'
#' @details
#' This function performs basic calculations to determine the components of the
#' formation and dissolution models for the network model estimation to be
#' conducted with \code{\link{netest}}. The inputs inputs for this function are
#' calculated externally to the package in a setup scenario file.
#'
#' @keywords msm
#'
#' @seealso
#' Network statistics calculated here are entered into \code{\link{base_nw_msm}}
#' to construct the base network, and then into the parameters in
#' \code{\link{param_msm}}.
#'
#' @export
#'
calc_nwstats_msm_whamp <- function(time.unit = 7,
                             num.B,
                             num.W,
                             num.H.KC,
                             num.B.KC,
                             num.O.KC,
                             num.H.OW,
                             num.B.OW,
                             num.O.OW,
                             num.H.EW,
                             num.B.EW,
                             num.O.EW,
                             agestr,
                             deg.mp,
                             deg.mp.H,
                             deg.mp.B,
                             deg.mp.O,
                             deg.m.region,
                             prop.deg.p.EW,
                             prop.deg.p.KC,
                             prop.deg.p.OW,
                             mdeg.inst,
                             qnts.18to49,
                             qnts.50to59,
                             inst.bho,
                             inst.region,
                             prop.hom.mpi.H,
                             prop.hom.mpi.B,
                             prop.hom.mpi.O,
                             sqrt.adiff.mpi,
                             prop.hom.region.mpi,
                             diss.main,
                             diss.pers,
                             durs.main,
                             durs.pers,
                             ages,
                             asmr.B,
                             asmr.W,
                             asmr.H..wa,
                             asmr.B..wa,
                             asmr.O..wa,
                             role.prob) {

  if (sum(deg.mp) != 1) {
    stop("deg.mp must sum to 1")
  }
  if (sum(deg.mp.H) != 1 | sum(deg.mp.B) !=1 | sum(deg.mp.O) !=1) {
    stop("deg.mp terms must sum to 1")
  }
  if (sum(deg.m.region) !=1) {
    stop("deg distribution must sum to 1.")
  }
  if (sum(inst.region) !=1) {
    stop("deg distribution must sum to 1.")
  }
  
  if (sum(agestr) !=1) {
    stop("agestr must sum to 1")
  }

  # total number and numbers by race/ethnicity, region, and age
  num <- num.H.KC + num.B.KC + num.O.KC + num.H.OW + num.B.OW + num.O.OW + num.H.EW + num.B.EW + num.O.EW
  num.H..wa <- sum(num.H.KC, num.H.OW, num.H.EW)
  num.B..wa <- sum(num.B.KC, num.B.OW, num.B.EW)
  num.O..wa <- sum(num.O.KC, num.O.OW, num.O.EW)

  num.KC <- sum(num.H.KC, num.B.KC, num.O.KC)
  num.OW <- sum(num.H.OW, num.B.OW, num.O.OW)
  num.EW <- sum(num.H.EW, num.B.EW, num.O.EW)
  
  # Degree by region
    mean.degm.EW <- ((sum(deg.mp[2,])*num)*deg.m.region[1])/num.EW
    mean.degm.KC <- ((sum(deg.mp[2,])*num)*deg.m.region[2])/num.KC
    mean.degm.OW <- ((sum(deg.mp[2,])*num)*deg.m.region[3])/num.OW
    
    deg.p.EW <- ((colSums(deg.mp)*num)*prop.deg.p.EW)/num.EW
    deg.p.KC <- ((colSums(deg.mp)*num)*prop.deg.p.KC)/num.KC
    deg.p.OW <- ((colSums(deg.mp)*num)*prop.deg.p.OW)/num.OW
    
    
    
  # Main partnerships -------------------------------------------------------

  # Persons in partnerships by casual degree
    totdeg.m.by.dp <- c(num * deg.mp[2, ])

  # Persons in partnerships by race/ethnicity (B, H, O)
    totdeg.m.by.race <- c(sum(num.B..wa * deg.mp.B[2,]),
                          sum(num.H..wa * deg.mp.H[2,]),
                          sum(num.O..wa * deg.mp.O[2,]))
  
  # Persons in partnerships by region (EW, KC, OW)
    totdeg.m.by.region <- c(num.EW * mean.degm.EW,
                            num.KC * mean.degm.KC,
                            num.OW * mean.degm.OW)

  # Number of partnerships
  edges.m <- (sum(totdeg.m.by.dp)) / 2

  # Race mixing: nodematch target stat: number of BB, HH, and OO partnerships
    edges.hom.m.hbo <- c((totdeg.m.by.race[1]*prop.hom.mpi.B[1] / 2), 
                     (totdeg.m.by.race[2]*prop.hom.mpi.H[1] / 2),
                     (totdeg.m.by.race[3]*prop.hom.mpi.O[1] / 2))
    
  # Sqrt absdiff term for age
  sqrt.adiff.m <- edges.m * sqrt.adiff.mpi[1]

  # Compile target stats: edges, nodefactor("deg.pers"), nodefactor("race"), nodefactor("region"), nodematch("race"), absdiff("sqrt.age")
          ##-- confirm which level to omit for nodefactor terms, confirm how specify target stats for offset terms
  stats.m <- c(edges.m, totdeg.m.by.dp[2:3], totdeg.m.by.race[1:2], totdeg.m.by.region[c(1,3)], edges.hom.m.hbo, sqrt.adiff.m)

  # Dissolution model
   ## Expected mortality is a weighted avg of the racial/ethnic- and age-specific mortality ratios
   exp.mort.H <- mean(asmr.H..wa[18:24])*agestr[1] + mean(asmr.H..wa[25:29])*agestr[2] + 
     mean(asmr.H..wa[30:34])*agestr[3] + mean(asmr.H..wa[35:39])*agestr[4] + 
     mean(asmr.H..wa[40:44])*agestr[5] + mean(asmr.H..wa[45:49])*agestr[6] +
     mean(asmr.H..wa[50:54])*agestr[7] + mean(asmr.H..wa[55:59])*agestr[8]
   
   exp.mort.B <- mean(asmr.B..wa[18:24])*agestr[1] + mean(asmr.B..wa[25:29])*agestr[2] + 
     mean(asmr.B..wa[30:34])*agestr[3] + mean(asmr.B..wa[35:39])*agestr[4] + 
     mean(asmr.B..wa[40:44])*agestr[5] + mean(asmr.B..wa[45:49])*agestr[6] +
     mean(asmr.B..wa[50:54])*agestr[7] + mean(asmr.B..wa[55:59])*agestr[8]
   
   exp.mort.O <- mean(asmr.O..wa[18:24])*agestr[1] + mean(asmr.O..wa[25:29])*agestr[2] + 
     mean(asmr.O..wa[30:34])*agestr[3] + mean(asmr.O..wa[35:39])*agestr[4] + 
     mean(asmr.O..wa[40:44])*agestr[5] + mean(asmr.O..wa[45:49])*agestr[6] +
     mean(asmr.O..wa[50:54])*agestr[7] + mean(asmr.O..wa[55:59])*agestr[8]
      
   exp.mort <- exp.mort.H*(num.H..wa / num) + exp.mort.B*(num.B..wa / num) + exp.mort.O*(num.O..wa / num)

  coef.diss.m <- dissolution_coefs(dissolution = diss.main,
                                   duration = durs.main / time.unit,
                                   d.rate = exp.mort)


  # Casual partnerships -----------------------------------------------------

  # Persons in partnerships by main degree
  totdeg.p.by.dm <- c(num * deg.mp[, 2] + num * deg.mp[, 3] * 2)

  # Persons in partnerships by race/ethnicity (B, H, O)
  totdeg.p.by.race <- c(sum(num.B..wa * deg.mp.B[, 2] + num.B..wa * deg.mp.B[, 3] * 2),
                        sum(num.H..wa * deg.mp.H[, 2] + num.H..wa * deg.mp.H[, 3] * 2),
                        sum(num.O..wa * deg.mp.O[, 2] + num.O..wa * deg.mp.O[, 3] * 2))
  
  # Persons in partnerships by region (EW, KC, OW)
  totdeg.p.by.region <- c(sum(num.EW * deg.p.EW[2], num.EW * deg.p.EW[3] *2),
                          sum(num.KC * deg.p.KC[2], num.KC * deg.p.KC[3] *2),
                          sum(num.OW * deg.p.OW[2], num.OW * deg.p.OW[3] *2))
  
  # Persons concurrent
  conc.p <- c(sum(deg.mp[, 3]) * num)
  
  # Number of partnerships
  edges.p <- sum(totdeg.p.by.dm) / 2

  # Race mixing: nodematch target stat: number of BB, HH, and OO partnerships
  edges.hom.p.hbo <- c((totdeg.p.by.race[1]*prop.hom.mpi.B[2] / 2), 
                   (totdeg.p.by.race[2]*prop.hom.mpi.H[2] / 2),
                   (totdeg.p.by.race[3]*prop.hom.mpi.O[2] / 2))
  
  # Regional mixing: nodematch target stat: number of within-region partnerships
  edges.hom.p.region <- edges.p * prop.hom.region.mpi[2]
    
  # Sqrt absdiff term for age
  sqrt.adiff.p <- edges.p * sqrt.adiff.mpi[2]

  # Compile target statistics: edges, nodefactor("deg.main"), concurrent, nodematch("race"), nodematch("region"), absdiff(sqrt.age)
    ##-- confirm which level to omit for nodefactor terms, confirm how specify target stats for offset terms
  stats.p <- c(edges.p, totdeg.p.by.dm[2], totdeg.p.by.race[1:2], totdeg.p.by.region[c(1,3)], conc.p, edges.hom.p.hbo, edges.hom.p.region, sqrt.adiff.p)
    

  # Dissolution model
  coef.diss.p <- dissolution_coefs(dissolution = diss.pers,
                                   duration = durs.pers / time.unit,
                                   d.rate = exp.mort)



  # Instant partnerships ----------------------------------------------------

  # Number of instant partnerships per time step, by main and casl degree
  num.inst <- num * deg.mp * mdeg.inst * time.unit
  
  # Number of instant parnterships by risk quantiles
  num.50to59 <- num*sum(agestr[7:8])
  num.18to49 <- num*sum(agestr[1:6])
  
  if (!is.na(qnts.18to49[1]) & !is.na(qnts.50to59[1])) {
      num.riskg.50to59 <- (0.25*num.50to59) * qnts.50to59 * time.unit
      num.riskg.18to49 <- (0.25*num.18to49) * qnts.18to49 * time.unit
      num.riskg <- c(num.riskg.50to59, num.riskg.18to49)
  }
  
  # Number of instant partnerships per time step by race/ethnicity
  totdeg.i.bho <- c(inst.bho[1]*num.B..wa, 
                    inst.bho[2]*num.H..wa, 
                    inst.bho[3]*num.O..wa)*time.unit
  
  # Number of instant partnerships per time step by region (EW, KC, OW)
  totdeg.i.region <- inst.region*sum(num.inst)
  
  # Total number of instant partnerships per time step
  totdeg.i <- sum(num.inst)

  # Number of partnerships
  edges.i <- sum(totdeg.i) / 2

  # Race mixing: nodematch target stat: number of BB, HH, and OO partnerships
  edges.hom.i.bho <- c((totdeg.i.bho[1]*prop.hom.mpi.B[3] / 2), 
                   (totdeg.i.bho[2]*prop.hom.mpi.H[3] / 2),
                   (totdeg.i.bho[3]*prop.hom.mpi.O[3] / 2))
  
  # Regional mixing: nodematch target stat: number of within-region partnerships
  edges.hom.i.region <- edges.i * prop.hom.region.mpi[3]
  
  # Sqrt absdiff term for age
  sqrt.adiff.i <- edges.i * sqrt.adiff.mpi[3]

  # Compile target stats: edges, nodefactor(c("deg.main", "deg.pers")), nodefactor("riskg") nodematch("race"), nodematch("region"), absdiff(sqrt.age)
    ##-- confirm which level to omit for nodefactor terms, confirm how specify target stats for offset terms
  if (!is.na(qnts.18to49[1]) & !is.na(qnts.50to59[1])) {
    stats.i <- c(edges.i, num.inst[-1], num.riskg[-8], totdeg.i.bho[1:2], totdeg.i.region[c(1,3)], edges.hom.i.bho, edges.hom.i.region, sqrt.adiff.i)
  } else {
    stats.i <- c(edges.i, num.inst[-1], totdeg.i.bho[1:2], totdeg.i.region[c(1,3)], edges.hom.i.hbo, edges.hom.i.region, sqrt.adiff.i)
  }




  # Compile results ---------------------------------------------------------
  out <- list()
  
  out$stats.m <- stats.m
  out$stats.p <- stats.p
  out$stats.i <- stats.i

  out$coef.diss.m <- coef.diss.m
  out$coef.diss.p <- coef.diss.p

  out$ages <- ages
  out$agestr <- agestr
  out$asmr.B <- asmr.B
  out$asmr.W <- asmr.W
  out$asmr.H..wa <- asmr.H..wa
  out$asmr.B..wa <- asmr.B..wa
  out$asmr.O..wa <- asmr.O..wa

  out$time.unit <- time.unit
  
  out$num.B <- num.B
  out$num.W <- num.W
  out$num.H.KC <- num.H.KC
  out$num.B.KC <- num.B.KC
  out$num.O.KC <- num.O.KC
  out$num.H.OW <- num.H.OW
  out$num.B.OW <- num.B.OW
  out$num.O.OW <- num.O.OW
  out$num.H.EW <- num.H.EW
  out$num.B.EW <- num.B.EW
  out$num.O.EW <- num.O.EW
  
  out$deg.mp <- deg.mp
  out$deg.mp.H <- deg.mp.H
  out$deg.mp.B <- deg.mp.B
  out$deg.mp.O <- deg.mp.O
  out$mean.degm.EW <- mean.degm.EW
  out$mean.degm.KC <- mean.degm.KC
  out$mean.degm.OW <- mean.degm.OW
  out$deg.p.EW <- deg.p.EW
  out$deg.p.KC <- deg.p.KC
  out$deg.p.OW <- deg.p.OW
 
  out$role.prob <- role.prob

  class(out) <- "nwstats"
  return(out)
  
}


#' @title Construct Base Network for Model Estimation and Simulation for the WHAMP model
#'
#' @description Initializes the base network for model estimation within
#'              \code{netest}.
#'
#' @param nwstats An object of class \code{nwstats}, as output from
#'        \code{\link{calc_nwstats_msm}}.
#'
#' @details
#' This function takes the output of \code{\link{calc_nwstats_msm}} and constructs
#' an empty network with the necessary attributes for race, region, square root of age,
#' and sexual role class. This base network is used for all three network
#' estimations.
#'
#' @seealso
#' The final vertex attributes on the network for cross-network degree are
#' calculated and set on the network with \code{\link{assign_degree}}.
#'
#' @keywords msm
#' @export
#'
base_nw_msm_whamp <- function(nwstats) {

  num.B <- nwstats$num.B #-- delete when finish debugging
  num.W <- nwstats$num.W #-- delete when finish debugging
  num.H.KC <- nwstats$num.H.KC
  num.B.KC <- nwstats$num.B.KC
  num.O.KC <- nwstats$num.O.KC
  num.H.OW <- nwstats$num.H.OW
  num.B.OW <- nwstats$num.B.OW
  num.O.OW <- nwstats$num.O.OW
  num.H.EW <- nwstats$num.H.EW
  num.B.EW <- nwstats$num.B.EW
  num.O.EW <- nwstats$num.O.EW
  agestr <- nwstats$agestr

  # Initialize network
  n <- num.H.KC + num.B.KC + num.O.KC + num.H.OW + num.B.OW + num.O.OW + num.H.EW + num.B.EW + num.O.EW
  nw <- network::network.initialize(n, directed = FALSE)

  # Calculate attributes
  race <- c(rep("B", num.B), rep("W", num.W)) #-- delete when finish debugging
  race <- sample(race) #-- delete when finish debugging
  
  race.region <- c(rep("H.KC", num.H.KC), rep("B.KC", num.B.KC), rep("O.KC", num.O.KC),  
                   rep("H.OW", num.H.OW), rep("B.OW", num.B.OW), rep("O.OW", num.O.OW), 
                   rep("H.EW", num.H.EW), rep("B.EW", num.B.EW), rep("O.EW", num.O.EW))
  race.region <- sample(race.region)
  
  race..wa <- rep(NA, n)
  race..wa[race.region %in% c("H.KC", "H.OW", "H.EW")] <- "H"
  race..wa[race.region %in% c("B.KC", "B.OW", "B.EW")] <- "B"
  race..wa[race.region %in% c("O.KC", "O.OW", "O.EW")] <- "O"
  
  region <- rep(NA, n)
  region[race.region %in% c("H.KC", "B.KC", "O.KC")] <- "KC"
  region[race.region %in% c("H.OW", "B.OW", "O.OW")] <- "OW"
  region[race.region %in% c("H.EW", "B.EW", "O.EW")] <- "EW"
  
  n.by.age <- table(apportion_lr(n, c("18-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59"), agestr))
  age <- c(sample(seq(18, 24.999, 1 / (365 / nwstats$time.unit)), n.by.age[[1]], TRUE),
           sample(seq(25, 29.999, 1 / (365 / nwstats$time.unit)), n.by.age[[2]], TRUE),
           sample(seq(30, 34.999, 1 / (365 / nwstats$time.unit)), n.by.age[[3]], TRUE),
           sample(seq(35, 39.999, 1 / (365 / nwstats$time.unit)), n.by.age[[4]], TRUE),
           sample(seq(40, 44.999, 1 / (365 / nwstats$time.unit)), n.by.age[[5]], TRUE),
           sample(seq(45, 49.999, 1 / (365 / nwstats$time.unit)), n.by.age[[6]], TRUE),
           sample(seq(50, 54.999, 1 / (365 / nwstats$time.unit)), n.by.age[[7]], TRUE),
           sample(seq(55, 59.999, 1 / (365 / nwstats$time.unit)), n.by.age[[8]], TRUE))
  
  sqrt.age <- sqrt(age)

  role <- sample(apportion_lr(n, c("I", "R", "V"), nwstats$role.prob)) 

  riskg.50to59 <- sample(apportion_lr(n*sum(agestr[7:8]), c("O1", "O2", "O3", "O4"), rep(0.25, 4)))
  riskg.18to49 <- sample(apportion_lr(n*sum(agestr[1:6]), c("Y1", "Y2", "Y3", "Y4"), rep(0.25, 4)))
  riskg <- rep(NA, n)
  riskg[age>=18 & age<50] <- riskg.18to49
  riskg[age>=50] <- riskg.50to59

  attr.names <- c("race", "race..wa", "region", "riskg", "sqrt.age", "role.class")
  attr.values <- list(race, race..wa, region, riskg, sqrt.age, role)
  nw <- network::set.vertex.attribute(nw, attr.names, attr.values)

  return(nw)
}


#' @title Assign Degree Vertex Attribute on Network Objects for the WHAMP model
#'
#' @description Assigns the degree vertex attributes on network objects
#'              conditional on their values from the other networks.
#'
#' @param nw Object of class \code{network} that is the target for the vertex
#'        attribute.
#' @param deg.type Type of degree to assign to \code{nw}, with options of
#'        \code{"pers"} to assign casual degree onto main network and
#'        \code{"main"} to assign main degree to casual network.
#' @param nwstats Object of class \code{nwstats}.
#'
#' @details
#' This function assigns the degree of other networks as a vertex attribute on the
#' target network given a bivariate degree mixing matrix of main, casual, and
#' one-partnerships contained in the \code{nwstats} data.
#'
#' @keywords msm
#' @export
#'
assign_degree_whamp <- function(nw, deg.type, nwstats) {

  if (!("network" %in% class(nw))) {
    stop("nw must be of class network")
  }
  if (!isTRUE(all.equal(sum(colSums(nwstats$deg.mp.H)), 1, tolerance = 5e-6))) {
    stop("H degree distributions do not sum to 1")
  }
  if (!isTRUE(all.equal(sum(colSums(nwstats$deg.mp.B)), 1, tolerance = 5e-6))) {
    stop("B degree distributions do not sum to 1")
  }
  if (!isTRUE(all.equal(sum(colSums(nwstats$deg.mp.O)), 1, tolerance = 5e-6))) {
    stop("O degree distributions do not sum to 1")
  }
  if (!isTRUE(all.equal(sum(colSums(nwstats$deg.mp)), 1, tolerance = 5e-6))) {
    stop("Degree distributions do not sum to 1")
  }
  

  if (deg.type == "main") {
    
    attr.name <- "deg.main"
    
    #Calculate expected mean main degree for each race by region combination assuming independence
    mdeg.main.H.KC <- sum(nwstats$deg.mp.H[2,])*nwstats$mean.degm.KC/sum(nwstats$deg.mp[2,])
    mdeg.main.B.KC <- sum(nwstats$deg.mp.B[2,])*nwstats$mean.degm.KC/sum(nwstats$deg.mp[2,])
    mdeg.main.O.KC <- sum(nwstats$deg.mp.O[2,])*nwstats$mean.degm.KC/sum(nwstats$deg.mp[2,])
    mdeg.main.H.OW <- sum(nwstats$deg.mp.H[2,])*nwstats$mean.degm.OW/sum(nwstats$deg.mp[2,])
    mdeg.main.B.OW <- sum(nwstats$deg.mp.B[2,])*nwstats$mean.degm.OW/sum(nwstats$deg.mp[2,])
    mdeg.main.O.OW <- sum(nwstats$deg.mp.O[2,])*nwstats$mean.degm.OW/sum(nwstats$deg.mp[2,])
    mdeg.main.H.EW <- sum(nwstats$deg.mp.H[2,])*nwstats$mean.degm.EW/sum(nwstats$deg.mp[2,])
    mdeg.main.B.EW <- sum(nwstats$deg.mp.B[2,])*nwstats$mean.degm.EW/sum(nwstats$deg.mp[2,])
    mdeg.main.O.EW <- sum(nwstats$deg.mp.O[2,])*nwstats$mean.degm.EW/sum(nwstats$deg.mp[2,])
 
    dist.H.KC <- c(1 - mdeg.main.H.KC, mdeg.main.H.KC)
    dist.B.KC <- c(1 - mdeg.main.B.KC, mdeg.main.B.KC)
    dist.O.KC <- c(1 - mdeg.main.O.KC, mdeg.main.O.KC)
    dist.H.OW <- c(1 - mdeg.main.H.OW, mdeg.main.H.OW)
    dist.B.OW <- c(1 - mdeg.main.B.OW, mdeg.main.B.OW)
    dist.O.OW <- c(1 - mdeg.main.O.OW, mdeg.main.O.OW)
    dist.H.EW <- c(1 - mdeg.main.H.EW, mdeg.main.H.EW)
    dist.B.EW <- c(1 - mdeg.main.B.EW, mdeg.main.B.EW)
    dist.O.EW <- c(1 - mdeg.main.O.EW, mdeg.main.O.EW)
    
    race..wa <- get.vertex.attribute(nw, "race..wa")
    region <- get.vertex.attribute(nw, "region")
    vH.KC <- which(race..wa == "H" & region == "KC")
    vB.KC <- which(race..wa == "B" & region == "KC")
    vO.KC <- which(race..wa == "O" & region == "KC")
    vH.OW <- which(race..wa == "H" & region == "OW")
    vB.OW <- which(race..wa == "B" & region == "OW")
    vO.OW <- which(race..wa == "O" & region == "OW")
    vH.EW <- which(race..wa == "H" & region == "EW")
    vB.EW <- which(race..wa == "B" & region == "EW")
    vO.EW <- which(race..wa == "O" & region == "EW")
    nH.KC <- length(vH.KC)
    nB.KC <- length(vB.KC)
    nO.KC <- length(vO.KC)
    nH.OW <- length(vH.OW)
    nB.OW <- length(vB.OW)
    nO.OW <- length(vO.OW)
    nH.EW <- length(vH.EW)
    nB.EW <- length(vB.EW)
    nO.EW <- length(vO.EW)
    
    num.degrees <- length(dist.H.KC)
    
    deg.H.KC <- apportion_lr(nH.KC, 0:(num.degrees - 1), dist.H.KC, shuffled = TRUE)
    deg.B.KC <- apportion_lr(nB.KC, 0:(num.degrees - 1), dist.B.KC, shuffled = TRUE)
    deg.O.KC <- apportion_lr(nO.KC, 0:(num.degrees - 1), dist.O.KC, shuffled = TRUE)
    deg.H.OW <- apportion_lr(nH.OW, 0:(num.degrees - 1), dist.H.OW, shuffled = TRUE)
    deg.B.OW <- apportion_lr(nB.OW, 0:(num.degrees - 1), dist.B.OW, shuffled = TRUE)
    deg.O.OW <- apportion_lr(nO.OW, 0:(num.degrees - 1), dist.O.OW, shuffled = TRUE)
    deg.H.EW <- apportion_lr(nH.EW, 0:(num.degrees - 1), dist.H.EW, shuffled = TRUE)
    deg.B.EW <- apportion_lr(nB.EW, 0:(num.degrees - 1), dist.B.EW, shuffled = TRUE)
    deg.O.EW <- apportion_lr(nO.EW, 0:(num.degrees - 1), dist.O.EW, shuffled = TRUE)
    
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.H.KC, v = vH.KC)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.B.KC, v = vB.KC)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.O.KC, v = vO.KC)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.H.OW, v = vH.OW)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.B.OW, v = vB.OW)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.O.OW, v = vO.OW)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.H.EW, v = vH.EW)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.B.EW, v = vB.EW)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.O.EW, v = vO.EW)
  }
  
  if (deg.type == "pers") {
    
    attr.name <- "deg.pers"
    
    #Calculate expected proportion with pers degree 1 and 2+ for each race by region combination assuming independence
    deg.pers1.H.KC <- sum(nwstats$deg.mp.H[,2])*nwstats$deg.p.KC[2]/sum(nwstats$deg.mp[,2])
    deg.pers1.B.KC <- sum(nwstats$deg.mp.B[,2])*nwstats$deg.p.KC[2]/sum(nwstats$deg.mp[,2])
    deg.pers1.O.KC <- sum(nwstats$deg.mp.O[,2])*nwstats$deg.p.KC[2]/sum(nwstats$deg.mp[,2])
    deg.pers1.H.OW <- sum(nwstats$deg.mp.H[,2])*nwstats$deg.p.OW[2]/sum(nwstats$deg.mp[,2])
    deg.pers1.B.OW <- sum(nwstats$deg.mp.B[,2])*nwstats$deg.p.OW[2]/sum(nwstats$deg.mp[,2])
    deg.pers1.O.OW <- sum(nwstats$deg.mp.O[,2])*nwstats$deg.p.OW[2]/sum(nwstats$deg.mp[,2])
    deg.pers1.H.EW <- sum(nwstats$deg.mp.H[,2])*nwstats$deg.p.EW[2]/sum(nwstats$deg.mp[,2])
    deg.pers1.B.EW <- sum(nwstats$deg.mp.B[,2])*nwstats$deg.p.EW[2]/sum(nwstats$deg.mp[,2])
    deg.pers1.O.EW <- sum(nwstats$deg.mp.O[,2])*nwstats$deg.p.EW[2]/sum(nwstats$deg.mp[,2])
    
    deg.pers2.H.KC <- sum(nwstats$deg.mp.H[,3])*nwstats$deg.p.KC[3]/sum(nwstats$deg.mp[,3])
    deg.pers2.B.KC <- sum(nwstats$deg.mp.B[,3])*nwstats$deg.p.KC[3]/sum(nwstats$deg.mp[,3])
    deg.pers2.O.KC <- sum(nwstats$deg.mp.O[,3])*nwstats$deg.p.KC[3]/sum(nwstats$deg.mp[,3])
    deg.pers2.H.OW <- sum(nwstats$deg.mp.H[,3])*nwstats$deg.p.OW[3]/sum(nwstats$deg.mp[,3])
    deg.pers2.B.OW <- sum(nwstats$deg.mp.B[,3])*nwstats$deg.p.OW[3]/sum(nwstats$deg.mp[,3])
    deg.pers2.O.OW <- sum(nwstats$deg.mp.O[,3])*nwstats$deg.p.OW[3]/sum(nwstats$deg.mp[,3])
    deg.pers2.H.EW <- sum(nwstats$deg.mp.H[,3])*nwstats$deg.p.EW[3]/sum(nwstats$deg.mp[,3])
    deg.pers2.B.EW <- sum(nwstats$deg.mp.B[,3])*nwstats$deg.p.EW[3]/sum(nwstats$deg.mp[,3])
    deg.pers2.O.EW <- sum(nwstats$deg.mp.O[,3])*nwstats$deg.p.EW[3]/sum(nwstats$deg.mp[,3])
    
    dist.H.KC <- c((1 - sum(deg.pers1.H.KC, deg.pers2.H.KC)), deg.pers1.H.KC, deg.pers2.H.KC)
    dist.B.KC <- c((1 - sum(deg.pers1.B.KC, deg.pers2.B.KC)), deg.pers1.B.KC, deg.pers2.B.KC)
    dist.O.KC <- c((1 - sum(deg.pers1.O.KC, deg.pers2.O.KC)), deg.pers1.O.KC, deg.pers2.O.KC)
    dist.H.OW <- c((1 - sum(deg.pers1.H.OW, deg.pers2.H.OW)), deg.pers1.H.OW, deg.pers2.H.OW)
    dist.B.OW <- c((1 - sum(deg.pers1.B.OW, deg.pers2.B.OW)), deg.pers1.B.OW, deg.pers2.B.OW)
    dist.O.OW <- c((1 - sum(deg.pers1.O.OW, deg.pers2.O.OW)), deg.pers1.O.OW, deg.pers2.O.OW)
    dist.H.EW <- c((1 - sum(deg.pers1.H.EW, deg.pers2.H.EW)), deg.pers1.H.EW, deg.pers2.H.EW)
    dist.B.EW <- c((1 - sum(deg.pers1.B.EW, deg.pers2.B.EW)), deg.pers1.B.EW, deg.pers2.B.EW)
    dist.O.EW <- c((1 - sum(deg.pers1.O.EW, deg.pers2.O.EW)), deg.pers1.O.EW, deg.pers2.O.EW)
    
    race..wa <- get.vertex.attribute(nw, "race..wa")
    region <- get.vertex.attribute(nw, "region")
    vH.KC <- which(race..wa == "H" & region == "KC")
    vB.KC <- which(race..wa == "B" & region == "KC")
    vO.KC <- which(race..wa == "O" & region == "KC")
    vH.OW <- which(race..wa == "H" & region == "OW")
    vB.OW <- which(race..wa == "B" & region == "OW")
    vO.OW <- which(race..wa == "O" & region == "OW")
    vH.EW <- which(race..wa == "H" & region == "EW")
    vB.EW <- which(race..wa == "B" & region == "EW")
    vO.EW <- which(race..wa == "O" & region == "EW")
    nH.KC <- length(vH.KC)
    nB.KC <- length(vB.KC)
    nO.KC <- length(vO.KC)
    nH.OW <- length(vH.OW)
    nB.OW <- length(vB.OW)
    nO.OW <- length(vO.OW)
    nH.EW <- length(vH.EW)
    nB.EW <- length(vB.EW)
    nO.EW <- length(vO.EW)
    
    num.degrees <- length(dist.H.KC)
    
    deg.H.KC <- apportion_lr(nH.KC, 0:(num.degrees - 1), dist.H.KC, shuffled = TRUE)
    deg.B.KC <- apportion_lr(nB.KC, 0:(num.degrees - 1), dist.B.KC, shuffled = TRUE)
    deg.O.KC <- apportion_lr(nO.KC, 0:(num.degrees - 1), dist.O.KC, shuffled = TRUE)
    deg.H.OW <- apportion_lr(nH.OW, 0:(num.degrees - 1), dist.H.OW, shuffled = TRUE)
    deg.B.OW <- apportion_lr(nB.OW, 0:(num.degrees - 1), dist.B.OW, shuffled = TRUE)
    deg.O.OW <- apportion_lr(nO.OW, 0:(num.degrees - 1), dist.O.OW, shuffled = TRUE)
    deg.H.EW <- apportion_lr(nH.EW, 0:(num.degrees - 1), dist.H.EW, shuffled = TRUE)
    deg.B.EW <- apportion_lr(nB.EW, 0:(num.degrees - 1), dist.B.EW, shuffled = TRUE)
    deg.O.EW <- apportion_lr(nO.EW, 0:(num.degrees - 1), dist.O.EW, shuffled = TRUE)
    
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.H.KC, v = vH.KC)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.B.KC, v = vB.KC)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.O.KC, v = vO.KC)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.H.OW, v = vH.OW)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.B.OW, v = vB.OW)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.O.OW, v = vO.OW)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.H.EW, v = vH.EW)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.B.EW, v = vB.EW)
    nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.O.EW, v = vO.EW)
    
  }


  return(nw)
}


# Het -----------------------------------------------------------------

#' @title Calculate Network Statistics
#'
#' @description This function calculates the target statistics for the formation
#'              and dissolution models estimated in \code{netest}.
#'
#' @param n Population size.
#' @param meandeg Mean degree.
#' @param prop.male Percent of the population that is male.
#' @param start.prev Starting HIV prevalence in the population.
#' @param part.dur Mean duration of partnerships.
#' @param time.unit Time unit used, relative to days.
#'
#' @keywords het
#' @export
#'
make_nw_het <- function(n = 10000,
                        meandeg = 0.8,
                        prop.male = 0.5,
                        start.prev = 0.05,
                        part.dur = 1000,
                        time.unit = 7) {

  nMale <- round(n * prop.male)

  male <- rep(0, n)
  male[sample(1:n, nMale)] <- 1

  # Set vertex attributes
  nw <- network.initialize(n = n, directed = FALSE)
  nw <- set.vertex.attribute(nw, attrname = "male", value = male)

  # Formation Model
  formation <- ~edges + offset(nodematch("male"))

  # Target stats
  edges.ts <- meandeg * (n/2)

  stats <- edges.ts

  # Dissolution model
  dissolution <- ~offset(edges)
  dur <- part.dur/time.unit
  d.rate <- time.unit * (((1 - start.prev) * 1/(55 - 18)/365) + (start.prev * 1/12/365))
  coef.diss <- dissolution_coefs(dissolution, duration = dur, d.rate = d.rate)

  out <- list()
  out$nw <- nw
  out$time.unit <- time.unit
  out$formation <- formation
  out$stats <- stats
  out$coef.diss <- coef.diss

  return(out)
}
