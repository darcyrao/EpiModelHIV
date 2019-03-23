
# MSM -----------------------------------------------------------------

#' @title Epidemic Model Parameters for the WHAMP model
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with \code{\link{netsim}} for EpiModelHIV
#'
#' @param nwstats Target statistics for the network model. An object of class
#'        \code{nwstats} output from \code{\link{calc_nwstats_msm}}.
#' @param change.size If \code{TRUE}, change the size of the network for epidemic simulations
#' @param new.size Size to which to scale the population for epidemic simulations
#' @param racedist Population distribution by race/ethnicity and region (King County Hispanic,
#'        black, other; other western WA Hispanic, black, and other, eastern WA Hispanic,
#'        black, and other)
#' @param agedist Target population age structure (proportion in each group 18-24, 25-29... 55-59)
#' @param nw_track Sets whether to run the script in the nwfeatures_msm_whamp function to
#'        track features of the sexual network (degree distribution, mixing, partnership age)
#'        over time, with options \code{TRUE} or \code{FALSE}
#' @param test.int Average interest interval in days for men who test regularly
#' @param testing.pattern Method for HIV testing, with options \code{"memoryless"}
#'        for constant hazard without regard to time since previous test, or
#'        \code{"interval"} deterministic fixed intervals.
#' @param test.window.int Length of the HIV test window period in days.
#' @param tt.traj.KC.B.prob Proportion of black MSM in King County who enter one of four
#'        testing/treatment trajectories: non-screener treated with partial viral suppression,
#'        non-screener treated with full viral suppression, regular screener treated with partial
#'        suppression, and regular screener treated with full suppression
#' @param tt.traj.KC.H.prob Proportion of Hispanic MSM in King County who enter into the four
#'        testing/treatment trajectories, as defined above.
#' @param tt.traj.KC.O.prob Proportion of other MSM in King County who enter into the four
#'        testing/treatment trajectories, as defined above.
#' @param tt.traj.OW.B.prob Proportion of black MSM in Western WA who enter one of four
#'        testing/treatment trajectories, as defined above
#' @param tt.traj.OW.H.prob Proportion of Hispanic MSM in Western WA who enter into the four
#'        testing/treatment trajectories, as defined above.
#' @param tt.traj.OW.O.prob Proportion of other MSM in Western WA who enter into the four
#'        testing/treatment trajectories, as defined above.
#' @param tt.traj.EW.B.prob Proportion of black MSM in Eastern WA who enter one of four
#'        testing/treatment trajectories, as defined above
#' @param tt.traj.EW.H.prob Proportion of Hispanic MSM in Eastern WA who enter into the four
#'        testing/treatment trajectories, as defined above.
#' @param tt.traj.EW.O.prob Proportion of other MSM in Eastern WA who enter into the four
#'        testing/treatment trajectories, as defined above.
#' @param tx.init.int.KC.B Days from diagnosis to treatment initiation for Black MSM in King County
#' @param tx.init.int.KC.H Days from diagnosis to treatment initiation for Hispanic MSM in King County
#' @param tx.init.int.KC.O Days from diagnosis to treatment initiation for Other MSM in King County
#' @param tx.init.int.OW.B Days from diagnosis to treatment initiation for Black MSM in other western WA
#' @param tx.init.int.OW.H Days from diagnosis to treatment initiation for Hispanic MSM in other western WA
#' @param tx.init.int.OW.O Days from diagnosis to treatment initiation for Other MSM in other western WA
#' @param tx.init.int.EW.B Days from diagnosis to treatment initiation for Black MSM in eastern WA
#' @param tx.init.int.EW.H Days from diagnosis to treatment initiation for Hispanic MSM in eastern WA
#' @param tx.init.int.EW.O Days from diagnosis to treatment initiation for Other MSM in eastern WA
#' @param tx.halt.full Probability per time step that a full suppressor who is
#'        currently on treatment will halt treatment.
#' @param tx.halt.part.rr Relative risk of halting treatment per time step for a partial suppressor
#'        currently on treatment, as compared to a full suppressor
#' @param tx.reinit.full.KC.B Probability per time step that black MSM in KC who are full suppressors and
#'        not currently on treatment but who have been in the past will re-initiate
#' @param tx.reinit.full.KC.H Probability per time step that Hispanic MSM in KC who are full suppressors and
#'        not currently on treatment but who have been in the past will re-initiate
#' @param tx.reinit.full.KC.O Probability per time step that other MSM in KC who are full suppressors and
#'        not currently on treatment but who have been in the past will re-initiate
#' @param tx.reinit.full.OW.B Probability per time step that black MSM in OW who are full suppressors and
#'        not currently on treatment but who have been in the past will re-initiate
#' @param tx.reinit.full.OW.H Probability per time step that Hispanic MSM in OW who are full suppressors and
#'        not currently on treatment but who have been in the past will re-initiate
#' @param tx.reinit.full.OW.O Probability per time step that other MSM in OW who are full suppressors and
#'        not currently on treatment but who have been in the past will re-initiate
#' @param tx.reinit.full.EW.B Probability per time step that black MSM in EW who are full suppressors and
#'        not currently on treatment but who have been in the past will re-initiate
#' @param tx.reinit.full.EW.H Probability per time step that Hispanic MSM in EW who are full suppressors and
#'        not currently on treatment but who have been in the past will re-initiate
#' @param tx.reinit.full.EW.O Probability per time step that other MSM in EW who are full suppressors and
#'        not currently on treatment but who have been in the past will re-initiate
#' @param tx.reinit.part.rr Relative risk of reinitiating treatment per time step for a partial suppressor
#'        currently on treatment, as compared to a full suppressor
#' @param max.time.off.tx.int Number of days off treatment for a before onset of AIDS, 
#'        including time before diagnosis.
#' @param sympt.onset.int Number of days from infection to symptom onset in the absence of ART
#' @param vl.acute.rise.int Number of days to peak viremia during acute infection.
#' @param vl.acute.peak Peak viral load (in log10 units) at the height of acute infection.
#' @param vl.acute.fall.int Number of days from peak viremia to set-point
#'        viral load during the acute infection period.
#' @param vl.set.point Set point viral load (in log10 units).
#' @param vl.aids.onset.int Number of days to AIDS for a treatment-naive
#'        patient.
#' @param vl.aids.int Time from AIDS onset to peak AIDS VL in untreated infection in days.
#' @param vl.aids.peak Peak VL in AIDS.
#' @param vl.full.supp Log10 viral load at full suppression on ART.
#' @param vl.part.supp Log10 viral load at partial suppression on ART.
#' @param full.supp.down.slope For full suppressors, number of log10 units that
#'        viral load falls per time step from treatment initiation or re-initiation
#'        until the level in \code{vl.full.supp}.
#' @param full.supp.up.slope For full suppressors, number of log10 units that
#'        viral load rises per time step from treatment halting until expected
#'        value.
#' @param part.supp.down.slope For partial suppressors, number of log10 units
#'        that viral load falls per time step from treatment initiation or
#'        re-initiation until the level in \code{vl.part.supp}.
#' @param part.supp.up.slope For partial suppressors, number of log10 units that
#'        viral load rises per time step from treatment halting until expected value.
#' @param growth.rate Rate of population growth per day.        
#' @param birth.age Age (in years) of new arrivals.
#' @param exit.age Age (in years) at which individuals age out of the network
#' @param asmr.rr.pos Relative increase in age-specific mortality for HIV-positive men
#'        (for age groups 18-44, 45-54, 55-59)
#' @param URAI.prob Probability of transmission for a man having unprotected
#'        receptive anal intercourse with an infected man at set point viral
#'        load.
#' @param UIAI.prob Probability of transmission for an uncircumcised man having
#'        unprotected insertive anal intercourse with an infected man at set
#'        point viral load.
#' @param acute.rr Relative risk of infection (compared to that predicted by
#'        elevated viral load) when positive partner is in the acute stage.
#' @param circ.rr Relative risk of infection from insertive anal sex when the
#'        negative insertive partner is circumcised.
#' @param condom.rr Relative risk of infection from anal sex when a condom is
#'        used.
#' @param tprob.scalar General relative scaler for all per-act probability of 
#'        transmission for model calibration.
#' @param disc.outset.main.prob Probability that an HIV-infected MSM will
#'        disclose his status at the start of a main partnership.
#' @param disc.at.diag.main.prob Probability that an MSM already in a main
#'        partnership will disclose at the time of diagnosis.
#' @param disc.post.diag.main.prob Probability that an HIV-infected MSM
#'        in a main partnership will disclose his status, assuming he didn't
#'        at the start of the partnership or at diagnosis.
#' @param disc.outset.pers.prob Probability that an HIV-infected MSM will
#'        disclose his status at the start of a casual partnership.
#' @param disc.at.diag.pers.prob Probability that an MSM already in a
#'        casual partnership will disclose at the time of diagnosis.
#' @param disc.post.diag.pers.prob Probability that an HIV-infected MSM
#'        in a casual partnership will disclose his status, assuming he
#'        didn't at the start of the partnership or at diagnosis.
#' @param disc.inst.prob Probability that an HIV-infected black MSM will
#'        disclose his status to a one-off partner.
#' @param circ.B.prob Probablity that a black new arrival in the population
#'        will be circumcised.
#' @param circ.H.prob Probablity that a Hispanic new arrival in the population
#'        will be circumcised.
#' @param circ.O.prob Probablity that an other race/ethnicity new arrival in 
#'        the population will be circumcised.
#' @param ccr5.B.prob Vector of length two of frequencies of the Delta 32
#'        mutation (homozygous and heterozygous, respectively) in the CCR5 gene
#'        among black MSM.
#' @param ccr5.H.prob Vector of length two of frequencies of the Delta 32
#'        mutation (homozygous and heterozygous, respectively) in the CCR5 gene
#'        among hispanic MSM.
#' @param ccr5.O.prob Vector of length two of frequencies of the Delta 32
#'        mutation (homozygous and heterozygous, respectively) in the CCR5 gene
#'        among other race/ethnicity MSM.
#' @param ccr5.heteroz.rr Relative risk of infection for men who are heterozygous
#'        in the CCR5 mutation.
#' @param base.ai.main.rate Expected coital frequency in main partnerships
#'        (acts per day).
#' @param base.ai.pers.rate Expected coital frequency in persistent partnerships
#'        (acts per day).
#' @param ai.scale General relative scaler for all act rates for model
#'        calibration.
#' @param cond.main.YY.prob Probability of condom use in a Young-Young (18-34) main
#'        partnership.
#' @param cond.main.other.prob Probability of condom use in main partnerships with
#'        other age combinations.
#' @param cond.pers.always.prob Fraction of men in persistent partnerships who always
#'        use condoms in those partnerships.
#' @param cond.pers.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in persistent partnerships.
#' @param cond.inst.always.prob Fraction of men in instant partnerships who always
#'        use condoms in those partnerships.
#' @param cond.inst.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in a instantaneous partnerships.
#' @param cond.always.prob.corr Correlation coefficient for probability of always
#'        using condoms in both persistent and instantaneous partnerships
#' @param cond.rr Condom probability scaler for model calibration purposes.
#' @param cond.diag.main.beta Beta multiplier for the log odds of using a
#'        condom in a main partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.main.beta Beta multiplier for the log odds of using a
#'        condom in a main partnership if the HIV-infected man has disclosed.
#' @param cond.diag.pers.beta Beta multiplier for the log odds of using a
#'        condom in a persistent partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.pers.beta Beta multiplier for the log odds of using a
#'        condom in a persistent partnership if the HIV-infected man has disclosed
#'        his status.
#' @param cond.diag.inst.beta Beta multiplier for the log odds of using a
#'        condom in an instantaneous partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.inst.beta Beta multiplier for the log odds of using a
#'        condom in an instantaneous partnership if the HIV-infected man has disclosed
#'        his status.
#' @param vv.iev.prob Probability that in a partnership of two versatile men,
#'        they will engage in intra-event versatility ("flipping") given that 
#'        they're having AI.
#'
#' @param prep.start Time step at which the PrEP intervention should start.
#' @param prep.class.prob The probability of adherence class in low adherence,
#'        medium adherence, or high adherence groups.
#' @param prep.class.hr The hazard ratio for infection per act associated with each
#'        level of adherence.
#' @param prep.coverage.init.region Vector of length 2 for the proportion of eligible men 
#'        (eligibility defined by CAI and discordant partnership status in line with 
#'        WA State guidelines) in King County and other counties in WA who are allowed 
#'        to start PrEP once they become eligible at the time step in which PrEP is initiated
#'        (as specified in \code{prep.start}).
#' @param prep.cov.method The method for calculating PrEP coverage, with options
#'        of \code{"curr"} to base the numerator on the number of people currently
#'        on PrEP and \code{"ever"} to base it on the number of people ever on
#'        PrEP.
#' @param prep.scaleup.rate Vector of length 2 for the rate at which PrEP is  
#'        scaled up per day from \code{prep.start} for each region (KC vs. other)
#' @param prep.cov.max.region Vector of length 2 for the maximum attainable PrEP coverage
#'        for each region (King County vs. other counties). Upon reaching these values, 
#'        PrEP uptake will stabilize.
#' @param prep.init.rate The rate at which persons initiate PrEP conditional on
#'        their eligibility, with 1 equal to instant start.
#' @param prep.tst.int Testing interval for those who are actively on PrEP. This
#'        overrides the mean testing interval parameters.
#' @param prep.risk.int Time window for assessment of risk eligibility for PrEP
#'        in days.
#' @param prep.risk.reassess.method Method for determining risk-based discontinuation
#'        of PrEP, with \code{"none"} for no risk-based discontinuation, \code{"inst"} 
#'        for reassessment every time step, and \code{"year"} for reassessment at yearly
#'        HIV screening visits.
#' @param prep.discont <- Proportion of PrEP users who will discontinue while still at risk.
#'        This will be used to assign an attribute to indicate men who will discontinue
#'        while they are still candidates for PrEP. The remainder will discontinue only
#'        if their risk changes such that they are no longer eligilbe for the intervention.
#' @param prep.discont.prob Probability per time step of halting PrEP for men who 
#'        discontinue.
#'
#' @param rcomp.prob Level of risk compensation from 0 to 1, where 0 is no risk
#'        compensation, 0.5 is a 50% reduction in the probability of condom use
#'        per act, and 1 is a complete cessation of condom use following PrEP
#'        initiation.
#' @param rcomp.adh.groups PrEP adherence groups for whom risk compensation
#'        occurs, as a vector with values 1, 2, 3 corresponding to
#'        low adherence, medium adherence, and high adherence to PrEP.
#' @param rcomp.main.only Logical, if risk compensation is limited to main
#'        partnerships only, versus all partnerships.
#' @param rcomp.discl.only Logical, if risk compensation is limited known-discordant
#'        partnerships only, versus all partnerships.
#' @param rcomp.discont Logical, if condom use stays at the level at which PrEP was 
#'        discontinued versus returning to pre-PrEP levels
#'
#' @param rgc.tprob Probability of rectal gonorrhea infection per act.
#' @param ugc.tprob Probability of urethral gonorrhea infection per act.
#' @param rct.tprob Probability of rectal chlamydia infection per act.
#' @param uct.tprob Probability of urethral chlamydia infection per act.
#' @param rgc.sympt.prob Probability of symptoms given infection with rectal
#'        gonorrhea.
#' @param ugc.sympt.prob Probability of symptoms given infection with urethral
#'        gonorrhea.
#' @param rct.sympt.prob Probability of symptoms given infection with rectal
#'        chlamydia.
#' @param uct.sympt.prob Probability of symptoms given infection with urethral
#'        chlamydia.
#' @param rgc.asympt.int Average duration in days of asymptomatic rectal gonorrhea.
#' @param ugc.asympt.int Average duration in days of asymptomatic urethral gonorrhea.
#' @param gc.tx.int Average duration in days of treated gonorrhea (both sites).
#' @param gc.ntx.int Average duration in days of untreated, symptomatic gonorrhea (both sites).
#'        If \code{NA}, uses site-specific durations for asymptomatic infections.
#' @param rct.asympt.int Average in days duration of asymptomatic rectal chlamydia.
#' @param uct.asympt.int Average in days duration of asymptomatic urethral chlamydia.
#' @param ct.tx.int Average in days duration of treated chlamydia (both sites).
#' @param ct.ntx.int Average in days duration of untreated, symptomatic chlamydia (both sites).
#'        If \code{NA}, uses site-specific durations for asymptomatic infections.
#' @param gc.prob.cease Probability of ceasing sexual activity during symptomatic
#'        infection with gonorrhea.
#' @param ct.prob.cease Probability of ceasing sexual activity during symptomatic
#'        infection with chlamydia.
#' @param gc.sympt.prob.tx Probability of treatment for symptomatic gonorrhea.
#' @param ct.sympt.prob.tx Probability of treatment for symptomatic chlamydia.
#' @param gc.asympt.prob.tx Probability of treatment for asymptomatic gonorrhea.
#' @param ct.asympt.prob.tx Probability of treatment for asymptomatic chlamydia.
#' @param prep.sti.screen.int Interval in days between STI screening at PrEP visits.
#' @param prep.sti.prob.tx Probability of treatment given positive screening during
#'        PrEP visit.
#' @param prep.continue.stand.tx Logical, if \code{TRUE} will continue standard
#'        STI treatment of symptomatic cases even after PrEP initiation.
#' @param sti.cond.rr Relative risk of STI infection (in either direction) given
#'        a condom used by the insertive partner.
#' @param hiv.rgc.rr Relative risk of HIV infection given current rectal gonorrhea.
#' @param hiv.ugc.rr Relative risk of HIV infection given current urethral gonorrhea.
#' @param hiv.rct.rr Relative risk of HIV infection given current rectal chlamydia.
#' @param hiv.uct.rr Relative risk of HIV infection given current urethral chlamydia.
#' @param hiv.dual.rr Additive proportional risk, from 0 to 1, for HIV infection
#'        given dual infection with both gonorrhea and chlamydia.
#'
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{param_msm}, which can be passed to
#' EpiModel function \code{netsim}.
#'
#' @keywords msm
#'
#' @export
#'
param_msm_whamp <- function(nwstats,
                          
                      change.size = FALSE,
                      new.size = 90000,
                      racedist = sumto1(c(0.0549, 0.0421, 0.4739, 0.0309, 0.0166, 0.2807, 0.0222, 0.0021, 0.0767)),
                      agedist = sumto1(c(0.1594, 0.1319, 0.1292, 0.1173, 0.1183, 0.1148, 0.1071, 0.122)),
                      nw_track = TRUE,
                      
                      test.int = 436,
                      testing.pattern = "interval",
                      test.window.int = 21,

                      tt.traj.KC.B.prob = c(0.080*(1-0.839), 0.080*(0.839), (1-0.080)*(1-0.839), (1-0.080)*(0.839)),
                      tt.traj.KC.H.prob = c(0.064*(1-0.902), 0.064*(0.902), (1-0.064)*(1-0.902), (1-0.064)*(0.902)),
                      tt.traj.KC.O.prob = c(0.104*(1-0.920), 0.104*(0.920), (1-0.104)*(1-0.920), (1-0.104)*(0.920)),
                      tt.traj.OW.B.prob = c(0.120*(1-0.811), 0.120*(0.811), (1-0.120)*(1-0.811), (1-0.120)*(0.811)),
                      tt.traj.OW.H.prob = c(0.102*(1-0.885), 0.102*(0.885), (1-0.102)*(1-0.885), (1-0.102)*(0.885)),
                      tt.traj.OW.O.prob = c(0.148*(1-0.905), 0.148*(0.905), (1-0.148)*(1-0.905), (1-0.148)*(0.905)),
                      tt.traj.EW.B.prob = c(0.161*(1-0.759), 0.161*(0.759), (1-0.161)*(1-0.759), (1-0.161)*(0.759)),
                      tt.traj.EW.H.prob = c(0.14*(1-0.847), 0.14*(0.847), (1-0.14)*(1-0.847), (1-0.14)*(0.847)),
                      tt.traj.EW.O.prob = c(0.192*(1-0.874), 0.192*(0.874), (1-0.192)*(1-0.874), (1-0.192)*(0.874)),
                      
                      tx.init.int.KC.B = 46,
                      tx.init.int.KC.H = 43,
                      tx.init.int.KC.O = 47,
                      tx.init.int.OW.B = 56,
                      tx.init.int.OW.H = 53,
                      tx.init.int.OW.O = 57,
                      tx.init.int.EW.B = 53,
                      tx.init.int.EW.H = 50,
                      tx.init.int.EW.O = 53,
                      
                      tx.halt.full = 0.002226782,
                      tx.halt.part.rr = 2,
                      tx.reinit.full.KC.B = 0.0244,
                      tx.reinit.full.KC.H = 0.0239,
                      tx.reinit.full.KC.O = 0.0275,
                      tx.reinit.full.OW.B = 0.0200,
                      tx.reinit.full.OW.H = 0.0195,
                      tx.reinit.full.OW.O = 0.0223,
                      tx.reinit.full.EW.B = 0.0219,
                      tx.reinit.full.EW.H = 0.0209,
                      tx.reinit.full.EW.O = 0.0237,
                      tx.reinit.part.rr = 0.5,

                      max.time.off.tx.int = 520 * 7,
                      sympt.onset.int = 2737.5,
                      vl.acute.rise.int = 45,
                      vl.acute.peak = 6.886,
                      vl.acute.fall.int = 45,
                      vl.set.point = 4.5,
                      vl.aids.onset.int = 520 * 7,
                      vl.aids.int = 52 * 2 * 7,
                      vl.aids.peak = 7,
                      vl.full.supp = 1.5,
                      vl.part.supp = 3.5,
                      full.supp.down.slope = 0.75,
                      full.supp.up.slope = 0.75,
                      part.supp.down.slope = 0.25,
                      part.supp.up.slope = 0.75,

                      growth.rate = 1.0001441 / 7,
                      birth.age = 18,
                      exit.age = 60,
                      asmr.rr.pos = c(3.791, 2.974, 1.984),

                      URAI.prob = 0.0082 * 1.09,
                      UIAI.prob = 0.0031 * 1.09,
                      acute.rr = 6,
                      circ.rr = 0.4,
                      condom.rr = 0.295,
                      tprob.scalar = 1,

                      disc.outset.main.prob = 1,
                      disc.at.diag.main.prob = 1,
                      disc.post.diag.main.prob = 0,
                      disc.outset.pers.prob = 0.5671,
                      disc.at.diag.pers.prob = 1,
                      disc.post.diag.pers.prob = 0,
                      disc.inst.prob = 0.4918,

                      circ.B.prob = 0.6449,
                      circ.H.prob = 0.4897,
                      circ.O.prob = 0.86,

                      ccr5.B.prob = c(0, 0.034),
                      ccr5.H.prob = c(0.003, 0.050),
                      ccr5.O.prob = c(0.017, 0.164), 
                      ccr5.heteroz.rr = 0.3,

                      base.ai.main.rate = 0.1864,
                      base.ai.pers.rate = 0.118,
                      ai.scale = 1, # set to 1 for now

                      cond.main.YY.prob = 0.3391,
                      cond.main.other.prob = 0.1399,
                      cond.pers.always.prob = 0.0924,
                      cond.pers.prob = 0.4001,
                      cond.inst.always.prob = 0.1982,
                      cond.inst.prob = 0.4138,
                      cond.always.prob.corr = 0.6009,
                      cond.rr = 1,
                      cond.diag.main.beta = -0.67,
                      cond.discl.main.beta = -0.85,
                      cond.diag.pers.beta = -0.67,
                      cond.discl.pers.beta = -0.85,
                      cond.diag.inst.beta = -0.67,
                      cond.discl.inst.beta = -0.85,

                      vv.iev.prob = 0.42,

                      prep.start = Inf, # Set to Inf for no PrEP
                      prep.class.prob = c(0.089, 0.127, 0.784),
                      prep.class.hr = c(0.69, 0.19, 0.05),
                      prep.coverage.init.region = c(0.4392,	0.2655),
                      prep.cov.method = "curr",
                      prep.scaleup.rate = c(0, 0), # Set to 0 for stable PrEP use
                      prep.cov.max.region = c(0.7497,	0.6232),
                      prep.init.rate = 1,
                      prep.tst.int = 90,
                      prep.risk.int = 365,
                      prep.risk.reassess.method = "inst", #year
                      prep.discont = 0.3,
                      prep.discont.prob = 0.017724,

                      rcomp.prob = 0,
                      rcomp.adh.groups = 1:3,
                      rcomp.main.only = FALSE,
                      rcomp.discl.only = FALSE,
                      rcomp.discont = FALSE,

                      rgc.tprob = 0.357698,
                      ugc.tprob = 0.248095,
                      rct.tprob = 0.321597,
                      uct.tprob = 0.212965,

                      rgc.sympt.prob = 0.076975,
                      ugc.sympt.prob = 0.824368,
                      rct.sympt.prob = 0.103517,
                      uct.sympt.prob = 0.885045,

                      rgc.asympt.int = 35.11851 * 7,
                      ugc.asympt.int = 35.11851 * 7,
                      gc.tx.int = 2 * 7,
                      gc.ntx.int = NA,

                      rct.asympt.int = 44.24538 * 7,
                      uct.asympt.int = 44.24538 * 7,
                      ct.tx.int = 2 * 7,
                      ct.ntx.int = NA,

                      gc.prob.cease = 0,
                      ct.prob.cease = 0,

                      gc.sympt.prob.tx = 0.90,
                      ct.sympt.prob.tx = 0.85,
                      gc.asympt.prob.tx = 0,
                      ct.asympt.prob.tx = 0,

                      prep.sti.screen.int = 182,
                      prep.sti.prob.tx = 1,
                      prep.continue.stand.tx = TRUE,

                      sti.cond.rr = 0.3,

                      hiv.rgc.rr = 2.780673,
                      hiv.ugc.rr = 1.732363,
                      hiv.rct.rr = 2.780673,
                      hiv.uct.rr = 1.732363,
                      hiv.dual.rr = 0.2,
                      ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  if (!(testing.pattern %in% c("memoryless", "interval"))) {
    stop("testing.pattern must be \"memoryless\" or \"interval\" ",
          call. = FALSE)
  }

  p$time.unit <- nwstats$time.unit

  intvars <- grep(names(p), pattern = ".int", fixed = TRUE)
  p[intvars] <- lapply(p[intvars], FUN = function(x) round(x / p$time.unit))

  ratevars <- grep(names(p), pattern = ".rate", fixed = TRUE)
  p[ratevars] <- lapply(p[ratevars], FUN = function(x) x * p$time.unit)

  p$role.prob <- nwstats$role.prob

  p$role.trans.matrix <- matrix(c(1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 1),
                                nrow = 3)

  p$riskh.start <- max(1, prep.start - p$prep.risk.int - 1)
  
  p$modes <- 1

  p$asmr.H..wa <- nwstats$asmr.H..wa
  p$asmr.B..wa <- nwstats$asmr.B..wa
  p$asmr.O..wa <- nwstats$asmr.O..wa
  
  
  p$nwstats <- NULL

  class(p) <- "param.net"
  return(p)
}


#' @title Epidemic Model Initial Conditions for the WHAMP model
#'
#' @description Sets the initial conditions for a stochastic epidemic models
#'              simulated with \code{\link{netsim}}.
#'
#' @param nwstats Target statistics for the network model. An object of class
#'        \code{nwstats} output from \code{\link{calc_nwstats_msm}}.
#' @param prev.H..wa Initial disease prevalence among Hispanic MSM.
#' @param prev.B..wa Initial disease prevalence among non-Hispanic black MSM.
#' @param prev.O..wa Initial disease prevalence among non-Hispanic other MSM.
#' @param prev.ugc Initial prevalence of urethral gonorrhea.
#' @param prev.rgc Initial prevalence of rectal gonorrhea.
#' @param prev.uct Initial prevalence of urethral chlamydia.
#' @param prev.rct Initial prevalence of rectal chlamydia.
#' @param ... Additional arguments passed to function.
#'
#' @return
#' A list object of class \code{init_msm}, which can be passed to EpiModel
#' function \code{\link{netsim}}.
#'
#' @keywords msm
#'
#' @export
init_msm_whamp <- function(nwstats,
                     prev.H..wa = 0.09, #-- Surveillance data on prevalent diagnoses suggest 7.7% prevalence for Other, 17.3% for Black, 12.0% for Hispanic. Note this underestimates true prevalence bc it only counts diagnosed cases. We set starting prevalence to slightly less than this for each group to reach equilibrium faster
                     prev.B..wa = 0.13,
                     prev.O..wa = 0.06,
                     prev.ugc = 0, #--This model will not represent STIs
                     prev.rgc = 0,
                     prev.uct = 0,
                     prev.rct = 0,
                     ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  p$num.H..wa <- nwstats$num.H.KC + nwstats$num.H.OW + nwstats$num.H.EW
  p$num.B..wa <- nwstats$num.B.KC + nwstats$num.B.OW + nwstats$num.B.EW
  p$num.O..wa <- nwstats$num.O.KC + nwstats$num.O.OW + nwstats$num.O.EW
  
  p$num.KC <- nwstats$num.H.KC + nwstats$num.B.KC + nwstats$num.O.KC
  p$num.OW <- nwstats$num.H.OW + nwstats$num.B.OW + nwstats$num.O.OW
  p$num.EW <- nwstats$num.H.EW + nwstats$num.B.EW + nwstats$num.O.EW
  
  p$ages <- nwstats$ages

  p$init.prev.age.slope.H..wa <- prev.H..wa / ((max(p$ages) - min(p$ages) + 1)/2) # Divide by half the age range to set the slope for prevalence by age
  p$init.prev.age.slope.B..wa <- prev.B..wa / ((max(p$ages) - min(p$ages) + 1)/2)
  p$init.prev.age.slope.O..wa <- prev.O..wa / ((max(p$ages) - min(p$ages) + 1)/2)
  
  p$nwstats <- NULL

  class(p) <- "init.net"
  return(p)
}


#' @title Epidemic Model Control Settings for the WHAMP model
#'
#' @description Sets the controls for stochastic network models simulated with
#'              \code{\link{netsim}}.
#'
#' @param simno Unique ID for the simulation run, used for file naming purposes
#'        if used in conjunction with the \code{EpiModelHPC} package.
#' @param nsims Number of simulations.
#' @param ncores Number of cores per run, if parallelization is used within the
#'        \code{EpiModelHPC} package.
#' @param nsteps Number of time steps per simulation.
#' @param start Starting time step for simulation, with default to 1 to run new
#'        simulation. This may also be set to 1 greater than the final time
#'        step of a previous simulation to resume the simulation with different
#'        parameters.
#' @param initialize.FUN Module function to use for initialization of the epidemic
#'        model.
#' @param aging.FUN Module function for aging.
#' @param deaths.FUN Module function for general and disease-realted deaths.
#' @param births.FUN Module function for births or entries into the population.
#' @param test.FUN Module function for diagnostic disease testing.
#' @param tx.FUN Module function for ART initiation and adherence.
#' @param prep.FUN Module function for PrEP initiation and utilization.
#' @param progress.FUN Module function for HIV disease progression.
#' @param vl.FUN Module function for HIV viral load evolution.
#' @param aiclass.FUN Module function for one-off AI risk class transitions.
#' @param roleclass.FUN Module function for transitions in sexual roles.
#' @param resim_nets.FUN Module function for network resimulation at each time
#'        step.
#' @param disclose.FUN Module function for HIV status disclosure.
#' @param acts.FUN Module function to simulate the number of sexual acts within
#'        partnerships.
#' @param condoms.FUN Module function to simulate condom use within acts.
#' @param position.FUN Module function to simulate sexual position within acts.
#' @param trans.FUN Module function to stochastically simulate HIV transmission
#'        over acts given individual and dyadic attributes.
#' @param stitrans.FUN Module function to simulate GC/CT transmission over current
#'        edgelist.
#' @param stirecov.FUN Module function to simulate recovery from GC/CT, heterogeneous
#'        by disease, site, symptoms, and treatment status.
#' @param stitx.FUN Module function to simulate treatment of GC/CT.
#' @param prev.FUN Module function to calculate prevalence summary statistics.
#' @param verbose.FUN Module function to print model progress to the console or
#'        external text files.
#' @param save.nwstats Calculate and save network statistics as defined in the
#'        \code{simnet} modules.
#' @param verbose If \code{TRUE}, print out simulation progress to the console
#'        if in interactive mode or text files if in batch mode.
#' @param verbose.int Integer specifying the interval between time steps at which
#'        progress is printed.
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{control_msm}, which can be passed to the
#' EpiModel function \code{netsim}.
#'
#' @keywords msm
#'
#' @export
control_msm_whamp <- function(simno = 1,
                        nsims = 1,
                        ncores = 1,
                        nsteps = 100,
                        start = 1,
                        initialize.FUN = initialize_msm_whamp,
                        aging.FUN = aging_msm,
                        deaths.FUN = deaths_msm_whamp,
                        births.FUN = births_msm_whamp,
                        test.FUN = test_msm_whamp,
                        tx.FUN = tx_msm_whamp,
                        progress.FUN = progress_msm_whamp,
                        vl.FUN = vl_msm_whamp,
                        aiclass.FUN = update_aiclass_msm_whamp,
                        roleclass.FUN = NULL,
                        resim_nets.FUN = simnet_msm_whamp,
                        disclose.FUN = disclose_msm_whamp,
                        nwfeatures.FUN = nwfeatures_msm_whamp,
                        acts.FUN = acts_msm_whamp,
                        condoms.FUN = condoms_msm_whamp,
                        position.FUN = position_msm_whamp,
                        prep.FUN = prep_msm_whamp,
                        trans.FUN = trans_msm_whamp,
                        stitrans.FUN = NULL,
                        stirecov.FUN = NULL,
                        stitx.FUN = NULL,
                        prev.FUN = prevalence_msm_whamp,
                        verbose.FUN = verbose_msm,
                        save.nwstats = FALSE,
                        verbose = TRUE,
                        verbose.int = 1,
                        ...) {

  formal.args <- formals(sys.function())
  dot.args <- list(...)
  p <- get_args(formal.args, dot.args)

  p$skip.check <- TRUE
  p$save.transmat <- FALSE

  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  bi.mods <- bi.mods[which(sapply(bi.mods, function(x) !is.null(eval(parse(text = x))),
                                  USE.NAMES = FALSE) == TRUE)]
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)

  p$save.other = c("attr", "temp", "el", "p", "cel.temp", "cel.complete")

  p$save.network = FALSE

  class(p) <- "control.net"
  return(p)
}



# HET -----------------------------------------------------------------


#' @title Parameters for Stochastic Network Model of HIV-1 Infection in
#'        Sub-Saharan Africa
#'
#' @description Sets the simulation parameters for the stochastic
#'              network model of HIV-1 Infection among Heterosexuals in
#'              Sub-Saharan Africa for the \code{EpiModelHIV} package.
#'
#' @param time.unit Unit of time relative to one day.
#'
#' @param acute.stage.mult Acute stage multiplier for increased infectiousness
#'        above impact of heightened viral load.
#' @param aids.stage.mult AIDS stage multiplier for increased infectiousness in
#'        AIDS above impact of heightened viral load.
#'
#' @param vl.acute.topeak Time in days to peak viremia during acute infection.
#' @param vl.acute.toset Time in days to viral set point following peak viremia.
#' @param vl.acute.peak Log 10 viral load at acute peak.
#' @param vl.setpoint Log 10 viral load at set point.
#' @param vl.aidsmax Maximum log 10 viral load during AIDS.
#'
#' @param cond.prob Probability of condoms per act with partners.
#' @param cond.eff Efficacy of condoms per act in HIV prevention.
#'
#' @param act.rate.early Daily per-partnership act rate in early disease.
#' @param act.rate.late Daily per-partnership act rate in late disease.
#' @param act.rate.cd4 CD4 count at which the \code{act.rate.late} applies.
#' @param acts.rand If \code{TRUE}, will draw number of total and unprotected
#'        acts from a binomial distribution parameterized by the \code{act.rate}.
#'
#' @param circ.prob.birth Proportion of men circumcised at birth.
#' @param circ.eff Efficacy of circumcision per act in HIV prevention.
#'
#' @param tx.elig.cd4 CD4 count at which a person becomes eligible for treatment.
#' @param tx.init.cd4.mean Mean CD4 count at which person presents for care.
#' @param tx.init.cd4.sd SD of CD4 count at which person presents for care.
#' @param tx.adhere.full Proportion of people who start treatment who are fully
#'        adherent.
#' @param tx.adhere.part Of the not fully adherent proportion, the percent of time
#'        they are on medication.
#' @param tx.vlsupp.time Time in days from treatment initiation to viral suppression.
#' @param tx.vlsupp.level Log 10 viral load level at suppression.
#' @param tx.cd4.recrat.feml Rate of CD4 recovery under treatment for males.
#' @param tx.cd4.recrat.male Rate of CD4 recovery under treatment for females.
#' @param tx.cd4.decrat.feml Rate of CD4 decline under periods of non-adherence
#'        for females.
#' @param tx.cd4.decrat.male Rate of CD4 decline under periods of non-adherence
#'        for males.
#' @param tx.coverage Proportion of treatment-eligible persons who have initiated
#'        treatment.
#' @param tx.prev.eff Proportional amount by which treatment reduces infectivity
#'        of infected partner.
#'
#' @param b.rate General entry rate per day for males and females specified.
#' @param b.rate.method Method for assigning birth rates, with options of "totpop"
#'        for births as a function of the total population size, "fpop" for births
#'        as a function of the female population size, and "stgrowth" for a constant
#'        stable growth rate.
#' @param b.propmale Proportion of entries assigned as male. If NULL, then set
#'        adaptively based on the proportion at time 1.
#'
#' @param ds.exit.age Age at which the age-specific ds.rate is set to 1, with NA
#'        value indicating no censoring.
#' @param ds.rate.mult Simple multiplier for background death rates.
#' @param di.cd4.aids CD4 count at which late-stage AIDS occurs and the risk of
#'        mortality is governed by \code{di.cd4.rate}.
#' @param di.cd4.rate Mortality in late-stage AIDS after hitting a nadir CD4 of
#'        \code{di.cd4.aids}.
#' @param ... additional arguments to be passed into model.
#'
#' @details This function sets the parameters for the models.
#'
#' @keywords het
#'
#' @export
#'
param_het <- function(time.unit = 7,

                      acute.stage.mult = 5,
                      aids.stage.mult = 1,

                      vl.acute.topeak = 14,
                      vl.acute.toset = 107,
                      vl.acute.peak = 6.7,
                      vl.setpoint = 4.5,
                      vl.aidsmax = 7,

                      cond.prob = 0.09,
                      cond.eff = 0.78,

                      act.rate.early = 0.362,
                      act.rate.late = 0.197,
                      act.rate.cd4 = 50,
                      acts.rand = TRUE,

                      circ.prob.birth = 0.9,
                      circ.eff = 0.53,

                      tx.elig.cd4 = 350,
                      tx.init.cd4.mean = 120,
                      tx.init.cd4.sd = 40,
                      tx.adhere.full = 0.76,
                      tx.adhere.part = 0.50,
                      tx.vlsupp.time = 365/3,
                      tx.vlsupp.level = 1.5,
                      tx.cd4.recrat.feml = 11.6/30,
                      tx.cd4.recrat.male = 9.75/30,
                      tx.cd4.decrat.feml = 11.6/30,
                      tx.cd4.decrat.male = 9.75/30,
                      tx.coverage = 0.3,
                      tx.prev.eff = 0.96,

                      b.rate = 0.03/365,
                      b.rate.method = "totpop",
                      b.propmale = NULL,

                      ds.exit.age = 55,
                      ds.rate.mult = 1,
                      di.cd4.aids = 50,
                      di.cd4.rate = 2/365,
                      ...) {

  ## Process parameters
  p <- list()
  formal.args <- formals(sys.function())
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg))
  }
  dot.args <- list(...)
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }


  ## trans.rate multiplier
  p$trans.rate <- p$trans.rate * p$trans.rate.mult


  ## Death rate transformations
  ltGhana <- EpiModelHIV::ltGhana
  ds.rates <- ltGhana[ltGhana$year == 2011, ]
  ds.rates$mrate <- ds.rates$mrate / 365
  if (is.numeric(ds.exit.age)) {
    ds.rates$mrate[ds.rates$agStart >= ds.exit.age] <- 1
  }
  ds.rates$reps <- ds.rates$agEnd - ds.rates$agStart + 1
  ds.rates$reps[ds.rates$agStart == 100] <- 1
  male <- rep(ds.rates$male, ds.rates$reps)
  mrate <- rep(ds.rates$mrate, ds.rates$reps)
  mrate <- pmin(1, mrate * ds.rate.mult)
  age <- rep(0:100, 2)
  ds.rates <- data.frame(male = male, age, mrate = mrate)
  ds.rates <- ds.rates[ds.rates$age != 0, ]
  p$ds.rates <- ds.rates

  ## Time unit scaling
  if (time.unit > 1) {

    ## Rates multiplied by time unit
    p$act.rate.early <- act.rate.early * time.unit
    p$act.rate.late <- act.rate.late * time.unit
    p$b.rate <- b.rate * time.unit
    p$ds.rates$mrate <- ifelse(p$ds.rates$mrate < 1,
                               p$ds.rates$mrate * time.unit,
                               p$ds.rates$mrate)

    p$dx.prob.feml <- p$dx.prob.feml * time.unit
    p$dx.prob.male <- p$dx.prob.male * time.unit
    p$tx.cd4.recrat.feml <- tx.cd4.recrat.feml * time.unit
    p$tx.cd4.recrat.male <- tx.cd4.recrat.male * time.unit
    p$tx.cd4.decrat.feml <- tx.cd4.decrat.feml * time.unit
    p$tx.cd4.decrat.male <- tx.cd4.decrat.male * time.unit
    p$di.cd4.rate <- di.cd4.rate * time.unit

    ## Intervals divided by time unit
    p$vl.acute.topeak <- vl.acute.topeak / time.unit
    p$vl.acute.toset <- vl.acute.toset / time.unit

    p$tx.vlsupp.time <- tx.vlsupp.time / time.unit

  }

  p$model <- "a2"

  class(p) <- "param.net"
  return(p)
}


#' @title Initial Conditions for Stochastic Network Model of HIV-1 Infection in
#'        Sub-Saharan Africa
#'
#' @description This function sets the initial conditions for the stochastic
#'              network models in the \code{epimethods} package.
#'
#' @param i.prev.male Prevalence of initially infected males.
#' @param i.prev.feml Prevalence of initially infected females.
#' @param ages.male initial ages of males in the population.
#' @param ages.feml initial ages of females in the population.
#' @param inf.time.dist Probability distribution for setting time of infection
#'        for nodes infected at T1, with options of \code{"geometric"} for randomly
#'        distributed on a geometric distribution with a probability of the
#'        reciprocal of the average length of infection, \code{"uniform"} for a
#'        uniformly distributed time over that same interval, or \code{"allacute"} for
#'        placing all infections in the acute stage at the start.
#' @param max.inf.time Maximum infection time in days for infection at initialization,
#'        used when \code{inf.time.dist} is \code{"geometric"} or \code{"uniform"}.
#' @param ... additional arguments to be passed into model.
#'
#' @details This function sets the initial conditions for the models.
#'
#' @keywords het
#'
#' @export
#'
init_het <- function(i.prev.male = 0.05,
                     i.prev.feml = 0.05,
                     ages.male = seq(18, 55, 7/365),
                     ages.feml = seq(18, 55, 7/365),
                     inf.time.dist = "geometric",
                     max.inf.time = 5 * 365,
                     ...) {

  ## Process parameters
  p <- list()
  formal.args <- formals(sys.function())
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg))
  }
  dot.args <- list(...)
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }


  ## Parameter checks
  if (!(inf.time.dist %in% c("uniform", "geometric", "allacute"))) {
    stop("inf.time.dist must be \"uniform\" or \"geometric\" or \"allacute\" ")
  }

  class(p) <- "init.net"
  return(p)
}


#' @title Control Settings for Stochastic Network Model of HIV-1 Infection in
#'        Sub-Saharan Africa
#'
#' @description This function sets the control settings for the stochastic
#'              network models in the \code{epimethods} package.
#'
#' @param simno Simulation ID number.
#' @param nsteps Number of time steps to simulate the model over in whatever unit
#'        implied by \code{time.unit}.
#' @param start Starting time step for simulation
#' @param nsims Number of simulations.
#' @param ncores Number of parallel cores to use for simulation jobs, if using
#'        the \code{EpiModel.hpc} package.
#' @param par.type Parallelization type, either of \code{"single"} for multi-core
#'        or \code{"mpi"} for multi-node MPI threads.
#' @param initialize.FUN Module to initialize the model at time 1.
#' @param aging.FUN Module to age active nodes.
#' @param cd4.FUN CD4 progression module.
#' @param vl.FUN HIV viral load progression module.
#' @param dx.FUN HIV diagnosis module.
#' @param tx.FUN HIV treatment module.
#' @param deaths.FUN Module to simulate death or exit.
#' @param births.FUN Module to simulate births or entries.
#' @param resim_nets.FUN Module to resimulate the network at each time step.
#' @param trans.FUN Module to simulate disease infection.
#' @param prev.FUN Module to calculate disease prevalence at each time step,
#'        with the default function of \code{\link{prevalence_het}}.
#' @param verbose.FUN Module to print simulation progress to screen, with the
#'        default function of \code{\link{verbose_het}}.
#' @param module.order A character vector of module names that lists modules the
#'        order in which they should be evaluated within each time step. If
#'        \code{NULL}, the modules will be evaluated as follows: first any
#'        new modules supplied through \code{...} in the order in which they are
#'        listed, then the built-in modules in their order of the function listing.
#'        The \code{initialize.FUN} will always be run first and the
#'        \code{verbose.FUN} always last.
#' @param save.nwstats Save out network statistics.
#' @param save.other Other list elements of dat to save out.
#' @param verbose If \code{TRUE}, print progress to console.
#' @param verbose.int Interval for printing progress to console.
#' @param skip.check If \code{TRUE}, skips the error check for parameter values,
#'        initial conditions, and control settings before running the models.
#' @param ... Additional arguments passed to the function.
#'
#' @details This function sets the parameters for the models.
#'
#' @keywords het
#'
#' @export
#'
control_het <- function(simno = 1,
                        nsteps = 100,
                        start = 1,
                        nsims = 1,
                        ncores = 1,
                        par.type = "single",
                        initialize.FUN = initialize_het,
                        aging.FUN = aging_het,
                        cd4.FUN = cd4_het,
                        vl.FUN = vl_het,
                        dx.FUN = dx_het,
                        tx.FUN = tx_het,
                        deaths.FUN = deaths_het,
                        births.FUN = births_het,
                        resim_nets.FUN = simnet_het,
                        trans.FUN = trans_het,
                        prev.FUN = prevalence_het,
                        verbose.FUN = verbose_het,
                        module.order = NULL,
                        save.nwstats = FALSE,
                        save.other = c("el", "attr"),
                        verbose = TRUE,
                        verbose.int = 1,
                        skip.check = TRUE,
                        ...) {

  p <- list()
  formal.args <- formals(sys.function())
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg))
  }
  dot.args <- list(...)
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }

  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  bi.mods <- bi.mods[which(sapply(bi.mods, function(x) !is.null(eval(parse(text = x))),
                                  USE.NAMES = FALSE) == TRUE)]
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names.dot.args, value = TRUE)

  p$save.transmat <- FALSE
  p$save.network <- FALSE

  class(p) <- "control.net"
  return(p)
}
