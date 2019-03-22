#' @title Network diagnostics.
#'
#' @description Module function for tracking features of the sexual newtork over time.
#'
#' @inheritParams aging_msm
#'
#' @details
#' At each time step, the cross-network degree distribution is calculated and stored, along with
#' statistics about race, region, and age mixing, and the duration of active ties. 
#'
#' @return
#' This function returns the \code{dat} object with the updated partnership list for diagnostics
#' on \code{dat$plist}.
#'
#' @keywords module WHAMP network diagnostics
#' @export
#'
nwfeatures_msm_whamp <- function(dat, at){
  
  if(dat$param$nw_track == TRUE){
    
    # Variables --------------------------------------------------------------
    
    # Attributes
    
    uid <- dat$attr$uid
    race..wa <- dat$attr$race..wa
    region <- dat$attr$region
    age <- floor(dat$attr$age)
    riskg <- dat$attr$riskg
    role <- dat$attr$role.class
    deg.main <- dat$attr$deg.main
    deg.pers <- dat$attr$deg.pers
    cel.temp <- dat$cel.temp
    cel.complete <- dat$cel.complete
    
    # Build edgelist with attributes
    el1 <- as.data.frame(dat$el[[1]])
    el1[,3] <- apply(el1[,c(1,2)], 1, FUN=max)
    el1[,4] <- apply(el1[,c(1,2)], 1, FUN=min)
    el1[,5]<-rep("main",length(el1[,1]))
    el1 <- el1[,3:5]
    
    el2 <- as.data.frame(dat$el[[2]])
    el2[,3] <- apply(el2[,c(1,2)], 1, FUN=max)
    el2[,4] <- apply(el2[,c(1,2)], 1, FUN=min)
    el2[,5]<-rep("pers",length(el2[,1]))
    el2 <- el2[,3:5]
    
    el3 <- el3 <- as.data.frame(dat$el[[3]])
    el3[,3] <- apply(el3[,c(1,2)], 1, FUN=max)
    el3[,4] <- apply(el3[,c(1,2)], 1, FUN=min)
    el3[,5]<-rep("inst",length(el3[,1]))
    el3 <- el3[ ,3:5]
    
    cel <- rbind(el1, el2, el3)
    cel[,4] <- race..wa[cel[,1]]
    cel[,5] <- race..wa[cel[,2]]
    cel[,6] <- region[cel[,1]]
    cel[,7] <- region[cel[,2]]
    cel[,8] <- age[cel[,1]]
    cel[,9] <- age[cel[,2]]
    cel[,10] <- riskg[cel[,1]]
    cel[,11] <- riskg[cel[,2]]
    cel[,12] <- role[cel[,1]]
    cel[,13] <- role[cel[,2]]
    
    ## Convert node IDs to UID.
    cel[,1] <- uid[cel[,1]]
    cel[,2] <- uid[cel[,2]]
    
    ##Assign rel IDs using UIDs    
    cel[,14] <- (cel[,1] * 1000000000) + cel[,2]
    colnames(cel) <- c("p1", "p2", "type", "race1","race2","region1","region2","age1","age2",
                       "riskg1","riskg2","role1","role2", "ID")
    
    ##Calculate degree distributions
    dat$epi$deg.main[at] <- mean(deg.main)
    dat$epi$deg.pers[at] <- mean(deg.pers)
    dat$epi$deg.inst[at] <- mean(get_degree(dat$el[[3]]))
    dat$epi$main0pers0[at] <- prop.table(table(dat$attr$deg.main, dat$attr$deg.pers))[1,1]
    dat$epi$main0pers1[at] <- prop.table(table(dat$attr$deg.main, dat$attr$deg.pers))[1,2]
    dat$epi$main0pers2[at] <- prop.table(table(dat$attr$deg.main, dat$attr$deg.pers))[1,3]
    dat$epi$main1pers0[at] <- prop.table(table(dat$attr$deg.main, dat$attr$deg.pers))[2,1]
    dat$epi$main1pers1[at] <- prop.table(table(dat$attr$deg.main, dat$attr$deg.pers))[2,2]
    dat$epi$main1pers2[at] <- prop.table(table(dat$attr$deg.main, dat$attr$deg.pers))[2,3]
    
    # Calculate the mean number of inst partnerships for groups defined by main and pers partnerships
    n.inst <- get_degree(dat$el[[3]])
    dat$epi$meaninst.00[at] <- mean(n.inst[deg.main == 0 & deg.pers == 0])
    dat$epi$meaninst.01[at] <- mean(n.inst[deg.main == 0 & deg.pers == 1])
    dat$epi$meaninst.02[at] <- mean(n.inst[deg.main == 0 & deg.pers == 2])
    dat$epi$meaninst.10[at] <- mean(n.inst[deg.main == 1 & deg.pers == 0])
    dat$epi$meaninst.11[at] <- mean(n.inst[deg.main == 1 & deg.pers == 1])
    dat$epi$meaninst.12[at] <- mean(n.inst[deg.main == 1 & deg.pers == 2])
    
    
    ##Calculate homophily and role mixing stats
    cel.m <- cel[cel$type %in% "main", ]
    cel.p <- cel[cel$type %in% "pers", ]
    cel.i <- cel[cel$type %in% "inst", ]
    
    # Main
    dat$epi$prop.hom.m.B[at] <- sum(cel.m$race1 %in% "B" & cel.m$race2 %in% "B") / 
                          sum(cel.m$race1 %in% "B" | cel.m$race2 %in% "B")
    dat$epi$prop.hom.m.H[at] <- sum(cel.m$race1 %in% "H" & cel.m$race2 %in% "H") / 
                          sum(cel.m$race1 %in% "H" | cel.m$race2 %in% "H")
    dat$epi$prop.hom.m.O[at] <- sum(cel.m$race1 %in% "O" & cel.m$race2 %in% "O") / 
                          sum(cel.m$race1 %in% "O" | cel.m$race2 %in% "O")
    
    dat$epi$prop.hom.m.KC[at] <- sum(cel.m$region1 %in% "KC" & cel.m$region2 %in% "KC") / 
                                      sum(cel.m$region1 %in% "KC" | cel.m$region2 %in% "KC")
    dat$epi$prop.hom.m.OW[at] <- sum(cel.m$region1 %in% "OW" & cel.m$region2 %in% "OW") / 
                                      sum(cel.m$region1 %in% "OW" | cel.m$region2 %in% "OW")
    dat$epi$prop.hom.m.EW[at] <- sum(cel.m$region1 %in% "EW" & cel.m$region2 %in% "EW") / 
                                      sum(cel.m$region1 %in% "EW" | cel.m$region2 %in% "EW")
    
    dat$epi$sqrtadiff.m[at] <- mean(abs(sqrt(cel.m$age1) - sqrt(cel.m$age2)))
    
    dat$epi$badroles.m[at] <- sum((cel.m$role1 %in% "R" & cel.m$role2 %in% "R") | 
                                    (cel.m$role1 %in% "I" & cel.m$role2 %in% "I"))
    
    # Pers
    dat$epi$prop.hom.p.B[at] <- sum(cel.p$race1 %in% "B" & cel.p$race2 %in% "B") / 
      sum(cel.p$race1 %in% "B" | cel.p$race2 %in% "B")
    dat$epi$prop.hom.p.H[at] <- sum(cel.p$race1 %in% "H" & cel.p$race2 %in% "H") / 
      sum(cel.p$race1 %in% "H" | cel.p$race2 %in% "H")
    dat$epi$prop.hom.p.O[at] <- sum(cel.p$race1 %in% "O" & cel.p$race2 %in% "O") / 
      sum(cel.p$race1 %in% "O" | cel.p$race2 %in% "O")
    
    dat$epi$prop.hom.p.KC[at] <- sum(cel.p$region1 %in% "KC" & cel.p$region2 %in% "KC") / 
      sum(cel.p$region1 %in% "KC" | cel.p$region2 %in% "KC")
    dat$epi$prop.hom.p.OW[at] <- sum(cel.p$region1 %in% "OW" & cel.p$region2 %in% "OW") / 
      sum(cel.p$region1 %in% "OW" | cel.p$region2 %in% "OW")
    dat$epi$prop.hom.p.EW[at] <- sum(cel.p$region1 %in% "EW" & cel.p$region2 %in% "EW") / 
      sum(cel.p$region1 %in% "EW" | cel.p$region2 %in% "EW")
    
    dat$epi$sqrtadiff.p[at] <- mean(abs(sqrt(cel.p$age1) - sqrt(cel.p$age2)))
    
    dat$epi$badroles.p[at] <- sum((cel.p$role1 %in% "R" & cel.p$role2 %in% "R") | 
                                    (cel.p$role1 %in% "I" & cel.p$role2 %in% "I"))
    
    # Inst
    dat$epi$prop.hom.i.B[at] <- sum(cel.i$race1 %in% "B" & cel.i$race2 %in% "B") / 
      sum(cel.i$race1 %in% "B" | cel.i$race2 %in% "B")
    dat$epi$prop.hom.i.H[at] <- sum(cel.i$race1 %in% "H" & cel.i$race2 %in% "H") / 
      sum(cel.i$race1 %in% "H" | cel.i$race2 %in% "H")
    dat$epi$prop.hom.i.O[at] <- sum(cel.i$race1 %in% "O" & cel.i$race2 %in% "O") / 
      sum(cel.i$race1 %in% "O" | cel.i$race2 %in% "O")
    
    dat$epi$prop.hom.i.KC[at] <- sum(cel.i$region1 %in% "KC" & cel.i$region2 %in% "KC") / 
      sum(cel.i$region1 %in% "KC" | cel.i$region2 %in% "KC")
    dat$epi$prop.hom.i.OW[at] <- sum(cel.i$region1 %in% "OW" & cel.i$region2 %in% "OW") / 
      sum(cel.i$region1 %in% "OW" | cel.i$region2 %in% "OW")
    dat$epi$prop.hom.i.EW[at] <- sum(cel.i$region1 %in% "EW" & cel.i$region2 %in% "EW") / 
      sum(cel.i$region1 %in% "EW" | cel.i$region2 %in% "EW")
    
    dat$epi$sqrtadiff.i[at] <- mean(abs(sqrt(cel.i$age1) - sqrt(cel.i$age2)))
    
    dat$epi$badroles.i[at] <- sum((cel.i$role1 %in% "R" & cel.i$role2 %in% "R") | 
                                    (cel.i$role1 %in% "I" & cel.i$role2 %in% "I"))
    
    ##Count duplicates (main and pers)
    relids.mp <- cel[cel$type %in% c("main", "pers"), 14] # main and persistent partnerships
    trelcount <- length(relids.mp)
    urelcount <- length(unique(relids.mp))
    dup_count <- trelcount - urelcount
    dat$epi$duplicates.mp[at] <- dup_count
    
    ##Count duplicates including inst
    relids.all <- cel[, 14]
    trelcount <- length(relids.all)
    urelcount <- length(unique(relids.all))
    dup_count <- trelcount - urelcount
    dat$epi$duplicates.mpi[at] <- dup_count
    
    ##drop duplicates for counting concurrency durations (their presence interferes with the merge)
    ##Low counts will have a finite inpact of duration calculations.
    drop <- which(duplicated(relids.all)==TRUE) 
    if (length(drop > 0)) {
      cel <- cel[-drop,]
      }
    
    ##Track rel start and end for main and pers
    if (at == 2) {
      cel.temp <- cel[cel$type %in% c("main", "pers"), ]
      cel.temp[,15] <- rep(at, length(cel.temp[,1])) 
      cel.temp[,16] <- rep(NA, length(cel.temp[,1])) 
      colnames(cel.temp) <- c("p1", "p2", "type", "race1","race2","region1","region2","age1","age2",
                            "riskg1","riskg2","role1","role2", "ID", "start", "end")
      cel.complete <- cel.temp[0,]
    }
    
    if(at > 2){
      ##Get the rels from cel.temp that are not in cel.
      ids.old <- unique(cel.temp[,"ID"])
      
      ids.new <- NULL
      for(i in 1:length(dat$temp$new.edges[,1])){
        ids.new[i] <- (max(dat$attr$uid[dat$temp$new.edges[i,1]], dat$attr$uid[dat$temp$new.edges[i,2]]) * 1000000000) 
        + min(dat$attr$uid[dat$temp$new.edges[i,1]], dat$attr$uid[dat$temp$new.edges[i,2]])
      }
      
      ids.cur <- unique(cel[,"ID"])
      ids.ongoing <- intersect(ids.old,ids.cur)
      ids.ended <- setdiff(ids.old,ids.ongoing) 
      
      ended <- cel.temp[cel.temp$ID %in% ids.ended,]
       if (nrow(ended) > 0) {
         ended$end <- at
       }
      cel.complete <- rbind(cel.complete, ended)
  
      
      cel.temp <- merge(x = cel, y = cel.temp, by = c("p1","p2"), all.x = TRUE)
      cel.temp <- cel.temp[,c(1:14,27,28)]
      
      
    }
    
    
    colnames(cel.temp) <- c("p1", "p2", "type", "race1","race2","region1","region2","age1","age2",
                            "riskg1","riskg2","role1","role2", "ID", "start", "end")
    new <- which(is.na(cel.temp$start) == TRUE)
    cel.temp$start[new] <- at
    
    dat$cel.temp <- cel.temp
    dat$cel.complete <- cel.complete
  }
  
  
  return(dat)
}