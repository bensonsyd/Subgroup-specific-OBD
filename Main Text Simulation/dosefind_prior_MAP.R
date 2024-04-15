##################################
# functions to find doses
##################################

## rule-based run-in dose finding
runin <- function(doses, FullDat, one.group){
  dose1<- dose0 <- doses[1]
  if (nrow(FullDat)==4 & sum(FullDat$DLT, na.rm=TRUE)== 0) {dose1<- dose0 <- doses[2]
  }else if (nrow(FullDat)==6 & sum(FullDat$DLT, na.rm=TRUE)== 1) {dose1<- dose0 <- doses[2]
  }else if (nrow(FullDat)==8){
    if(one.group == F){
      if (sum(FullDat$DLT[FullDat$Subgroup==0], na.rm=TRUE)== 0) {dose0<-doses[2]}
      if (sum(FullDat$DLT[FullDat$Subgroup==1], na.rm=TRUE)== 0) {dose1<-doses[2]}
    }
  }else if (nrow(FullDat)==12){
    if(one.group == F){
      if (sum(FullDat$DLT[FullDat$Subgroup==0], na.rm=TRUE)== 1) {dose0<-doses[2]}
      if (sum(FullDat$DLT[FullDat$Subgroup==1], na.rm=TRUE)== 1) {dose1<-doses[2]}
    }else if(one.group == T){
      if(sum(FullDat$DLT, na.rm=TRUE) == 2){dose1 <- dose0 <- doses[2]}
    }
  }
  return(c(dose0, dose1))
}

## model-based dose finding for two groups
dose.find.twogroup <- function(FullDat, Design.eff, Design.eff.c, PriorSpec.tox, toxicity_target, PriorSpec.eff, tox.threshold, eff.threshold, doses, res, inclusion.threshold, tox.level, eff.level, u00, u11){
  posterior.tox <- logit.spike(DLT~Dose.c + Subgroup.c + Dose.c:Subgroup.c, 
                               prior = PriorSpec.tox, data=FullDat, niter = 20000, na.action = na.exclude, ping = 0)
  pr.tox<-predict.logit.spike(posterior.tox, newdata = as.matrix(Design.eff.c)[,c(1,2,4,5)], 
                              burn=5000)
  admiss.tox <- rowMeans(pr.tox<tox.threshold)>tox.level
  #pr.tox <- rowMeans(pr.tox)
  tox.map <- rep(NA, nrow(pr.tox))
  for(i in 1:nrow(pr.tox)){
    tox.map[i] <- map_estimate(pr.tox[i,])
  }
  pr.tox <- tox.map
  
  if(!is.null(res)){
    res$pip.alpha2 <- summary(posterior.tox, burn=5000, order=F)$coefficients[rownames(summary(posterior.tox, burn=5000, order=F)$coefficients)=='Subgroup.c'][5]
    res$pip.alpha3 <- summary(posterior.tox, burn=5000, order=F)$coefficients[rownames(summary(posterior.tox, burn=5000, order=F)$coefficients)=='Dose.c:Subgroup.c'][5]
    res$pip.subgroup.tox <- mean(posterior.tox$beta[5001:20000,3] != 0 | posterior.tox$beta[5001:20000,4] != 0)
  }
  
  if (is.null(toxicity_target)){
    posterior.eff <- logit.spike(eff~ Dose.c + I(Dose.c^2) + Subgroup.c + Dose.c:Subgroup.c + I(Dose.c^2):Subgroup.c, 
                                 prior = PriorSpec.eff, data=FullDat, niter = 20000, na.action = na.exclude, ping = 0)
    pr.eff<-predict.logit.spike(posterior.eff, newdata = as.matrix(Design.eff.c), burn = 5000)
    admiss.eff <- rowMeans(pr.eff>eff.threshold)>eff.level
    #pr.eff <- rowMeans(pr.eff)
    eff.map <- rep(NA, nrow(pr.eff))
    for(i in 1:nrow(pr.eff)){
      eff.map[i] <- map_estimate(pr.eff[i,])
    }
    pr.eff <- eff.map
    if(!is.null(res)){
      res$pip.beta3 <- summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='Subgroup.c'][5]
      res$pip.beta4 <- summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='Dose.c:Subgroup.c'][5]
      res$pip.beta5 <- summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='I(Dose.c^2):Subgroup.c'][5]
      res$pip.subgroup.eff <- mean(posterior.eff$beta[5001:20000,4] != 0 | posterior.eff$beta[5001:20000,5] != 0 | posterior.eff$beta[5001:20000,6] != 0)
    }
    
    utility <- utility_group(u00=u00, u11=u11, ptox=pr.tox, peff=pr.eff)
    
    # A dose must meet the toxicity threshold and EITHER meet the efficacy threshold or not have been tried yet
    if(!is.null(res)){
      # Group Z=0
      dose0_options <-
        which(admiss.tox[Design.eff$Subgroup==0]==T & admiss.eff[Design.eff$Subgroup==0]==T & res[3:(2+length(doses))] >0)
      res$obd0.bayes <- ifelse(length(dose0_options) == 0, 0, doses[dose0_options][which.max(utility[Design.eff$Subgroup==0][dose0_options])])
      
      # Group Z=1
      dose1_options <-
        which(admiss.tox[Design.eff$Subgroup==1]==T & admiss.eff[Design.eff$Subgroup==1]==T & res[(2+length(doses)+1):(2+length(doses)+length(doses))] >0)
      res$obd1.bayes <- ifelse(length(dose1_options) == 0, 0, doses[dose1_options][which.max(utility[Design.eff$Subgroup==1][dose1_options])])
    }else{
      # Group Z=0
      dose0_options <- which(admiss.tox[Design.eff$Subgroup==0]==T & 
                               (admiss.eff[Design.eff$Subgroup==0]==T | !(doses %in% FullDat$Dose[FullDat$Subgroup==0])))
      if(length(dose0_options) != 0){
        dose0 <- doses[dose0_options][which.max(utility[Design.eff$Subgroup==0][dose0_options])]
      }else{
        dose0 <- NA
        print("no dose available for Z=0")
        print(pr.tox[Design.eff$Subgroup==0])
        print(pr.eff[Design.eff$Subgroup==0])
      }
      
      # Group Z=1
      dose1_options <- which(admiss.tox[Design.eff$Subgroup==1]==T & 
                               (admiss.eff[Design.eff$Subgroup==1]==T | !(doses %in% FullDat$Dose[FullDat$Subgroup==1])))
      if(length(dose1_options) != 0){
        dose1 <- doses[dose1_options][which.max(utility[Design.eff$Subgroup==1][dose1_options])]
      }else{
        dose1 <- NA
        print("no dose available for Z=1")
        print(pr.tox[Design.eff$Subgroup==1])
        print(pr.eff[Design.eff$Subgroup==1])
      }
    }
  }
  
  #### If you want to target a certain toxicity rate (like Cotteril and Jaki)
  if (!is.null(toxicity_target)){
    if(!is.null(res)){
      # Group Z=0
      dose0_options <- which(admiss.tox[Design.eff$Subgroup==0]==T & res[3:(2+length(doses))] >0)
      res$ttd0.bayes <- ifelse(length(dose0_options)==0, 0, doses[dose0_options][which.min(abs(pr.tox[Design.eff$Subgroup==0][dose0_options]-toxicity_target))])
      
      # Group Z=1
      dose1_options <- which(admiss.tox[Design.eff$Subgroup==1]==T & res[(2+length(doses)+1):(2+length(doses)+length(doses))] >0)
      res$ttd1.bayes <- ifelse(length(dose1_options)==0, 0, doses[dose1_options][which.min(abs(pr.tox[Design.eff$Subgroup==1][dose1_options]-toxicity_target))])
    }else{
      # Group Z=0
      dose0_options <- which(admiss.tox[Design.eff$Subgroup==0]==T)
      if(length(dose0_options) != 0){
        dose0 <- doses[dose0_options][which.min(abs(pr.tox[Design.eff$Subgroup==0][dose0_options]-toxicity_target))]
      }else{
        dose0 <- NA
        print("no dose available for Z=0")
        print(pr.tox[Design.eff$Subgroup==0])
      }
      
      # Group Z=1
      dose1_options <- which(admiss.tox[Design.eff$Subgroup==1]==T)
      if(length(dose1_options) != 0){
        dose1 <- doses[dose1_options][which.min(abs(pr.tox[Design.eff$Subgroup==1][dose1_options]-toxicity_target))]
      }else{
        dose1 <- NA
        print("no dose available for Z=1")
        print(pr.tox[Design.eff$Subgroup==1])
      }
    }
  }
  
  if(is.null(res)){
    # Don't allow escalation by more than 1 step above tried doses 
    if(length(dose0_options) != 0 & length(FullDat$Dose[FullDat$Subgroup == 0]) > 0){
      if (max(FullDat$Dose[FullDat$Subgroup == 0], na.rm=TRUE) != doses[length(doses)] & 
          dose0 > doses[which(doses == max(FullDat$Dose[FullDat$Subgroup == 0], na.rm=TRUE))+1]){ 
        dose0 <- doses[which(doses == max(FullDat$Dose[FullDat$Subgroup == 0], na.rm=TRUE))+1]
      }
    }
    if(length(dose1_options) != 0 & length(FullDat$Dose[FullDat$Subgroup == 1])>0){
      if (max(FullDat$Dose[FullDat$Subgroup == 1], na.rm=TRUE) != doses[length(doses)] & 
          dose1 > doses[which(doses == max(FullDat$Dose[FullDat$Subgroup == 1], na.rm=TRUE))+1]){
        dose1 <- doses[which(doses == max(FullDat$Dose[FullDat$Subgroup == 1], na.rm=TRUE))+1]
      }
    }
    return(c(dose0, dose1))
  }else{
    res$alpha_0.bayes <-   summary(posterior.tox, burn=5000, order=F)$coefficients[rownames(summary(posterior.tox, burn=5000, order=F)$coefficients)=='(Intercept)'][1]
    res$alpha_1.bayes <-   summary(posterior.tox, burn=5000, order=F)$coefficients[rownames(summary(posterior.tox, burn=5000, order=F)$coefficients)=='Dose.c'][1]
    res$alpha_2.bayes <- res$alpha_3.bayes <- 0
    if(summary(posterior.tox, burn=5000, order=F)$coefficients[rownames(summary(posterior.tox, burn=5000, order=F)$coefficients)=='Subgroup.c'][5]>inclusion.threshold){res$alpha_2.bayes <- summary(posterior.tox, burn=5000, order=F)$coefficients[rownames(summary(posterior.tox, burn=5000, order=F)$coefficients)=='Subgroup.c'][1]}
    if(summary(posterior.tox, burn=5000, order=F)$coefficients[rownames(summary(posterior.tox, burn=5000, order=F)$coefficients)=='Dose.c:Subgroup.c'][5]>inclusion.threshold){res$alpha_3.bayes <- summary(posterior.tox, burn=5000, order=F)$coefficients[rownames(summary(posterior.tox, burn=5000, order=F)$coefficients)=='Dose.c:Subgroup.c'][1]}
    
    if (is.null(toxicity_target)){
      res$beta_0.bayes <-   summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='(Intercept)'][1]
      res$beta_1.bayes <-   summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='Dose.c'][1]
      res$beta_2.bayes <-   summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='I(Dose.c^2)'][1]
      res$beta_3.bayes <- res$beta_4.bayes <- res$beta_5.bayes <- 0
      if(summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='Subgroup.c'][5]>inclusion.threshold){res$beta_3.bayes <- summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='Subgroup.c'][1]}
      if(summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='Dose.c:Subgroup.c'][5]>inclusion.threshold){res$beta_4.bayes <- summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='Dose.c:Subgroup.c'][1]}
      if(summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='I(Dose.c^2):Subgroup.c'][5]>inclusion.threshold){res$beta_5.bayes <- summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='I(Dose.c^2):Subgroup.c'][1]}
    }
    return(res)
  }
}

## model-based dose finding for one group
dose.find.onegroup <- function(FullDat, Design.eff, Design.eff.c, PriorSpec.tox, toxicity_target, PriorSpec.eff, tox.threshold, eff.threshold, doses, res, tox.level, eff.level, u00, u11){
  posterior.tox <- logit.spike(DLT~Dose.c, 
                               prior = PriorSpec.tox, data=FullDat, niter = 20000, na.action = na.exclude, ping = 0)
  pr.tox<-predict.logit.spike(posterior.tox, newdata = as.matrix(Design.eff.c)[,c(1,2)], 
                              burn=5000)
  admiss.tox <- rowMeans(pr.tox<tox.threshold)>tox.level
  #pr.tox <- rowMeans(pr.tox)
  tox.map <- rep(NA, nrow(pr.tox))
  for(i in 1:nrow(pr.tox)){
    tox.map[i] <- map_estimate(pr.tox[i,])
  }
  pr.tox <- tox.map
  
  if (is.null(toxicity_target)){
    posterior.eff <- logit.spike(eff~ Dose.c + I(Dose.c^2), 
                                 prior = PriorSpec.eff, data=FullDat, niter = 20000, na.action = na.exclude, ping = 0)
    pr.eff<-predict.logit.spike(posterior.eff, newdata = as.matrix(Design.eff.c), burn = 5000)
    admiss.eff <- rowMeans(pr.eff>eff.threshold)>eff.level
    #pr.eff <- rowMeans(pr.eff)
    eff.map <- rep(NA, nrow(pr.eff))
    for(i in 1:nrow(pr.eff)){
      eff.map[i] <- map_estimate(pr.eff[i,])
    }
    pr.eff <- eff.map
    utility <- utility_group(u00=u00, u11=u11, ptox=pr.tox, peff=pr.eff)
    
    if(!is.null(res)){
      if(all(FullDat$Subgroup == 0)){
        dose_options <-
          which(admiss.tox==T & admiss.eff==T & res[3:(2+length(doses))] >0)
      }else{
        dose_options <-
          which(admiss.tox==T & admiss.eff==T & res[(2+length(doses)+1):(2+length(doses)*2)] >0)
      }
      res$obd0.bayes <- res$obd1.bayes <- ifelse(length(dose_options) == 0, 0, doses[dose_options][which.max(utility[dose_options])])
    }else{
      # A dose must meet the toxicity threshold and EITHER meet the efficacy threshold or not have been tried yet
      dose_options <- which(admiss.tox==T & 
                              (admiss.eff==T | !(doses %in% FullDat$Dose)))
      if(length(dose_options) != 0){
        dose <- doses[dose_options][which.max(utility[dose_options])]
      }else{
        dose <- NA
        print("no dose available")
        print(pr.tox)
        print(pr.eff)
      }
    }
  }
  
  #### If you want to target a certain toxicity rate (like Cotteril and Jaki)
  if (!is.null(toxicity_target)){
    if(!is.null(res)){
      if(all(FullDat$Subgroup == 0)){
        dose_options <- which(admiss.tox==T & res[3:(2+length(doses))] >0)
      }else if(all(FullDat$Subgroup == 1)){
        dose_options <- which(admiss.tox==T & res[(2+length(doses)+1):(2+length(doses)*2)] >0)
      }
      
      res$ttd0.bayes <- res$ttd1.bayes <- ifelse(length(dose_options)==0, 0, doses[dose_options][which.min(abs(pr.tox[dose_options]-toxicity_target))])
    }else{
      dose_options <- which(admiss.tox==T)
      if(length(dose_options) != 0){
        dose <- doses[dose_options][which.min(abs(pr.tox[dose_options]-toxicity_target))]
      }else{
        dose <- NA
        print("no dose available")
        print(pr.tox)
      }
    }
  }
  
  if(is.null(res)){
    if(length(dose_options) != 0 & length(FullDat$Dose) > 0){
      # Don't allow escalation by more than 1 step above tried doses 
      if (max(FullDat$Dose, na.rm=TRUE) != doses[length(doses)] & 
          dose > doses[which(doses == max(FullDat$Dose, na.rm=TRUE))+1]){ 
        dose <- doses[which(doses == max(FullDat$Dose, na.rm=TRUE))+1]
      }
    }
    dose0 <- dose1 <- dose
    return(c(dose0, dose1))
  }else{
    res$alpha_0.bayes <-   summary(posterior.tox, burn=5000, order=F)$coefficients[rownames(summary(posterior.tox, burn=5000, order=F)$coefficients)=='(Intercept)'][1]
    res$alpha_1.bayes <-   summary(posterior.tox, burn=5000, order=F)$coefficients[rownames(summary(posterior.tox, burn=5000, order=F)$coefficients)=='Dose.c'][1]
    res$alpha_2.bayes <- res$alpha_3.bayes <- 0
    
    if (is.null(toxicity_target)){
      res$beta_0.bayes <-   summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='(Intercept)'][1]
      res$beta_1.bayes <-   summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='Dose.c'][1]
      res$beta_2.bayes <-   summary(posterior.eff, burn=5000, order=F)$coefficients[rownames(summary(posterior.eff, burn=5000, order=F)$coefficients)=='I(Dose.c^2)'][1]
      res$beta_3.bayes <- res$beta_4.bayes <- res$beta_5.bayes <- 0
    }
    return(res)
  }
}
