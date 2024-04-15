library(BoomSpikeSlab)
library(parallel)
library(bayestestR)

source("helperfuns_anyorder.R")
source("dosefind_prior_MAP.R")

simtrial<- function(ncohorts, doses=c(.2,.4,.6,.8,1), starting_data=data.frame(),
                    rule_based_runin = T, 
                    alpha_0 = -4, alpha_1=2.5, alpha_2=1.5, alpha_3 = .9,
                    beta_0 =-2.5, beta_1=8, beta_2= -4.5, beta_3=0, beta_4=0, beta_5=0, 
                    u10=100,u11=60, u00=40,u01=0, p.eff.true=NULL, p.tox.true=NULL, hyperparams,
                    PIP.tox =c(1,1,.5,.5),
                    prior.information.weight.tox = 1,
                    PIP.eff = c(1,1,1,.5,.5,.5),
                    prior.information.weight.eff = 1, toxicity_target = NULL, tox.threshold=.3,
                    eff.threshold=.35, inclusion.threshold=.3,
                    print.last.spikeslab=F, one.group = F, correlation=NULL,
                    tox.level=.1,  eff.level=.1){
  
  if (!is.null(toxicity_target)) {if (toxicity_target>tox.threshold) stop('Target toxicity must be lower than toxicity threshold.')}
  # Build the priors
  if(one.group == F){
    PriorX.tox <- diag(4) # Setting the covariance between the slabs = 0
    PriorX.tox[1,1] <- 4/hyperparams$b0t_sd
    PriorX.tox[2,2] <- 4/hyperparams$b1t_sd
    PriorX.tox[3,3] <- 4/hyperparams$b2t_sd
    PriorX.tox[4,4] <- 4/hyperparams$b3t_sd
    names(PriorX.tox) <- c("Intercept", "Dose", "Subgroup", "Subgroup*Dose")
    PriorSpec.tox <- LogitZellnerPrior( 
      predictors = as.matrix(PriorX.tox),
      successes = NULL,
      trials = NULL,
      prior.success.probability = .5,
      prior.information.weight = prior.information.weight.tox,
      diagonal.shrinkage = 0,
      optional.coefficient.estimate = c(hyperparams$b0t_mean, hyperparams$b1t_mean, hyperparams$b2t_mean, hyperparams$b3t_mean),
      max.flips = -1,
      prior.inclusion.probabilities = PIP.tox)
    
    PriorX.eff <- diag(6)
    PriorX.eff[1,1] <- sqrt(24)/hyperparams$b0e_sd
    PriorX.eff[2,2] <- sqrt(24)/hyperparams$b1e_sd
    PriorX.eff[3,3] <- sqrt(24)/hyperparams$b2e_sd
    PriorX.eff[4,4] <- sqrt(24)/hyperparams$b3e_sd
    PriorX.eff[5,5] <- sqrt(24)/hyperparams$b4e_sd
    PriorX.eff[6,6] <- sqrt(24)/hyperparams$b5e_sd
    names(PriorX.eff) <- c("Intercept", "Dose", "Dose^2", "Subgroup", "Subgroup*Dose", "Subgroup*Dose^2")
    PriorSpec.eff <- LogitZellnerPrior( 
      predictors = as.matrix(PriorX.eff),
      successes = NULL,
      trials = NULL,
      prior.success.probability = .5,
      prior.information.weight = prior.information.weight.eff,
      diagonal.shrinkage = 0,
      optional.coefficient.estimate = c(hyperparams$b0e_mean, hyperparams$b1e_mean, hyperparams$b2e_mean, hyperparams$b3e_mean, hyperparams$b4e_mean, hyperparams$b5e_mean),
      max.flips = -1,
      prior.inclusion.probabilities = PIP.eff)
    
    # All possible covariate combinations
    Design.eff <- data.frame(base::cbind(rep(1,length(doses)*2), rep(doses, each =2),
                                   rep(doses^2, each =2),
                                   rep(c(0,1), length.out = length(doses)*2),
                                   rep(doses, each =2)*rep(c(0,1), length.out=length(doses)*2),
                                   rep(doses^2, each =2)*rep(c(0,1), length.out=length(doses)*2)))
    names(Design.eff) <- c("Intercept", "Dose", "Dose^2", "Subgroup", "Subgroup*Dose", "Subgroup*Dose^2")
    
    Design.eff.c <- data.frame(base::cbind(rep(1,length(doses)*2), rep(log(doses) - (1/length(doses))*sum(log(doses)), each =2),
                                   rep((log(doses) - (1/length(doses))*sum(log(doses)))^2, each =2),
                                   rep(c(-1,1), length.out = length(doses)*2),
                                   rep(log(doses) - (1/length(doses))*sum(log(doses)), each =2)*rep(c(-1,1), length.out=length(doses)*2),
                                   rep((log(doses) - (1/length(doses))*sum(log(doses)))^2, each =2)*rep(c(-1,1), length.out=length(doses)*2)))
    names(Design.eff.c) <- c("Intercept", "Dose", "Dose^2", "Subgroup", "Subgroup*Dose", "Subgroup*Dose^2")
    
    #Find the true OBDs
    if(is.null(p.eff.true) | is.null(p.tox.true)){
      if(is.null(correlation)){
        pr.tox <- c(rbind(prob_logis(alpha_0 + alpha_1*doses),
                          prob_logis(alpha_0 + alpha_1*doses+ alpha_2 + alpha_3*doses)))
        pr.eff <- c(rbind(prob_logis(beta_0 + beta_1*doses + beta_2*doses^2),
                          prob_logis(beta_0 + beta_1*doses + beta_2*doses^2+
                                       beta_3+ beta_4*doses  + beta_5*doses^2)))
        utility <- utility_group(u11=u11, u00=u00,ptox=pr.tox, peff=pr.eff)
        
      }else{
        util0 <- rep(NA, length(doses))
        util1 <- rep(NA, length(doses))
        for(m in 1:length(doses)){
          margprob0 <- c(prob_logis(alpha_0+alpha_1*doses[m]), prob_logis(beta_0+beta_1*doses[m]+beta_2*doses[m]^2))
          p00_0 <- ifelse((1-margprob0[1])*(1-margprob0[2]) + correlation*sqrt(margprob0[1]*margprob0[2]*(1-margprob0[1])*(1-margprob0[2]))>0, (1-margprob0[1])*(1-margprob0[2]) + correlation*sqrt(margprob0[1]*margprob0[2]*(1-margprob0[1])*(1-margprob0[2])), 0)
          p10_0 <- ifelse(1-margprob0[1]-p00_0>0, 1-margprob0[1]-p00_0, 0)
          p01_0 <- ifelse(1-margprob0[2]-p00_0>0, 1-margprob0[2]-p00_0, 0)
          p11_0 <- ifelse(p00_0 + margprob0[1] + margprob0[2] -1>0, p00_0 + margprob0[1] + margprob0[2] -1, 0)
          util0[m] <- p00_0*u00 + p10_0*u10 + p01_0*u01 + p11_0*u11
          
          margprob1 <- c(prob_logis(alpha_0+alpha_1*doses[m]+alpha_2+alpha_3*doses[m]), prob_logis(beta_0+beta_1*doses[m]+beta_2*doses[m]^2 + beta_3 + beta_4*doses[m] +beta_5*doses[m]^2))
          p00_1 <- ifelse((1-margprob1[1])*(1-margprob1[2]) + correlation*sqrt(margprob1[1]*margprob1[2]*(1-margprob1[1])*(1-margprob1[2]))>0, (1-margprob1[1])*(1-margprob1[2]) + correlation*sqrt(margprob1[1]*margprob1[2]*(1-margprob1[1])*(1-margprob1[2])), 0)
          p10_1 <- ifelse(1-margprob1[1]-p00_1>0, 1-margprob1[1]-p00_1, 0)
          p01_1 <- ifelse(1-margprob1[2]-p00_1>0, 1-margprob1[2]-p00_1, 0)
          p11_1 <- ifelse(p00_1 + margprob1[1] + margprob1[2] -1>0, p00_1 + margprob1[1] + margprob1[2] -1, 0)
          util1[m] <- p00_1*u00 + p10_1*u10 + p01_1*u01 + p11_1*u11
        }
        utility <- c(rbind(util0, util1))
      }
    }else{
      pr.tox <- c(rbind(p.tox.true[1:length(doses)],
                        p.tox.true[(length(doses)+1):(2*length(doses))]))
      pr.eff <- c(rbind(p.eff.true[1:length(doses)],
                        p.eff.true[(length(doses)+1):(2*length(doses))]))
      utility <- utility_group(u11=u11, u00=u00,ptox=pr.tox, peff=pr.eff)
    }
    
    trueOBD0 <- doses[which(utility[Design.eff$Subgroup==0] ==
                              max(utility[Design.eff$Subgroup==0]))]
    trueOBD1 <- doses[which(utility[Design.eff$Subgroup==1] ==
                              max(utility[Design.eff$Subgroup==1]))]
    
  }else if(one.group == T){
    PriorX.tox <- diag(2)
    PriorX.tox[1,1] <- sqrt(8)/hyperparams$b0t_sd
    PriorX.tox[2,2] <- sqrt(8)/hyperparams$b1t_sd
    names(PriorX.tox) <- c("Intercept", "Dose")
    PriorSpec.tox <- LogitZellnerPrior( 
      predictors = as.matrix(PriorX.tox),
      successes = NULL,
      trials = NULL,
      prior.success.probability = .5,
      prior.information.weight = prior.information.weight.tox,
      diagonal.shrinkage = 0,
      optional.coefficient.estimate = c(hyperparams$b0t_mean, hyperparams$b1t_mean),
      max.flips = -1,
      prior.inclusion.probabilities = PIP.tox[1:2])
    
    PriorX.eff <- diag(3)
    PriorX.eff[1,1] <- sqrt(12)/hyperparams$b0e_sd
    PriorX.eff[2,2] <- sqrt(12)/hyperparams$b1e_sd
    PriorX.eff[3,3] <- sqrt(12)/hyperparams$b2e_sd
    names(PriorX.eff) <- c("Intercept", "Dose", "Dose^2")
    PriorSpec.eff <- LogitZellnerPrior( 
      predictors = as.matrix(PriorX.eff),
      successes = NULL,
      trials = NULL,
      prior.success.probability = .5,
      prior.information.weight = prior.information.weight.eff,
      diagonal.shrinkage = 0,
      optional.coefficient.estimate = c(hyperparams$b0e_mean, hyperparams$b1e_mean, hyperparams$b2e_mean),
      max.flips = -1,
      prior.inclusion.probabilities = PIP.eff[1:3])
    
    # Design matrix
    Design.eff <- data.frame(base::cbind(rep(1,length(doses)), doses,
                                   doses^2))
    names(Design.eff) <- c("Intercept", "Dose", "Dose^2")
    
    Design.eff.c <- data.frame(cbind(rep(1,length(doses)), log(doses) - (1/length(doses))*sum(log(doses)),
                                   (log(doses) - (1/length(doses))*sum(log(doses)))^2))
    names(Design.eff.c) <- c("Intercept", "Dose", "Dose^2")
    
    if(is.null(p.tox.true) | is.null(p.eff.true)){
      if(is.null(correlation)){
        pr.tox <- prob_logis(alpha_0 + alpha_1*doses)
        pr.eff <- prob_logis(beta_0 + beta_1*doses + beta_2*doses^2)
        utility <- utility_group(u11=u11, u00=u00,ptox=pr.tox, peff=pr.eff)
        
      }else{
        utility <- rep(NA, length(doses))
        for(m in 1:length(doses)){
          margprob <- c(prob_logis(alpha_0+alpha_1*doses[m]), prob_logis(beta_0+beta_1*doses[m]+beta_2*doses[m]^2))
          p00 <- ifelse((1-margprob[1])*(1-margprob[2]) + correlation*sqrt(margprob[1]*margprob[2]*(1-margprob[1])*(1-margprob[2]))>0, (1-margprob[1])*(1-margprob[2]) + correlation*sqrt(margprob[1]*margprob[2]*(1-margprob[1])*(1-margprob[2])), 0)
          p10 <- ifelse(1-margprob[1]-p00>0, 1-margprob[1]-p00, 0)
          p01 <- ifelse(1-margprob[2]-p00>0, 1-margprob[2]-p00, 0)
          p11 <- ifelse(p00 + margprob[1] + margprob[2] -1>0, p00 + margprob[1] + margprob[2] -1, 0)
          utility[m] <- p00*u00 + p10*u10 + p01*u01 + p11*u11
          
        }
      }
    }else{
      pr.tox <- p.tox.true
      pr.eff <- p.tox.true
      utility <- utility_group(u11=u11, u00=u00,ptox=pr.tox, peff=pr.eff)
    }
    
    trueOBD <- doses[which(utility ==
                             max(utility))]
    trueOBD0 <- trueOBD
    trueOBD1 <- trueOBD
    
  }
  
  # Make a data frame to store the results
  if (is.null(toxicity_target)){ res <- data.frame(obd0.bayes=NA, obd1.bayes=NA)
  }else {res <- data.frame(ttd0.bayes=NA, ttd1.bayes=NA)}
  for(i in 1:length(doses)){
    res <- cbind(res, NA)
    colnames(res)[i+2] <- paste0("N",i,"_0",collapse="")
  }
  for(i in 1:length(doses)){
    res <- cbind(res, NA)
    colnames(res)[i+2+length(doses)] <- paste0("N",i,"_1",collapse="")
  }
  for(i in 1:length(doses)){
    res <- cbind(res, NA)
    colnames(res)[i+2+length(doses)*2] <- paste0("d",i,"_0",collapse="")
  }
  for(i in 1:length(doses)){
    res <- cbind(res, NA)
    colnames(res)[i+2+length(doses)*3] <- paste0("d",i,"_1",collapse="")
  }
  for(i in 1:length(doses)){
    res <- cbind(res, NA)
    colnames(res)[i+2+length(doses)*4] <- paste0("e",i,"_0",collapse="")
  }
  for(i in 1:length(doses)){
    res <- cbind(res, NA)
    colnames(res)[i+2+length(doses)*5] <- paste0("e",i,"_1",collapse="")
  }
  res <- cbind(res, alpha_0.bayes = NA, alpha_1.bayes=NA, alpha_2.bayes=NA, alpha_3.bayes = NA, # Estimated coefficients
               beta_0.bayes = NA, beta_1.bayes=NA, beta_2.bayes= NA, beta_3.bayes=NA, beta_4.bayes=NA, beta_5.bayes=NA,
               pip.alpha2=NA, pip.alpha3=NA, pip.beta3=NA, pip.beta4=NA, pip.beta5=NA, pip.subgroup.tox=NA, pip.subgroup.eff=NA)
  
  dose1<- dose0 <- doses[1]
  
  # Simulate a single trial 
  FullDat <- data.frame()
  for(i in seq(from=1, to=ncohorts)){
    FullDat_add <- data_gen(FullDat=FullDat, starting_data=starting_data, alpha_0=alpha_0, alpha_1=alpha_1, alpha_2=alpha_2,
                            alpha_3=alpha_3, beta_0=beta_0, beta_1=beta_1, beta_2=beta_2, beta_3=beta_3, beta_4=beta_4, beta_5=beta_5, dose0=dose0, dose1=dose1, correlation=correlation, p.tox.true=p.tox.true, p.eff.true=p.eff.true, doses=doses)
    FullDat_add$Subgroup.c <- ifelse(FullDat_add$Subgroup == 0, -1, FullDat_add$Subgroup)
    FullDat_add$Dose.c <- log(FullDat_add$Dose) - (1/length(doses))*sum(log(doses))
    FullDat <- rbind(FullDat, FullDat_add)
    
    # Determine the dose for the next cohort
    if (rule_based_runin == T & length(unique(FullDat$Dose))==1 & nrow(FullDat)<=12){ # Rule-based run-in
      dosing <- runin(doses=doses, FullDat=FullDat, one.group=one.group)
    }else {
      if(one.group == F){
        dosing <- dose.find.twogroup(FullDat=FullDat, Design.eff=Design.eff, Design.eff.c=Design.eff.c, PriorSpec.tox=PriorSpec.tox, toxicity_target=toxicity_target, PriorSpec.eff=PriorSpec.eff, tox.threshold=tox.threshold, eff.threshold=eff.threshold, doses=doses, res=NULL, inclusion.threshold=inclusion.threshold, tox.level=tox.level, eff.level=eff.level, u00=u00, u11=u11)
      }else if(one.group == T){
        dosing <- dose.find.onegroup(FullDat=FullDat, Design.eff=Design.eff, Design.eff.c=Design.eff.c, PriorSpec.tox=PriorSpec.tox, toxicity_target=toxicity_target, PriorSpec.eff=PriorSpec.eff, tox.threshold=tox.threshold, eff.threshold=eff.threshold, doses=doses, res=NULL, tox.level=tox.level, eff.level=eff.level, u00=u00, u11=u11)
      }
    }
    # If neither group has any safe and effective doses, the trial ends
    if(is.na(dosing[1]) & is.na(dosing[2])){break}
    if(is.na(dosing[1]) & table(FullDat$Subgroup)[2] >= ncohorts){break}
    if(is.na(dosing[2]) & table(FullDat$Subgroup)[1] >= ncohorts){break}
    
    dose0 <- dosing[1]
    dose1 <- dosing[2]
    
  }
  
  #### Summarize the trial and get the final OBDs   
  if (print.last.spikeslab == T) {
    print(summary(posterior.tox, burn = 5000, order = F))
    if (is.null(toxicity_target)){print(summary(posterior.eff, burn = 5000, order = F))}
  }
  FullDat$Dose.factor <- factor(FullDat$Dose, levels = doses) 
  FullDat$Subgroup.factor <- factor(FullDat$Subgroup, levels=0:1)
  # Make dose a factor so that the tables contain all doses, even if not all were used
  # Count number of people, DLTs, efficacy events per subgroup/dose
  res[,3:(length(doses)*2+2)] <- c(table(FullDat$Dose.factor, FullDat$Subgroup.factor))
  if (sum(FullDat$DLT, na.rm=TRUE) == 0) {res[,(length(doses)*2+2+1):(length(doses)*4+2)] <- 0
  }else {res[,(length(doses)*2+2+1):(length(doses)*4+2)] <- c(table(FullDat$Dose.factor[FullDat$DLT==1], FullDat$Subgroup.factor[FullDat$DLT==1]))}
  if (sum(FullDat$eff, na.rm=TRUE) == 0) {res[,(length(doses)*4+2+1):(length(doses)*6+2)] <- 0
  }else {res[,(length(doses)*4+2+1):(length(doses)*6+2)] <- c(table(FullDat$Dose.factor[FullDat$eff==1], FullDat$Subgroup.factor[FullDat$eff==1]))}

  # If we never got past rule-based run-in, don't estimate final OBDs
  if  (rule_based_runin == T & length(unique(FullDat$Dose))==1 & nrow(FullDat)<=12) return(list(res=res, FullDat = FullDat))
  
  if(one.group == F){
    res <- dose.find.twogroup(FullDat=FullDat, Design.eff=Design.eff, Design.eff.c=Design.eff.c, PriorSpec.tox=PriorSpec.tox, toxicity_target=toxicity_target, PriorSpec.eff=PriorSpec.eff, tox.threshold=tox.threshold, eff.threshold=eff.threshold, doses=doses, res=res, inclusion.threshold=inclusion.threshold, tox.level=tox.level, eff.level=eff.level, u00=u00, u11=u11)
  }else if(one.group == T){
    res <- dose.find.onegroup(FullDat=FullDat, Design.eff=Design.eff, Design.eff.c=Design.eff.c, PriorSpec.tox=PriorSpec.tox, toxicity_target=toxicity_target, PriorSpec.eff=PriorSpec.eff, tox.threshold=tox.threshold, eff.threshold=eff.threshold, doses=doses, res=res, tox.level=tox.level, eff.level=eff.level, u00=u00, u11=u11)
  }
  return(list(res=res, FullDat = FullDat))    
}
