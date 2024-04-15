##################################
# helper functions to calculate pr.tox, pr.eff and utility and generate data
##################################
prob_logis<-function(lc){ 
  exp(lc)/(1+exp(lc))
}

# Calculate the utility given probability of efficacy and toxicity for a dose or 
# several doses (if peff, ptox are vectors)
utility_group <- function(u10=100,u11=60, u00=40,u01=0, ptox, peff){
  ut <- u10*peff*(1-ptox) + u11*peff*ptox + u00*(1-peff)*(1-ptox) + u01*(1-peff)*(ptox)
  return(ut)
}

data_gen <- function(FullDat, starting_data, alpha_0, alpha_1, alpha_2, alpha_3, beta_0, beta_1, beta_2, beta_3, beta_4, beta_5, dose0, dose1, correlation, p.tox.true, p.eff.true, doses){
  # If we haven't used all the rows of starting_data, use them
  if (nrow(FullDat) < nrow(starting_data)) {
    FullDat_add <-  starting_data[nrow(FullDat)+1:2,] 
    if(sum(is.na(FullDat_add[,2])) > 0){
      for(q in 1:2){
        if(FullDat_add[q,1] == 0){
          FullDat_add[q,2] <- dose0
          if(is.null(p.tox.true) | is.null(p.eff.true)){
            if(is.null(correlation)){
              if(!is.na(dose0)){FullDat_add[q,3] <- rbinom(1,1,prob_logis(alpha_0+alpha_1*FullDat_add[q,2]))}
              if(!is.na(dose0)){FullDat_add[q,4] <- rbinom(1,1,prob_logis(beta_0+beta_1*FullDat_add[q,2]+beta_2*FullDat_add[q,2]^2))}
            }else{
              if(!is.na(dose0)){
                margprob0 <- c(prob_logis(alpha_0+alpha_1*FullDat_add[q,2]), prob_logis(beta_0+beta_1*FullDat_add[q,2]+beta_2*FullDat_add[q,2]^2))
                
                p00_0 <- ifelse((1-margprob0[1])*(1-margprob0[2]) + correlation*sqrt(margprob0[1]*margprob0[2]*(1-margprob0[1])*(1-margprob0[2]))>0, (1-margprob0[1])*(1-margprob0[2]) + correlation*sqrt(margprob0[1]*margprob0[2]*(1-margprob0[1])*(1-margprob0[2])), 0)
                p10_0 <- ifelse(1-margprob0[1]-p00_0>0, 1-margprob0[1]-p00_0, 0)
                p01_0 <- ifelse(1-margprob0[2]-p00_0>0, 1-margprob0[2]-p00_0, 0)
                p11_0 <- ifelse(p00_0 + margprob0[1] + margprob0[2] -1>0, p00_0 + margprob0[1] + margprob0[2] -1, 0)
                
                outcome <- which(rmultinom(1,1,prob=c(p00_0, p10_0, p01_0, p11_0)) == 1)
                
                FullDat_add[q,3] <- ifelse(outcome == 1 | outcome == 2, 0, 1)
                FullDat_add[q,4] <- ifelse(outcome == 1 | outcome == 3, 0, 1)
              }
            }
          }else{
            if(!is.na(dose0)){FullDat_add[q,3] <- rbinom(1,1,p.tox.true[1:length(doses)][which(doses == FullDat_add[q,2])])}
            if(!is.na(dose0)){FullDat_add[q,4] <- rbinom(1,1,p.eff.true[1:length(doses)][which(doses == FullDat_add[q,2])])}
          }
          
        }else if(FullDat_add[q,1] == 1){
          FullDat_add[q,2] <- dose1
          if(is.null(p.tox.true) | is.null(p.eff.true)){
            if(is.null(correlation)){
              if(!is.na(dose1)){FullDat_add[q,3] <- rbinom(1,1,prob_logis(alpha_0+alpha_1*FullDat_add[q,2]+alpha_2+alpha_3*FullDat_add[q,2]))}
              if(!is.na(dose1)){FullDat_add[q,4] <- rbinom(1,1,prob_logis(beta_0+beta_1*FullDat_add[q,2]+beta_2*FullDat_add[q,2]^2 + beta_3 + beta_4*FullDat_add[q,2] +beta_5*FullDat_add[q,2]^2))}
            }else{
              if(!is.na(dose1)){
                margprob1 <- c(prob_logis(alpha_0+alpha_1*FullDat_add[q,2]+alpha_2+alpha_3*FullDat_add[q,2]), prob_logis(beta_0+beta_1*FullDat_add[q,2]+beta_2*FullDat_add[q,2]^2 + beta_3 + beta_4*FullDat_add[q,2] +beta_5*FullDat_add[q,2]^2))
                
                p00_1 <- ifelse((1-margprob1[1])*(1-margprob1[2]) + correlation*sqrt(margprob1[1]*margprob1[2]*(1-margprob1[1])*(1-margprob1[2]))>0, (1-margprob1[1])*(1-margprob1[2]) + correlation*sqrt(margprob1[1]*margprob1[2]*(1-margprob1[1])*(1-margprob1[2])), 0)
                p10_1 <- ifelse(1-margprob1[1]-p00_1>0, 1-margprob1[1]-p00_1, 0)
                p01_1 <- ifelse(1-margprob1[2]-p00_1>0, 1-margprob1[2]-p00_1, 0)
                p11_1 <- ifelse(p00_1 + margprob1[1] + margprob1[2] -1>0, p00_1 + margprob1[1] + margprob1[2] -1, 0)
                
                outcome <- which(rmultinom(1,1,prob=c(p00_1, p10_1, p01_1, p11_1)) == 1)
                
                FullDat_add[q,3] <- ifelse(outcome == 1 | outcome == 2, 0, 1)
                FullDat_add[q,4] <- ifelse(outcome == 1 | outcome == 3, 0, 1)
              }
            }
          }else{
            if(!is.na(dose1)){FullDat_add[q,3] <- rbinom(1,1,p.tox.true[(length(doses)+1):(2*length(doses))][which(doses == FullDat_add[q,2])])}
            if(!is.na(dose1)){FullDat_add[q,4] <- rbinom(1,1,p.eff.true[(length(doses)+1):(2*length(doses))][which(doses == FullDat_add[q,2])])}
          }
          
        }
      }
    }
  }else{ # Otherwise simulate data 
    FullDat_add<- data.frame(Subgroup = c(0,1), 
                             Dose = c(dose0,dose1),
                             DLT = rep(NA,2),
                             eff = rep(NA,2))
    if(is.null(p.tox.true) | is.null(p.eff.true)){
      if(is.null(correlation)){
        if(!is.na(dose0)){FullDat_add$DLT[1] <- rbinom(1,1,prob_logis(alpha_0+alpha_1*FullDat_add$Dose[1]))}
        if(!is.na(dose1)){FullDat_add$DLT[2] <- rbinom(1,1,prob_logis(alpha_0+alpha_1*FullDat_add$Dose[2]+alpha_2+alpha_3*FullDat_add$Dose[2]))}
        if(!is.na(dose0)){FullDat_add$eff[1] <- rbinom(1,1,prob_logis(beta_0+beta_1*FullDat_add$Dose[1]+beta_2*FullDat_add$Dose[1]^2))}
        if(!is.na(dose1)){FullDat_add$eff[2] <- rbinom(1,1,prob_logis(beta_0+beta_1*FullDat_add$Dose[2]+beta_2*FullDat_add$Dose[2]^2 + beta_3 + beta_4*FullDat_add$Dose[2] +beta_5*FullDat_add$Dose[2]^2))}
      }else{
        if(!is.na(dose0)){
          margprob0 <- c(prob_logis(alpha_0+alpha_1*FullDat_add$Dose[1]), prob_logis(beta_0+beta_1*FullDat_add$Dose[1]+beta_2*FullDat_add$Dose[1]^2))
          
          p00_0 <- ifelse((1-margprob0[1])*(1-margprob0[2]) + correlation*sqrt(margprob0[1]*margprob0[2]*(1-margprob0[1])*(1-margprob0[2]))>0, (1-margprob0[1])*(1-margprob0[2]) + correlation*sqrt(margprob0[1]*margprob0[2]*(1-margprob0[1])*(1-margprob0[2])), 0)
          p10_0 <- ifelse(1-margprob0[1]-p00_0>0, 1-margprob0[1]-p00_0, 0)
          p01_0 <- ifelse(1-margprob0[2]-p00_0>0, 1-margprob0[2]-p00_0, 0)
          p11_0 <- ifelse(p00_0 + margprob0[1] + margprob0[2] -1>0, p00_0 + margprob0[1] + margprob0[2] -1, 0)
          
          outcome <- which(rmultinom(1,1,prob=c(p00_0, p10_0, p01_0, p11_0)) == 1)
          
          FullDat_add$DLT[1] <- ifelse(outcome == 1 | outcome == 2, 0, 1)
          FullDat_add$eff[1] <- ifelse(outcome == 1 | outcome == 3, 0, 1)
        }
        if(!is.na(dose1)){
          margprob1 <- c(prob_logis(alpha_0+alpha_1*FullDat_add$Dose[2]+alpha_2+alpha_3*FullDat_add$Dose[2]), prob_logis(beta_0+beta_1*FullDat_add$Dose[2]+beta_2*FullDat_add$Dose[2]^2 + beta_3 + beta_4*FullDat_add$Dose[2] +beta_5*FullDat_add$Dose[2]^2))
          
          p00_1 <- ifelse((1-margprob1[1])*(1-margprob1[2]) + correlation*sqrt(margprob1[1]*margprob1[2]*(1-margprob1[1])*(1-margprob1[2]))>0, (1-margprob1[1])*(1-margprob1[2]) + correlation*sqrt(margprob1[1]*margprob1[2]*(1-margprob1[1])*(1-margprob1[2])), 0)
          p10_1 <- ifelse(1-margprob1[1]-p00_1>0, 1-margprob1[1]-p00_1, 0)
          p01_1 <- ifelse(1-margprob1[2]-p00_1>0, 1-margprob1[2]-p00_1, 0)
          p11_1 <- ifelse(p00_1 + margprob1[1] + margprob1[2] -1>0, p00_1 + margprob1[1] + margprob1[2] -1, 0)
          
          outcome <- which(rmultinom(1,1,prob=c(p00_1, p10_1, p01_1, p11_1)) == 1)
          
          FullDat_add$DLT[2] <- ifelse(outcome == 1 | outcome == 2, 0, 1)
          FullDat_add$eff[2] <- ifelse(outcome == 1 | outcome == 3, 0, 1)
        }
      }
    }else{
      if(!is.na(dose0)){FullDat_add$DLT[1] <- rbinom(1,1,p.tox.true[1:length(doses)][which(doses == FullDat_add$Dose[1])])}
      if(!is.na(dose1)){FullDat_add$DLT[2] <- rbinom(1,1,p.tox.true[(length(doses)+1):(2*length(doses))][which(doses == FullDat_add$Dose[2])])}
      if(!is.na(dose0)){FullDat_add$eff[1] <- rbinom(1,1,p.eff.true[1:length(doses)][which(doses == FullDat_add$Dose[1])])}
      if(!is.na(dose1)){FullDat_add$eff[2] <- rbinom(1,1,p.eff.true[(length(doses)+1):(2*length(doses))][which(doses == FullDat_add$Dose[2])])}
    }
    
  }
  return(FullDat_add)
}


