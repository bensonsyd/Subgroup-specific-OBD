library(gtools)

twogroup.pp <- function(doses, hyperparams, sd, n, seed, u01, u00, u11, u10){
  set.seed(seed)
  
  pr <- list(
    b0t_mean = hyperparams$b0t_mean, b0t_sd = hyperparams$b0t_sd,
    b1t_mean = hyperparams$b1t_mean, b1t_sd = hyperparams$b1t_sd,
    b2t_mean = 0, b2t_sd = sd,
    b3t_mean = 0, b3t_sd = sd,
    b0e_mean = hyperparams$b0e_mean, b0e_sd = hyperparams$b0e_sd,
    b1e_mean = hyperparams$b1e_mean, b1e_sd = hyperparams$b1e_sd,
    b2e_mean = hyperparams$b2e_mean, b2e_sd = hyperparams$b2e_sd,
    b3e_mean = 0, b3e_sd = sd,
    b4e_mean = 0, b4e_sd = sd,
    b5e_mean = 0, b5e_sd = .2
  )
  
  utility0 <- matrix(rep(NA, n*length(doses)), ncol=length(doses))
  utility1 <- matrix(rep(NA, n*length(doses)), ncol=length(doses))
  
  for(i in 1:length(doses)){
    
    t0 <- rnorm(n, pr$b0t_mean, sd=pr$b0t_sd)
    t1 <- rnorm(n, pr$b1t_mean, sd=pr$b1t_sd)
    t2 <- rnorm(n, pr$b2t_mean, pr$b2t_sd)
    t3 <- rnorm(n, pr$b3t_mean, pr$b3t_sd)
    
    logitt0 <- t0 + t1*doses[i] - t2 - t3*doses[i]
    logitt1 <- t0 + t1*doses[i] + t2 + t3*doses[i]
    
    e0 <- rnorm(n, pr$b0e_mean, sd=pr$b0e_sd)
    e1 <- rnorm(n, pr$b1e_mean, sd=pr$b1e_sd)
    e2 <- rnorm(n, pr$b2e_mean, sd=pr$b2e_sd)
    e3 <- rnorm(n, pr$b3e_mean, pr$b3e_sd)
    e4 <- rnorm(n, pr$b4e_mean, pr$b4e_sd)
    e5 <- rnorm(n, pr$b5e_mean, pr$b5e_sd)
    
    logite0 <- e0 + e1*doses[i] + e2*doses[i]^2 - e3 - e4*doses[i] - e5*doses[i]^2
    logite1 <- e0 + e1*doses[i] + e2*doses[i]^2 + e3 + e4*doses[i] + e5*doses[i]^2
    
    pt0 <- exp(logitt0)/(1+exp(logitt0))
    pt1 <- exp(logitt1)/(1+exp(logitt1))
    pe0 <- exp(logite0)/(1+exp(logite0))
    pe1 <- exp(logite1)/(1+exp(logite1))
    
    utility0[,i] <- u01*(1-pe0)*pt0 + u11*pe0*pt0 + u00*(1-pe0)*(1-pt0) + u10*pe0*(1-pt0)
    utility1[,i] <- u01*(1-pe1)*pt1 + u11*pe1*pt1 + u00*(1-pe1)*(1-pt1) + u10*pe1*(1-pt1)
    
  }
  
  O0 <- table(factor(apply(utility0, 1, which.max), levels=1:length(doses)))/n
  O1 <- table(factor(apply(utility1, 1, which.max), levels=1:length(doses)))/n
  
  prior.predictive <- list(subgroup1=O0, subgroup2=O1)
  res <- list(hyperparams=pr, prior.predictive=prior.predictive)
  return(res)
  
}

onegroup.pp <- function(doses, hyperparams, n, seed, u01, u00, u11, u10){
  set.seed(seed)
  
  pr <- list(
    b0t_mean = hyperparams$alpha_mean, b0t_sd = hyperparams$alpha_sd,
    b1t_mean = hyperparams$beta_mean, b1t_sd = hyperparams$beta_sd,
    b0e_mean = hyperparams$gamma_mean, b0e_sd = hyperparams$gamma_sd,
    b1e_mean = hyperparams$zeta_mean, b1e_sd = hyperparams$zeta_sd,
    b2e_mean = hyperparams$eta_mean, b2e_sd = hyperparams$eta_sd
  )
  
  utility <- matrix(rep(NA, n*length(doses)), ncol=length(doses))
  
  for(i in 1:length(doses)){
    
    t0 <- rnorm(n, pr$b0t_mean, sd=pr$b0t_sd)
    t1 <- rnorm(n, pr$b1t_mean, sd=pr$b1t_sd)
    
    logitt <- t0 + t1*doses[i]
    
    e0 <- rnorm(n, pr$b0e_mean, sd=pr$b0e_sd)
    e1 <- rnorm(n, pr$b1e_mean, sd=pr$b1e_sd)
    e2 <- rnorm(n, pr$b2e_mean, sd=pr$b2e_sd)
    
    logite <- e0 + e1*doses[i] + e2*doses[i]^2
    
    pt <- exp(logitt)/(1+exp(logitt))
    pe <- exp(logite)/(1+exp(logite))
    
    utility[,i] <- u01*(1-pe)*pt + u11*pe*pt + u00*(1-pe)*(1-pt) + u10*pe*(1-pt)
    
  }
  
  O <- table(factor(apply(utility, 1, which.max), levels=1:length(doses)))/n
  
  prior.predictive <- O
  res <- list(hyperparams=pr, prior.predictive=prior.predictive)
  return(res)
  
}



