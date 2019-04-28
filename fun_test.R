library(survival)

pt <- function(gamma, X=NULL){
  # Calculate the probability based on logistic model
  # Input gamma: coef in logistic model
  
  if(length(gamma) == 1){
    p_t <- exp(gamma)
  }else{
    n <- dim(X)[1]
    XX <- cbind(rep(1,n), X)
    p_t <- exp(XX %*% gamma)
  }
  p_t <- p_t /(1 + p_t)
  return(p_t)
}


score_test1 <- function(gamma){
  # Calculate the proposed score test statistic given gamma
  # Input gamma: coef in logistic model
  # Output score: The score value in the numerator of proposed score statistic
  # Output Vs: The variance of the averaged score; n*Vs is the denominator 
  #             of the proposed score test statistic
  # Output stat: The value of the proposed statistic given gamma
  # Output phistar: The phistar in the denominator of the proposed statistic given gamma
  p <- dim(X)[2]
  rk_tildeT <- rank(Y[, 1])
  order_tildeT <- order(Y[, 1])
  sort_delta <- Y[, 2]
  sort_delta[rk_tildeT] <- Y[, 2]
  fit1 <- coxph(Y ~ X)
  beta_hat <- summary(fit1)$coef[, 1]
  
  pi <- drop(pt(gamma, X))
  eXbeta <- drop(exp(X %*% beta_hat))
  X_sort <- X
  ztilde <- pi * (NET %*% X)
  zstar <- drop(ztilde %*% beta_hat)
  
  
  eXbeta[rk_tildeT] <- eXbeta
  zstar[rk_tildeT] <- zstar
  X_sort[rk_tildeT, ] <- X_sort
  ztilde[rk_tildeT, ] <- ztilde
  
  eXbeta_sum <- sum(eXbeta) - cumsum(eXbeta) + eXbeta
  zstar_eXbeta_sum <- sum(zstar * eXbeta) - cumsum(zstar * eXbeta) + zstar * eXbeta
  
  score <- sum(sort_delta * (zstar - zstar_eXbeta_sum / eXbeta_sum))
  eXbeta_X <- eXbeta * X_sort 
  num_beta2 <- t(colSums(eXbeta_X) - t(apply(eXbeta_X, 2, cumsum))) + eXbeta_X
  num_beta2_2 <- matrix(t(apply(num_beta2, 1, tcrossprod)), ncol=p^2)
  
  num_beta1 <- matrix(eXbeta * t(apply(X_sort, 1, tcrossprod)), ncol=p^2)
  num_beta1 <- t(colSums(num_beta1) - t(apply(num_beta1, 2, cumsum))) + num_beta1
  I22 <- matrix(colSums(sort_delta * (num_beta1 / eXbeta_sum - num_beta2_2 / (eXbeta_sum^2))),
                p, p )
  
  num_betarho01 <- eXbeta_X * zstar
  num_betarho01 <- t(colSums(num_betarho01) - t(apply(num_betarho01, 2, cumsum))) + num_betarho01
  num_betarho02 <- eXbeta*ztilde
  num_betarho02 <- t(colSums(num_betarho02) - t(apply(num_betarho02, 2, cumsum))) + num_betarho02
  num_betarho03 <- num_beta2 * zstar_eXbeta_sum
  
  I12 <- - drop(colSums(sort_delta * (ztilde - num_betarho01 / eXbeta_sum - num_betarho02 / eXbeta_sum + num_betarho03 / eXbeta_sum^2)))
  
  lam0_hat <- sort_delta / eXbeta_sum
  s11 <- sort_delta * (zstar - zstar_eXbeta_sum / eXbeta_sum)
  s12 <- zstar * eXbeta * cumsum(lam0_hat)
  s13 <- cumsum(zstar_eXbeta_sum / eXbeta_sum * lam0_hat) * eXbeta
  phistar_1 <- s11 - s12 + s13
  
  s21 <- sort_delta * (X_sort - num_beta2 / eXbeta_sum)
  s22 <- X_sort * eXbeta * cumsum(lam0_hat)
  s23 <- eXbeta * apply(num_beta2 / eXbeta_sum * lam0_hat, 2, cumsum)
  phistar_2 <- s21 - s22 + s23
  phistar <- drop(phistar_1 -  I12 %*% solve(I22,t(phistar_2)))
  Vs <- mean(phistar^2)
  stat <- score^2 / (n * Vs)
  return(list(score=score, Vs=Vs, stat=stat, phistar=phistar))
}

critical <- function(i){
  # Calculate one sample of the pertubed test statistic T_n^*
  # Input i: index of the number of samples.
  rn <- rnorm(n)
  sample <- max((colSums(rn * allphi))^2 / (allVs * n))
  return(sample)
}

score_test <- function(low, up, ngrid, nsample){ 
  # Calculate the proposed score test statistic 
  # Input low: lower bound of grid for gamma
  # Input up: upper bound of grid for gamma
  # Input ngrid: number of points to search for each gamma.
  # Input nsample: number of samples generated for T_n^*.
  # Output teststat: value of proposed statistic
  # Output critical: The 0.95 critical value
  # Output reject: The test result; whether the hypothesis is rejected.
  
  ## Grid for gamma. 
  ## This is an example of dim=2. More dimensions can be easily extended.
  grid_gamma1 <- seq(low,up,length=ngrid)
  grid_gamma2 <- seq(low,up,length=ngrid)
  grid_gamma3 <- seq(low,up,length=ngrid)
  n1 <- length(grid_gamma1)
  n2 <- length(grid_gamma2)
  n3 <- length(grid_gamma3)
  
  test <- foreach(a=rep(grid_gamma1, each=n2*n3), 
                b=rep(rep(grid_gamma2, each=n3), n1),
                d=rep(grid_gamma3, n1*n2), 
                .packages=c("survival"), 
                .export=c("score_test1", "critical", "X", "Y", "NET", "n")) %dopar% 
    score_test1(c(a, b, d))
  test_result <- max(unlist(lapply(test,`[[`,"stat")))
  allphi <- matrix(unlist(lapply(test,`[[`,"phistar")), nrow=n) 
  allVs <- unlist(lapply(test,`[[`,"Vs"))
  
  samples <- unlist(foreach(i=1:nsample, .export=c("critical", "n")) %dopar% critical(i))
  critical <- quantile(samples, prob=0.95)
  reject <- (test_result > critical)
  stopCluster(cl)
  return(list(teststat=test_result, critical=critical, reject=reject))
}
