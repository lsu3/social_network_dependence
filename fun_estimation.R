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


Lambda0t <- function(A_t1, rho0_t, beta_t){
  # Calculate the baseline cumulative harzard function given A, rho0 and beta
  # Input A_t1, rho0_t, beta_t: values of A, rho0, beta at the (t)th iteration
  rk_tildeT <- rank(Y[, 1])
  order_tildeT <- order(Y[, 1])
  sort_delta <- Y[, 2]
  sort_delta[rk_tildeT] <- Y[, 2]
  
  Xbeta_t <- X %*% beta_t
  den <- exp(Xbeta_t) * ((1 - A_t1) + A_t1 * exp(rho0_t * (NET %*% Xbeta_t)))
  den[rk_tildeT] <- den
  den <- sum(den) - cumsum(den) + den
  lam0 <- sort_delta / den
  Lambda0 <- cumsum(lam0)
  lam0[order_tildeT] <- lam0
  Lambda0[order_tildeT] <- Lambda0
  return(list(Lambda0=Lambda0,lam0=lam0))
}

At <- function(beta_t, rho0_t, gamma_t, Lambda0_t){
  # Calculate the posterior probability of xi, which is denoted by A
  # Input beta_t, rho0_t, gamma_t, Lambda0_t: values of other parameters at the 
  #                                           (t)th interation.
  delta <- Y[, 2]
  p_t <- pt(gamma_t, X)
  dependence <- rho0_t * (NET %*% X %*% beta_t)
  num <- exp(delta * dependence-exp(X %*% beta_t + dependence) * Lambda0_t) * p_t
  den <- num + exp(-exp(X %*% beta_t) * Lambda0_t) * (1 - p_t)
  A_t <- num / den
  return(A_t)  
}

pl1 <- function(rho0, beta, A_t){
  # Calculate the log profile likelihood function (defined on Page 10 in the paper)
  rk_tildeT <- rank(Y[, 1])
  order_tildeT <- order(Y[, 1])
  delta <- Y[, 2]
  sort_delta <- Y[, 2]
  sort_delta[rk_tildeT] <- Y[, 2]
  
  Xbeta <- X %*% beta
  eXbeta <- exp(Xbeta)
  dependence <- rho0 * (NET %*% X %*% beta)
  logterm <- eXbeta *(1 - A_t + A_t * exp(dependence))
  logterm[rk_tildeT] <- logterm
  logterm <- sum(logterm)-cumsum(logterm) + logterm
  logterm[order_tildeT] <- logterm
  pl1 <- sum(delta * (Xbeta + A_t * dependence - log(logterm)))
  return(pl1)  
}


Q <- function(beta, rho0, gamma, lam0, Lambda0, A_t){
  # Calculate the Q function, which the posterior expectation function in E-step.
  # Definition can be found on Page 9 in the paper.
  delta <- Y[, 2]
  XX <- cbind(rep(1,n), X)
  l1 <- ifelse(delta==0, 0, delta * log(lam0))
  l1 <- l1 + delta * (X %*% beta + rho0 * A_t * (NET %*% X %*% beta)) -
    exp(X %*% beta) * Lambda0 * (1 - A_t + A_t * exp(rho0 * (NET %*% X %*% beta)))
  if(length(gamma)==1){
    l2 <- A_t * gamma - log(1 + exp(gamma))
  }else{
    l2 <- A_t * (XX %*% gamma ) - log(1 + exp(XX %*% gamma))
  }
  l <- l1 + l2
  return(sum(l))
}

logL <- function(beta, rho0, gamma, lam0, Lambda0){
  # Calculate the observed log likelihood function
  delta <- Y[, 2]
  p_t <- pt(gamma, X)
  l1 <- ifelse(delta == 0, 0, delta * log(lam0))
  l1 <- l1 + delta * (X %*% beta + rho0 * (NET %*% X %*% beta)) - 
    exp(X %*% beta) * Lambda0 * (exp(rho0 * (NET %*% X %*% beta)))
  l2 <- ifelse(delta == 0, 0, delta * log(lam0))
  l2 <- l2 + delta * (X %*% beta) - exp(X %*% beta) * Lambda0
  return(sum(log((exp(l1) * p_t + exp(l2) * (1 - p_t)))))
}

EM <- function(rho0, beta, gamma, A_t, maxiter=1000, tol1=1e-6, tol2=1e-6){
  # The proposed EM algorithm to get the parameter estimates.
  # Input rho0, beta, gamma, A_t: Initial values for the prameters.
  k <-length(gamma)
  p <- length(beta)
  
  rk_tildeT <- rank(Y[, 1])
  
  for(i in 1:maxiter){
    Lambda0_t <- Lambda0t(A_t, rho0, beta)
    A_t <- At(beta, rho0, gamma, Lambda0_t$Lambda0)
    
    # update gamma
    if(length(gamma) == 1){
      logfit <- glm(A_t~1, family="quasibinomial")
    }else{
      logfit <- glm(A_t~X, family="quasibinomial")
    }
    
    gamma_t <- summary(logfit)$coeff[ ,1]
    
    
    #update rho0
    rho0_t <- optimize(pl1, interval=c(-1, 1), beta=beta, A_t=A_t, maximum=T)$maximum
    
    # update beta
    C_t <- A_t * exp(rho0_t * (NET %*% X %*% beta)) + 1-A_t
    cox <- coxph(Y ~ X + offset(log(C_t)))
    beta_t <- summary(cox)$coef[, 1]
    
    if(max(abs(gamma_t-gamma))<tol1 &
         abs(rho0_t-rho0)<tol2 & 
         max(abs(beta_t-beta))<tol2){
      log_likelihood <- logL(beta_t, rho0_t, gamma_t,  Lambda0_t$lam0,  Lambda0_t$Lambda0)
      break
    }
    else{
      beta <- beta_t
      rho0 <- rho0_t
      gamma <- gamma_t
      log_likelihood <- logL(beta_t, rho0_t, gamma_t,  Lambda0_t$lam0,  Lambda0_t$Lambda0)
    }
  }
  return(list(beta=beta_t, rho0=rho0_t, gamma=gamma_t,  A_t=A_t, iteration=i,
              log_likelihood=log_likelihood))  
}

