EM1 <- function(rho0, beta, gamma, A_t){
  # One-step estimate of the proposed EM algorithm, used in function gradM.
  # update gamma
  if(length(gamma)==1){
    logfit <- glm(A_t~1, family="quasibinomial")
  }else{
    logfit <- glm(A_t~X, family="quasibinomial")
  }
  
  gamma_t <- summary(logfit)$coeff[ ,1]
  
  #update rho0
  rho0_t <- optimize(pl1, interval=c(-1, 1), beta=beta, A_t=A_t, maximum=T)$maximum
  # update beta
  C_t <- A_t * exp(rho0_t * (NET %*% X %*% beta)) + 1 - A_t
  cox <- coxph(Y ~ X + offset(log(C_t)))
  beta_t <- summary(cox)$coef[, 1]
  
  return(list(beta=beta_t, rho0=rho0_t, gamma=gamma_t))  
}

EM_var <- function(rho0_hat, beta_hat, gamma_hat, A_t1, maxiter=300, tol=1e-6){
  # Calculate the variable value of A and Lambda0 used in the one-step EM.
  # This function iterates between A and Lambda0 till convergence, given other parameters.
  Lambda0_t1 <- Lambda0t(A_t1, rho0_hat, beta_hat)$Lambda0
  for(i in 1:maxiter){
    A_t <- At(beta_hat, rho0_hat, gamma_hat, Lambda0_t1)
    Lambda0_t <- Lambda0t(A_t, rho0_hat, beta_hat)$Lambda0
    if(max(abs(Lambda0_t-Lambda0_t1)) < tol){
      break
    }
    else{
      Lambda0_t1 <- Lambda0_t
    }
  }
  
  lam_hat <- Lambda0t(A_t, rho0_hat, beta_hat)
  Lambda0_hat <- lam_hat$Lambda0
  lam0_hat <- lam_hat$lam0
  l <- Q(beta_hat, rho0_hat, gamma_hat, lam0_hat, Lambda0_hat, A_t)
  return(list(l=l, A_t=A_t))
}

LapQ <- function(rho0_hat, beta_hat, gamma_hat, A_t1, k, p){
  # Calculate the negative of the second derivative of function g, defined on Page 12
  rk_tildeT <- rank(Y[, 1])
  order_tildeT <- order(Y[, 1])
  delta <- Y[, 2]
  sort_delta <- Y[, 2]
  sort_delta[rk_tildeT] <- Y[, 2]
  # the laplacian matrix is in the order of rho0, beta, gamma
  Lap <- matrix(0, p+k+1, p+k+1)
  Xbeta <- X %*% beta_hat
  Hj <- drop(NET %*% Xbeta)
  Gj <- NET %*% X
  Dj <- drop(exp(Xbeta) * (1 - A_t1))
  Ej <- drop(exp(Xbeta) * A_t1 * exp(rho0_hat * Hj))
  Mj <- X + rho0_hat * Gj
  
  Hj[rk_tildeT] <- Hj
  Gj[rk_tildeT, ] <- Gj
  Dj[rk_tildeT] <- Dj
  Ej[rk_tildeT] <- Ej
  Mj[rk_tildeT, ] <- Mj
  X_sort <- X
  X_sort[rk_tildeT, ] <- X
  A_sort <- drop(A_t1)
  A_sort[rk_tildeT] <- A_sort
  
  C <- sum(Dj + Ej) - cumsum(Dj + Ej) +(Dj + Ej)
  # second derivative wrt rho0
  num_rho0 <- (sum(Ej * Hj^2)-cumsum(Ej * Hj^2)+(Ej * Hj^2)) * C-(sum(Ej * Hj) - cumsum(Ej * Hj)+(Ej * Hj))^2
  Lap[1,1] <- -sum(sort_delta * num_rho0 / C^2)
  
  # second derivative wrt beta
  DX_EM <- Dj * X_sort + Ej * Mj
   
  num_beta2 <- t(colSums(DX_EM) - t(apply(DX_EM, 2, cumsum))) + DX_EM
  num_beta2 <- matrix(t(apply(num_beta2, 1, tcrossprod)), ncol=p^2)
  
  num_beta1 <- matrix(Dj * t(apply(X_sort, 1, tcrossprod)) + Ej * t(apply(Mj, 1, tcrossprod)), ncol=p^2)
  num_beta1 <- t(colSums(num_beta1) - t(apply(num_beta1, 2, cumsum))) + num_beta1
  Lap[2:(p+1), 2:(p+1)] <- matrix(-colSums(sort_delta * (num_beta1 * C-num_beta2) / C^2), p, p)
  
  
  # derivative wrt beta and rho0
  num_betarho01 <- Ej * (Hj * Mj + Gj)
  num_betarho01 <- t(colSums(num_betarho01) - t(apply(num_betarho01, 2, cumsum))) + num_betarho01
  num_betarho02 <- Dj * X_sort + Ej * Mj
  num_betarho02 <- t(colSums(num_betarho02) - t(apply(num_betarho02, 2, cumsum))) + num_betarho02
  num_betarho03 <- sum(Ej * Hj) - cumsum(Ej * Hj) + Ej * Hj
  
  Lap[2:(p+1), 1] <- Lap[1, 2:(p+1)] <- 
    drop(colSums(sort_delta * A_sort * Gj - sort_delta * (num_betarho01 * C - num_betarho02 * num_betarho03) / C^2))
  XX <- cbind(rep(1,n), X)
  # second derivative wrt gamma
  if(length(gamma_hat) == 1){
    e_power <- rep(exp(gamma_hat), n)
  }else{
    e_power <- drop(exp(XX %*% gamma_hat))
  }
  Lap[(p+2):(p+1+k), (p+2):(p+1+k)] <- matrix(-colSums((e_power / (1 + e_power)^2) * t(apply(XX, 1, tcrossprod))), k, k)
  
  return(Lap)
}

gradM <- function(THETA_hat, A_t1, d, k, p){
  # Calculate the derivative of function M defined on Page 12
  #THETA_hat is in the order rho0, beta, gamma
  gradM <- matrix(0, p+k+1, p+k+1)
  for(j in 1:(p+k+1)){
    THETA_plus <- THETA_hat
    THETA_plus[j] <- THETA_hat[j] + d
    A_t <- EM_var(THETA_plus[1], THETA_plus[2:(p+1)], THETA_plus[(p+2):(p+k+1)], A_t1)$A_t
    
    M_plus <- EM1(THETA_plus[1], THETA_plus[2:(p+1)], THETA_plus[(p+2):(p+k+1)], A_t)
    M_plus <- c(M_plus$rho0, M_plus$beta, M_plus$gamma)
    gradM[, j] <- (M_plus-THETA_hat)/d
  }
  return(gradM)  
}

sdhat <- function(rho0_hat, beta_hat, gamma_hat, A_t, d){
  # Calculate the standard errors of parameters.
  k <- length(gamma_hat)
  p <- length(beta_hat)
  lap <- LapQ(rho0_hat, beta_hat,gamma_hat, A_t, k, p)
  lapl <- lap%*%(diag(p+k+1)-gradM(c(rho0_hat, beta_hat, gamma_hat), A_t, d, k, p))
  sd_hat <- sqrt(diag(solve(-lapl)))
  return(sd_hat)
}