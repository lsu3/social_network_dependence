library(ergm)
library(network)
library(Matrix)
library(igraph)
library(survival)

scb_network <- function(p, p_prime, n){
  # Generate network from stochastic block model
  # p: list of within group probabilities for each block
  # p_prime: between group probability
  # n: list of # of nodes for each block
  if(length(p) != length(n)){
    stop("Length of p and n do not match!")
  }
  n_block <- length(p)
  n_total <- sum(n)
  
  blocks <- list()
  for(i in 1:n_block){
    blocks[[i]] <- as.matrix(network(n[i], directed=FALSE, density=p[i]))  
  }
  
  for (i in 1:n_block){
    blocks[[i]][blocks[[i]]==0] <- 1000 # 1000 is an indicator
  }
  NET <- bdiag(blocks)
  NET <- as.matrix(NET)
  NET[NET==0] <- rbinom(sum(NET==0,na.rm=TRUE), 1, p_prime)
  NET[NET==1000] <- 0
  
  # Make NET symmetric
  NET_diag <- diag(NET)
  NET[lower.tri(NET, diag=T)] <- 0
  NET <- NET + t(NET) + diag(NET_diag)
  
  return(NET)
}

# Example:
# W <- scb_network(p=c(0.01, 0.1, 0.05, 0.2, 0.1), p_prime=1e-4, n=c(500, 500, 400, 400, 200))


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

generate_data <- function(rho0, gamma, BETA, lambda0, c, NET){
  # Generate data used in simulation
  X1 <- rbinom(n, 1, 0.5)
  X2 <- runif(n, -1, 1)
  X <- cbind(X1, X2)
  XX <- cbind(rep(1,n), X)
  p <- pt(gamma, X)
  xi <- rbinom(n, size = 1, p)
  rho <- rho0 * xi 
  lambda1 <- lambda0 * exp(X %*% BETA + rho * (NET %*% X %*% BETA))
  t <- rexp(n, rate = lambda1)
  C <- runif(n, 0, c)
  tildeT <- apply(cbind(C, t), 1, 'min')
  delta <- (t < C)
  N <- sum(delta==1) # Number of observed events
  Y <- Surv(tildeT, delta)
  return(list(X=X, XX=XX, p=p, xi=xi, tildeT=tildeT, Y=Y, delta=delta, censor_rate=(n-N)/n))
}
