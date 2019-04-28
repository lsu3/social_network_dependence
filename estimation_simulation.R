#==============================
#   Set up working directory
#==============================
code.dir <- './'
source(paste0(code.dir, "data_generator.R"))
source(paste0(code.dir, "fun_estimation.R"))
source(paste0(code.dir, "fun_variance.R"))

#===============================
#  Set up simulation settings
#===============================
# Generate network
n <- 2000
NET <- scb_network(p=c(0.01, 0.1, 0.05, 0.2, 0.1), p_prime=1e-4, n=c(500, 500, 400, 400, 200))

# set up parameters for responses
rho0 <- 0.1
gamma <- c(0, 1, -1) # the first is intercept
BETA <- c(1, -1)
lambda0 <- 0.5
q <- length(gamma)
c <- 2.5  # change this to control censor proportion

# simulation setup
S <- 1 # number of replicates
rho0_ini <- 0
gamma_ini <- c(0, 0, 0) # initial of p_i^t are all 0.5
k <- length(gamma_ini)
d <-5/n
beta_hat <- matrix(0, S, q-1)
rho0_hat <- rep(0, S)
gamma_hat <- matrix(0, S, q)
truep <- matrix(0, S, n)
truexi <- matrix(0, S, n)
phat <- matrix(0, S, n)
Athat <- matrix(0, S, n)
THETA <- c(rho0, BETA, gamma)
cover <- rep(0, 2*q)
sd_hat <- matrix(0, S, 2*q)
censor_rate <- rep(0, S)
iteration <- rep(0, S)

#===============================
#  Start simulation
#===============================
set.seed(0506)
start <- Sys.time()
for (s in 1:S){
  print(s)
  data <- generate_data(rho0, gamma, BETA, lambda0, c, NET)
  X <- data$X
  XX <- data$XX
  truep[s, ] <- data$p
  truexi[s, ] <- data$xi
  tildeT <- data$tildeT
  Y <- data$Y
  delta <- data$delta
  rk_tildeT <- rank(tildeT)
  order_tildeT <- order(tildeT)
  sort_delta <- delta
  sort_delta[rk_tildeT] <- delta
  censor_rate[s] <- data$censor_rate
  
  fit1 <- coxph(Y ~ X)
  beta_ini <- summary(fit1)$coef[, 1]
  p <- length(beta_ini)
  result <- EM(rho0_ini, beta_ini, gamma_ini, rep(0.5, n))
  beta_hat[s, ] <- result$beta
  rho0_hat[s] <- result$rho0
  gamma_hat[s, ] <- result$gamma
  iteration[s] <- result$iteration
  
  result_sd <- sdhat(result$rho0, result$beta, 
                     result$gamma, result$A_t, d)
  
  sd_hat[s, ] <- result_sd
  
  THETA_hat <- c(rho0_hat[s], beta_hat[s, ], gamma_hat[s, ])
  CI_upper <- THETA_hat + 1.96*(sd_hat[s, ])
  CI_lower <- THETA_hat - 1.96*(sd_hat[s, ])
  cover <- cover + as.numeric((THETA >= CI_lower)&(THETA <= CI_upper))
  print(cover)  
}
end <- Sys.time()
print(end-start)

#===============================
#  Results summary
#===============================
c(mean(rho0_hat), colMeans(beta_hat), colMeans(gamma_hat))
c(sd(rho0_hat), sqrt(apply(beta_hat, 2, "var")), sqrt(apply(gamma_hat, 2, "var")))
colMeans(sd_hat)
cover / S
mean(censor_rate)

# Average number of iterations
mean(iteration)
