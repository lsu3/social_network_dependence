#==============================
#     Load in packages
#==============================
library(parallel)
library(foreach)
library(doParallel)
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)

#==============================
#   Set up working directory
#==============================
code.dir <- './'
source(paste0(code.dir, "fun_test.R"))
source(paste0(code.dir, "data_generator.R"))

#===============================
#  Set up simulation settings
#===============================
# Generate network
n <- 2000
NET <- scb_network(p=c(0.01, 0.1, 0.05, 0.2, 0.1), p_prime=1e-4, n=c(500, 500, 400, 400, 200))

# set up parameters for responses
rho0 <-  0.02
gamma <- c(0, 1, -1) # the first is intercept
BETA <- c(1, -1)
lambda0 <- 0.5
q <- length(gamma)
c <- 9 # change this to control censor proportion

# simulation setup
S <- 100 # number of replicates
censor_rate <- rep(0, S)
total_reject <- 0
all_test <- rep(0, S)

# Grid of gramma
grid_gamma1 <- seq(-2,2,length=9)
grid_gamma2 <- seq(-2,2,length=9)
grid_gamma3 <- seq(-2,2,length=9)
n1 <- length(grid_gamma1)
n2 <- length(grid_gamma2)
n3 <- length(grid_gamma3)

#===============================
#  Start simulation
#===============================
set.seed(1222)
start_time <- Sys.time()
for (s in 1:S){
  print(s)
  data <- generate_data(rho0, gamma, BETA, lambda0, c, NET)
  X <- data$X
  XX <- data$XX
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
  
  test <- foreach(a=rep(grid_gamma1, each=n2*n3), 
                b=rep(rep(grid_gamma2, each=n3), n1),
                d=rep(grid_gamma3, n1*n2),
                .packages=c("survival")) %dopar% 
    score_test1(c(a, b, d))
  test_result <- max(unlist(lapply(test,`[[`,"stat")))
  all_test[s] <- test_result
  allphi <- matrix(unlist(lapply(test,`[[`,"phistar")), nrow=n) # each column is for each theta
  allVs <- unlist(lapply(test,`[[`,"Vs"))
  
  samples<- unlist(foreach(i=1:1000) %dopar% critical(i))
  reject <- as.numeric(test_result > quantile(samples, prob=0.95))
  total_reject <- total_reject + reject
  print(paste("test_result", test_result))
  print(paste("quantile", quantile(samples, prob=0.95)))
  print(paste("total_reject",total_reject))
}
end_time <- Sys.time()
stopCluster(cl)

print(end_time - start_time)


