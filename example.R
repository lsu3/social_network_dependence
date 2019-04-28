#==============================
#     Load in packages
#==============================
library(parallel)
library(foreach)
library(doParallel)
library(survival)

#==============================
#   Set up working directory
#==============================
code.dir = './'
source(paste0(code.dir, "fun_test.R"))
source(paste0(code.dir, "fun_estimation.R"))
source(paste0(code.dir, "fun_variance.R"))
source(paste0(code.dir, "data_prep.R"))

#===============================
#   Read in Data
#===============================
network <- read.table("network.csv", sep=",", header=T)  # read in the network data
network <- network[, -1] # delete the row number

attr <- read.table("attributes.csv", sep=",", header=T) # read in the covariates data
attr <- attr[, -1] # delete the row number
n <- dim(attr)[1]  

NET <- adj(n, network)
data <- getXY(attr)
X <- data$X
Y <- data$Y

#===============================
#   Test
#===============================
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)
# grid for gamma
grid_gamma1 <- seq(-2,2,length=9)
grid_gamma2 <- seq(-2,2,length=9)
grid_gamma3 <- seq(-2,2,length=9)
n1 <- length(grid_gamma1)
n2 <- length(grid_gamma2)
n3 <- length(grid_gamma3)
# get test statistic
test <- foreach(a=rep(grid_gamma1, each=n2*n3), 
              b=rep(rep(grid_gamma2, each=n3), n1),
              d=rep(grid_gamma3, n1*n2), 
              .packages=c("survival")) %dopar% 
  score_test1(c(a, b, d))
test_result <- max(unlist(lapply(test,`[[`,"stat")))
# get resampling critical value
allphi <- matrix(unlist(lapply(test,`[[`,"phistar")), nrow=n) 
allVs <- unlist(lapply(test,`[[`,"Vs"))
samples <- unlist(foreach(i=1:1000) %dopar% critical(i))
critical <- quantile(samples, prob=0.95)
reject <- (test_result > critical)
stopCluster(cl)

print(paste("test stat", test_result))
print(paste("critical value", critical))
print(paste("reject?", reject))
# > print(paste("test stat", test_result))
# [1] "test stat 51.5756628715865"
# > print(paste("critical value", critical))  # random result
# [1] "critical value 7.47247864792162"
# > print(paste("reject?", reject))
# [1] "reject? TRUE"

#===============================
#   Parameter Estimation
#===============================
## Set up for parameter estimation
rho0_ini <- 0
gamma_ini <- c(0, 0, 0) 
fit1 <- coxph(Y ~ X)
beta_ini <- summary(fit1)$coef[, 1]

## Run parameter estimation
result_parameter <- EM(rho0_ini, beta_ini, gamma_ini, rep(0.5, n))
print(result_parameter$rho0)
print(result_parameter$beta)
print(result_parameter$gamma)
# > print(result_parameter$rho0)
# [1] 0.04769068
# > print(result_parameter$beta)
# XX1        XX2 
# 0.9896974 -0.9912256 
# > print(result_parameter$gamma)
# (Intercept)         XX1         XX2 
# -0.2698215   0.8200203  -1.1638705 

#=====================================
#   Variance Estimation
#=====================================
## Set up for variance estimation
d <- 5/n
result_sd <- sdhat(result_parameter$rho0, result_parameter$beta, 
                   result_parameter$gamma, result_parameter$A_t, d)
print(result_sd)
# > print(result_sd)
# [1] 0.007665225 0.098483325 0.085008835 0.284343432 0.396145569 0.344841522
