library(survival)
adj <- function(n, df){
  NET <- matrix(0, n, n)
  for(i in 1:dim(df)[1]){
    NET[df[i, 1], df[i, 2]] <- 1
  }
  return(NET)
}

getXY <- function(df){
  X <- as.matrix(df[, c("X1", "X2")])
  Y <- Surv(attr$tildeT, attr$delta)
  return(list(X=X, Y=Y))
}