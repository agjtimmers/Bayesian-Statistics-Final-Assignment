##############################
## Convergence Diagnostics ##
#############################

## Autocorrelation plots ##
#make function that calculates the autocorrelations for me
autocorrelation <- function(x,lag=50){
  x  <- as.matrix(as.ts(x))
  meanx  <- mean(x)
  n  <- as.integer(nrow(x))
  ac <- as.numeric(sapply(0:lag,function(k)
    sum((x[1:(n-k)] - meanx) * (x[(k+1):n] - meanx)) / sum((x - meanx)^2)))
  return(ac)
}

## Running mean plots ##
#make function that calculates the running means for me
runmean <- function(x, n = (n.iters - burnin)) {
  meanx <- numeric(n) #storage
  for (k in 1:n) {
    meanx[k] <- mean(x[1:k]) #mean, add one iteration at a time
  }
  return(meanx)
}

## MC errors ##
#make function that calculates the MC errors for me
mc_error <- function(sd, n.iters) {
  mc_error <- sd / sqrt(n.iters)
  return(mc_error)
}