### POSTERIOR PREDICTIVE CHECK ###
#function to compute mode, used in the function below
mode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

#normality of the residuals -> mean = median = mode, discrepancy measure is difference mean and mode
#note this function can only be used on a model with two predictors, because it was not written to be generally applicable!
ppp <- function(rchains, Y, X1, X2, n = length(X1)) { #input is dataframe with the predicted parameters of the number of chains bound by the rows, 
  #a vector containing the dependent variable, two vectors, one for each predictor, and the number of cases (set to default on the number of entries in X1)
  obs_mode.mean.dif <- rep(0, nrow(rchains)) #storage 
  sim_mode.mean.dif <- rep(0, nrow(rchains)) #storage 
  #sample replicated datasets and compute discrepancy measure 
  for (i in 1:nrow(rchains)){
    y <- rchains[i,1] + X1 * rchains[i,2] + X2 * rchains[i,3] #predicted y values
    residuals <- Y - y #observed residuals
    obs_mode.mean.dif[i] <- mean(residuals) - mode(residuals) #difference mean and mode observed
    ppp.y <- rnorm(n, rchains[i,1] + X1 * rchains[i,2] + X2 * rchains[i,3], sqrt(rchains[i,5])) #generate datasets
    ppp.residuals <- ppp.y - (rchains[i,1] + X1 * rchains[i,2] + X2 * rchains[i,3]) #simulated residuals 
    sim_mode.mean.dif[i] <- mean(ppp.residuals) - mode(ppp.residuals) #difference mean and mode simulated
  }
  ppp.value <- mean(ifelse(abs(sim_mode.mean.dif) > abs(obs_mode.mean.dif), 1, 0)) #compare absolute values!
  return(ppp.value)
}


