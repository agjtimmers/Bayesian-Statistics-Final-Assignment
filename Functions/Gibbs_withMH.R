#################################################
## Gibbs sampler with Metropolis-Hastings step ##
#################################################

#make function for Gibbs sampler with Metropolis-Hastings step 
Gibbs_withMH <- function(Y, X1, X2, X3, n, n.iters = 10000, burnin = 1000, #iterations and burn-in in default
                         mu0 = 0.001, s20 = 1000, #set all hyperparameters already vague as default
                         mu1 = 0.001, s21 = 1000, df = 1,
                         mu2 = 0.001, s22 = 1000, 
                         mu3 = 0.001, s23 = 1000, 
                         a0 = 0.001, b0 = 0.001,
                         beta0, beta1, beta2, beta3, s2){
  #get everything ready
  parameters <- matrix(0, n.iters, 5) #storage
  colnames(parameters) <- c("b0", "b1", "b2", "b3", "variance") #names
  parameters[1,] <- c(beta0, beta1, beta2, beta3, s2) #first row is initial values
  
  if (beta3 == 0) {
    for(i in 2:n.iters){
      #conditional posterior beta0
      m0.nom   <- sum(Y - X1*beta1 - X2*beta2)/s2 + mu0/s20 #nominator 
      m0.denom <- n/s2 + 1/s20 #denominator
      m0.post  <- m0.nom/m0.denom #and divide to get the mean
      s0.post  <- 1/m0.denom #get the variance to put in (sqrt)
      beta0    <- rnorm(1, m0.post, sqrt(s0.post)) #and draw the next value
      
      #dependent metropolis-hastings beta1
      currentb1  <- parameters[i-1, 2] #current value
      proposedb1 <- rnorm(1, currentb1, 1) #proposed value normal draw
      dcur       <- exp((-(currentb1^2)  * (sum(X1^2) / (2*s2))) + (currentb1  * (sum(X1 * (Y - beta0 - X2*beta2)) / s2))) * (1 + ((currentb1  - mu1)^2 / (df * (s21^2))))^(-((df + 1)/2))
      dpro       <- exp((-(proposedb1^2) * (sum(X1^2) / (2*s2))) + (proposedb1 * (sum(X1 * (Y - beta0 - X2*beta2)) / s2))) * (1 + ((proposedb1 - mu1)^2 / (df * (s21^2))))^(-((df + 1)/2))
      alpha      <- (dpro / dcur) * (currentb1 / proposedb1) #acceptance probability
      ref        <- runif(1, 0, 1) #reference to which I compare
      beta1      <- ifelse(ref <= alpha, proposedb1, currentb1) #update or retain current value
      
      #conditional posterior beta2
      m2.nom   <- sum(X2*(Y - beta0 - X1*beta1))/s2 + mu2/s22
      m2.denom <- sum(X2^2)/s2 + 1/s22
      b2.post  <- m2.nom/m2.denom
      s2.post  <- 1/m2.denom
      beta2    <- rnorm(1, b2.post, sqrt(s2.post))
      
      #conditional posterior variance
      alpha <- n/2 + a0 #formula alpha slide
      beta  <- sum((Y - beta0 - X1*beta1 - X2*beta2)^2)/2 + b0 #formula beta slide
      s2    <- 1/rgamma(1, shape = alpha, rate = beta) #and fill in to draw
      
      #and keep it
      parameters[i,] <- c(beta0, beta1, beta2, beta3, s2)
    }
  }
  else {
    for(i in 2:n.iters){
      m0.nom   <- sum(Y - X1*beta1 - X2*beta2 - X3*beta3)/s2 + mu0/s20 #nominator 
      m0.denom <- n/s2 + 1/s20 #denominator
      m0.post  <- m0.nom/m0.denom #and divide to get the mean
      s0.post  <- 1/m0.denom #get the variance to put in (sqrt)
      beta0    <- rnorm(1, m0.post, sqrt(s0.post)) #and draw the next value
      
      #dependent metropolis-hastings for beta1
      currentb1  <- parameters[i-1, 2] #current value
      proposedb1 <- rnorm(1, currentb1, 1) #proposed value normal draw
      dcur       <- exp((-(currentb1^2)  * (sum(X1^2) / (2*s2))) + (currentb1  * (sum(X1 * (Y - beta0 - X2*beta2 - X3*beta3)) / s2))) * (1 + ((currentb1  - mu1)^2 / (df * (s21^2))))^(-((df + 1)/2))
      dpro       <- exp((-(proposedb1^2) * (sum(X1^2) / (2*s2))) + (proposedb1 * (sum(X1 * (Y - beta0 - X2*beta2 - X3*beta3)) / s2))) * (1 + ((proposedb1 - mu1)^2 / (df * (s21^2))))^(-((df + 1)/2))
      alpha      <- (dpro / dcur) * (currentb1 / proposedb1) #acceptance probability
      ref        <- runif(1, 0, 1) #reference to which we compare
      beta1      <- ifelse(ref <= alpha, proposedb1, currentb1) #update or retain current value
      
      #conditional posterior beta2
      m2.nom   <- sum(X2*(Y - beta0 - X1*beta1 - X3*beta3))/s2 + mu2/s22
      m2.denom <- sum(X2^2)/s2 + 1/s22
      b2.post  <- m2.nom/m2.denom
      s2.post  <- 1/m2.denom
      beta2    <- rnorm(1, b2.post, sqrt(s2.post))
      
      #conditional posterior beta2
      m3.nom   <- sum(X3*(Y - beta0 - X1*beta1 - X2*beta2))/s2 + mu3/s23
      m3.denom <- sum(X3^2)/s2 + 1/s23
      b3.post  <- m3.nom/m3.denom
      s3.post  <- 1/m3.denom
      beta3    <- rnorm(1, b3.post, sqrt(s3.post))
      
      #conditional posterior variance
      alpha  <- n/2 + a0 #formula alpha slide
      beta   <- sum((Y - beta0 - X1*beta1 - X2*beta2 - X3*beta3)^2)/2 + b0 #formula beta slide
      s2     <- 1/rgamma(1, shape = alpha, rate = beta) #and fill in to draw
      
      #and keep it
      parameters[i,] <- c(beta0, beta1, beta2, beta3, s2)
    }
  } 
  #remove burn-in period
  parameters_burnin <<- parameters[-c(1:burnin),]
  return(parameters_burnin)
}
