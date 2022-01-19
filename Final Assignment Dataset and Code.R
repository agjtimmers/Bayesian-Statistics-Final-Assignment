##########################################
## Final Assignment Bayesian Statistics ##
## Annemarie Timmers (6238106) ###########
##########################################

### PREPARATION ###
#load the libraries
lapply(c("mlbench", "psych", "ggplot2", "cowplot", "reshape2", "bain", "BayesFactor"), require, character.only = TRUE)
#and the dataset
data("BostonHousing")

### INTRODUCTION DATASET ###
## Table 1 - Descriptive statistics ##
describe(BostonHousing[, c(14, 6, 13)])

## Table 2 - Bivariate associations ##
cor(BostonHousing[, c(14, 6, 13)]) 

### ESTIMATION ###
## Function ##
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

## Prepare input ##
#center variables
BostonHousing$rm.cent <- scale(BostonHousing$rm, scale = F)
BostonHousing$lstat.cent <- scale(BostonHousing$lstat, scale = F)
#assign input for the function
Y  <- as.matrix(BostonHousing$medv) #dependent variable
X1 <- as.matrix(BostonHousing$rm.cent) #predictor 1
X2 <- as.matrix(BostonHousing$lstat.cent) #predictor 2
n  <- length(BostonHousing$medv) #number of observations

## Run two chains ##
set.seed(10) #set a seed so I get the some results when I run it another time
chain1 <- Gibbs_withMH(Y, X1, X2, X3 = 0, n, #no need for the default parameters
                       beta0 = 5, beta1 = 0.1, beta2 = 0.1, beta3 = 0, s2 = 100)                       
chain2 <- Gibbs_withMH(Y, X1, X2, X3 = 0, n,
                       beta0 = 50, beta1 = 2, beta2 = 2, beta3 = 0, s2 = 50) #different initial values
#look at the results quickly, this is also used to calculate the Gelman-Rubin statistics
est_ch1 <- data.frame(est = apply(chain1[, c(1:3, 5)], 2, mean), sd = apply(chain1[, c(1:3, 5)], 2, sd)) #seems good
est_ch2 <- data.frame(est = apply(chain2[, c(1:3, 5)], 2, mean), sd = apply(chain2[, c(1:3, 5)], 2, sd)) #seems good

## Calculate acceptance ratios ##
n.iters <- 10000 #10,000
burnin  <- 1000
acceptance_ratios <- c(length(unique(chain1[,2])) / (n.iters - burnin),  length(unique(chain2[,2])) / (n.iters - burnin)) #good, like a wanted between .2 and .5

## Preparation convergence plots and results ## 
#get the chains in a dataframe together and order (names and stuff) (used for convergence plots)
chains <- as.data.frame(cbind(chain1, chain2)) 
colnames(chains) <- c("b0_ch1", "b1_ch1", "b2_ch1", "b3_ch1", "s2_ch1", "b0_ch2", "b1_ch2", "b2_ch2", "b3_ch2", "s2_ch2")
chains$iteration <- 1:nrow(chains) #add column indicating the iteration for the plots

#when there is convergence, it can be assumed they sample from the same distribution and can be combined (also used for covergence plots)
rchains <- data.frame(rbind(chain1, chain2), 
                      chain = c(rep(1, (n.iters - burnin)), rep(2, (n.iters - burnin))),
                      iteration = rep(1:(n.iters - burnin), 2))
est <- data.frame(est = apply(rchains[, c(1:3, 5)], 2, mean)) #posterior means of the combined chains, use later

### INTERPRETATION OF ESTIMATES AND INTERVALS ###
## Table 3 - Results Model 1 ## 
b  <- c(est[2,], est[3,]) #prep for standardized coefficients in output
sy <- sd(BostonHousing$medv) 
sx <- c(sd(BostonHousing$rm.cent), sd(BostonHousing$lstat.cent)) 
output <- data.frame(mean = round(apply(rchains[, c(1:3, 5)], 2, mean), 3), 
                     beta = c("", round(b * (sx/sy), 3), ""),
                     sd = round(c(apply(rchains[, c(1:3, 5)], 2, sd)), 3),
                     naive.se = round(apply(rchains[, c(1:3, 5)], 2, sd)/sqrt(n.iters), 3),
                     cred.low = round(apply(rchains[, c(1:3, 5)], 2, quantile, probs = 0.025), 3),
                     cred.high = round(apply(rchains[, c(1:3, 5)], 2, quantile, probs = 0.975), 3))

### CONVERGENCE ###
## Trace plots ##
#trace plot beta0 for both chains
traceb0 <- ggplot(chains, aes(x = iteration)) +
  geom_line(aes(y = b0_ch1), color = "red") + 
  geom_line(aes(y = b0_ch2), color = "blue") +
  ggtitle(bquote(bold("Trace plot of ")~beta[0])) + ylab(bquote(beta[0])) +
  theme(plot.title = element_text(hjust = 0.5)) #looks good

#trace plot beta1 for both chains
traceb1 <- ggplot(chains, aes(x = iteration)) +
  geom_line(aes(y = b1_ch1), color = "red") + 
  geom_line(aes(y = b1_ch2), color = "blue") +
  ggtitle(bquote(bold("Trace plot of ")~beta[1])) + ylab(bquote(beta[1])) +
  theme(plot.title = element_text(hjust = 0.5)) #looks good

#trace plot beta2 for both chains
traceb2<- ggplot(chains, aes(x = iteration)) +
  geom_line(aes(y = b2_ch1), color = "red") + 
  geom_line(aes(y = b2_ch2), color = "blue") +
  ggtitle(bquote(bold("Trace plot of ")~beta[2])) + ylab(bquote(beta[2])) +
  theme(plot.title = element_text(hjust = 0.5)) #looks good

#trace plot s2 for both chains
traces2 <- ggplot(chains, aes(x = iteration)) +
  geom_line(aes(y = s2_ch1), color = "red") + 
  geom_line(aes(y = s2_ch2), color = "blue") +
  ggtitle(bquote(bold("Trace plot of ")~sigma^2)) + ylab(bquote(sigma^2)) +
  theme(plot.title = element_text(hjust = 0.5)) #looks good

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

#calculate the autocorrelations and put them in a dataframe
autocor <- apply(chains[1:10], 2, autocorrelation) 
autocor <- data.frame(lag = 1:51, autocor)

#plot autocorrelation beta0
autb0 <- ggplot(data = autocor, mapping = aes(x = lag)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(y = b0_ch1, xend = lag, yend = 0), color = "red") + 
  geom_segment(mapping = aes(y = b0_ch2, xend = lag, yend = 0), color = "blue") +
  ggtitle(bquote(bold("Autocorrelation  plot of ")~beta[0])) + ylab(bquote(beta[0])) +
  theme(plot.title = element_text(hjust = 0.5)) #looks good

#plot autocorrelation beta1
autb1 <- ggplot(data = autocor, mapping = aes(x = lag)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(y = b1_ch1, xend = lag, yend = 0), color = "red") + 
  geom_segment(mapping = aes(y = b1_ch2, xend = lag, yend = 0), color = "blue") +
  ggtitle(bquote(bold("Autocorrelation plot of ")~beta[1])) + ylab(bquote(beta[1])) +
  theme(plot.title = element_text(hjust = 0.5)) #looks good

#plot autocorrelation beta2
autb2 <- ggplot(data = autocor, mapping = aes(x = lag)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(y = b2_ch1, xend = lag, yend = 0), color = "red") + 
  geom_segment(mapping = aes(y = b2_ch2, xend = lag, yend = 0), color = "blue") +
  ggtitle(bquote(bold("Autocorrelation  plot of ")~beta[2])) + ylab(bquote(beta[2])) +
  theme(plot.title = element_text(hjust = 0.5)) #looks good

#plot autocorrelation s2
auts2 <- ggplot(data = autocor, mapping = aes(x = lag)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(y = s2_ch1, xend = lag, yend = 0), color = "red") + 
  geom_segment(mapping = aes(y = s2_ch2, xend = lag, yend = 0), color = "blue") +
  ggtitle(bquote(bold("Autocorrelation  plot of ")~sigma^2)) + ylab(bquote(sigma^2)) +
  theme(plot.title = element_text(hjust = 0.5)) #looks good

## Trace and autocorrelation plots in one figure ##
plot_grid(traceb0, autb0, traceb1, autb1, traceb2, autb2, traces2, auts2, 
          align = "hv", nrow = 4, rel_widths = c(1/2, 1/2, 1/2, 1/2), rel_heights = c(1/4, 1/4, 1/4, 1/4)) #tip: click zoom

## Density plots ##
#density plot beta0 for both chains
chains_b0 <- as.data.frame(chains[, c(1, 6, 11)])
densb0 <- ggplot(melt(chains_b0, id.var='iteration', variable.name='chain'), aes(x=value, fill=chain)) +
  geom_density(alpha = 0.35, size = 1) +
  ylab(bquote("density "~beta[0])) + 
  scale_fill_manual(values = c("red", "blue")) #looks good

#density plot beta1 for both chains
chains_b1 <- as.data.frame(chains[, c(2, 7, 11)])
densb1 <- ggplot(melt(chains_b1, id.var='iteration', variable.name='chain'), aes(x=value, fill=chain)) +
  geom_density(alpha = 0.35, size = 1) +
  ylab(bquote("density "~beta[1])) + 
  scale_fill_manual(values = c("red", "blue")) #looks good

#density plot beta2 for both chains
chains_b2 <- as.data.frame(chains[, c(3, 8, 11)])
densb2 <- ggplot(melt(chains_b2, id.var='iteration', variable.name='chain'), aes(x=value, fill=chain)) +
  geom_density(alpha = 0.35, size = 1) +
  ylab(bquote("density "~beta[2])) + 
  scale_fill_manual(values = c("red", "blue")) #looks good

#density plot s2 for both chains
chains_s2 <- as.data.frame(chains[, c(5, 10, 11)])
denss2 <- ggplot(melt(chains_s2, id.var='iteration', variable.name='chain'), aes(x=value, fill=chain)) +
  geom_density(alpha = 0.35, size = 1) +
  ylab(bquote("density "~sigma^2)) + 
  scale_fill_manual(values = c("red", "blue")) #looks good

## Running mean plots ##
#make function that calculates the running means for me
runmean <- function(x, n = (n.iters - burnin)) {
  meanx <- numeric(n) #storage
  for (k in 1:n) {
    meanx[k] <- mean(x[1:k]) #mean, add one iteration at a time
  }
  return(meanx)
}

#calculate the running means and put them in a dataframe
runningmean  <- apply(chains[1:10], 2, runmean) 
runningmeans <- data.frame(iteration = 1:nrow(chains), runningmean)

#running mean plot beta0 for both chains
runb0 <- ggplot(runningmeans, aes(x = iteration)) +
  geom_line(aes(y = b0_ch1), color = "red") + 
  geom_line(aes(y = b0_ch2), color = "blue") +
  ggtitle(bquote(bold("Running mean plot of ")~beta[0])) + ylab(bquote("running mean of "~beta[0])) +
  theme(plot.title = element_text(hjust = 0.5)) #looks good

#running mean plot beta1 for both chains
runb1 <- ggplot(runningmeans, aes(x = iteration)) +
  geom_line(aes(y = b1_ch1), color = "red") + 
  geom_line(aes(y = b1_ch2), color = "blue") +
  ggtitle(bquote(bold("Running mean plot of ")~beta[1])) + ylab(bquote("running mean of "~beta[1])) +
  theme(plot.title = element_text(hjust = 0.5)) #looks okay

#running mean plot beta2 for both chains
runb2 <- ggplot(runningmeans, aes(x = iteration)) +
  geom_line(aes(y = b2_ch1), color = "red") + 
  geom_line(aes(y = b2_ch2), color = "blue") +
  ggtitle(bquote(bold("Running mean plot of ")~beta[2])) + ylab(bquote("running mean of "~beta[2])) +
  theme(plot.title = element_text(hjust = 0.5)) #looks okay

#running mean plot beta0 for both chains
runs2 <- ggplot(runningmeans, aes(x = iteration)) +
  geom_line(aes(y = s2_ch1), color = "red") + 
  geom_line(aes(y = s2_ch2), color = "blue") +
  ggtitle(bquote(bold("Running mean plot of ")~sigma^2)) + ylab(bquote("running mean of "~sigma^2)) +
  theme(plot.title = element_text(hjust = 0.5)) #looks good

## Density and running mean plots in one figure ##
plot_grid(densb0, runb0, densb1, runb1, densb2, runb2, denss2, runs2, 
          align = "h", nrow = 4, rel_widths = c(1/2, 1/2, 1/2, 1/2), rel_heights = c(1/4, 1/4, 1/4, 1/4)) #tip: click zoom

## Gelman-Rubin statistics ##
n.iters_burnin <- n.iters - burnin #prep, this is needed in all calculations
B <- n.iters_burnin * ((est_ch1[,1] - est[,1])^2 + (est_ch2[,1] - est[,1])^2) #between chains variance
s2j_ch1 <- (1/(n.iters_burnin - 1)) * sum((chain1 - est_ch1[,1])^2) #variance within chain 1
s2j_ch2 <- (1/(n.iters_burnin - 1)) * sum((chain2 - est_ch2[,1])^2) #variance within chain 2
W <- (s2j_ch1 + s2j_ch2) / 2 #pooled
R <- ((((n.iters_burnin - 1) / n.iters_burnin) * W) + ((1 / n.iters_burnin) * B)) / W #all just short of 1, good

## MC errors ##
output[1,]$sd / sqrt(n.iters) #intercept, good below 5%
output[2,]$sd / sqrt(n.iters) #beta1, good below 5%
output[3,]$sd / sqrt(n.iters) #beta2, good below 5%
output[4,]$sd / sqrt(n.iters) #variance, good below 5%

### POSTERIOR PREDICTIVE CHECK ###
#normality of the residuals -> mean = median = mode, discrepancy measure is difference mean and mode
obs_mode.mean.dif <- rep(0, nrow(rchains)) #storage 
sim_mode.mean.dif <- rep(0, nrow(rchains)) #storage 
#function to compute mode
mode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
#sample replicated datasets and compute discrepancy measure 
set.seed(10) #reproducibility
for (i in 1:nrow(rchains)){
  y <- rchains[i,1] + X1 * rchains[i,2] + X2 * rchains[i,3] #predicted y values
  residuals <- Y - y #observed residuals
  obs_mode.mean.dif[i] <- mean(residuals) - mode(residuals) #difference mean and mode observed
  ppp.y <- rnorm(n, rchains[i,1] + X1 * rchains[i,2] + X2 * rchains[i,3], sqrt(rchains[i,5])) #generate datasets
  ppp.residuals <- ppp.y - (rchains[i,1] + X1 * rchains[i,2] + X2 * rchains[i,3]) #simulated residuals 
  sim_mode.mean.dif[i] <- mean(ppp.residuals) - mode(ppp.residuals) #difference mean and mode simulated
}
ppp.value <- mean(ifelse(abs(sim_mode.mean.dif) > abs(obs_mode.mean.dif), 1, 0)) #compare absolute values!

### MODEL COMPARISON WIT THE DIC ###
## Run competing model ##
BostonHousing$rm2 <- BostonHousing$rm.cent^2 #make variable where rm is squared
X3 <- as.matrix(BostonHousing$rm2) #prep, assign predictor 3
set.seed(10) #reproducibility
#and run the function again
chain1B <- Gibbs_withMH(Y, X1, X2, X3, n, #no need for the default parameters
                        beta0 = 5, beta1 = 0.1, beta2 = 0.1, beta3 = 0.1, s2 = 100)
chain2B <- Gibbs_withMH(Y, X1, X2, X3, n, 
                        beta0 = 50, beta1 = 2, beta2 = 2, beta3 = 2, s2 = 50) #different initial values
#I looked at the results per chain and their acceptance ratios, but those are not mentioned in the report, so deleted

#get output needed to get results and compute DIC
rchainsB <- data.frame(rbind(chain1B, chain2B), 
                       chain = c(rep(1, (n.iters - burnin)), rep(2, (n.iters - burnin))),
                       iteration = rep(1:(n.iters - burnin), 2))
estB <- data.frame(est = apply(rchainsB[, 1:5], 2, mean)) #dataframe single column with posterior means to use

## Table 4 - Results Model 2 ##
bB  <- c(estB[2,], estB[3,], estB[4,]) #prep for standardized coefficients
syB <- sd(BostonHousing$medv) 
sxB <- c(sd(BostonHousing$rm.cent), sd(BostonHousing$lstat.cent), sd(BostonHousing$rm2)) 
outputB <- data.frame(mean = round(apply(rchainsB[, c(1:5)], 2, mean), 3), #posterior mean
                      beta = c("", round(bB * (sxB/syB), 3), ""), #standardized coefficient
                      sd = round(c(apply(rchainsB[, c(1:5)], 2, sd)), 3), #posterior sd
                      naive.se = round(apply(rchainsB[, c(1:5)], 2, sd)/sqrt(n.iters), 3), 
                      cred.low = round(apply(rchainsB[, c(1:5)], 2, quantile, probs = 0.025), 3), #low credible interval
                      cred.high = round(apply(rchainsB[, c(1:5)], 2, quantile, probs = 0.975), 3)) #high credible interval

## DIC Model 1 ##
y <- est[1,] + X1 * est[2,] + X2 * est[3,] #predicted y values
deviance <- -2 * sum(dnorm(Y, mean = y, sd = sqrt(est[4,]), log = T)) #Dhat
mean_deviance <- rep(0, nrow(rchains)) #storage
#calculate y, residuals and the deviance for each row in rchains (so all iterations in both chains combined)
for (i in 1:nrow(rchains)) {
  y <- rchains[i,1] + X1 * rchains[i,2] + X2 * rchains[i,3]
  mean_deviance[i] <- -2 * sum(dnorm(Y, mean = y, sd = sqrt(rchains[i,5]), log = T)) #mean is Dbar
}
penalty <- mean(mean_deviance) - deviance #Dbar - Dhat
DIC <- deviance + 2*penalty #3173.55

## DIC Model 2 ##
yB <- estB[1,] + X1 * estB[2,] + X2 * estB[3,] + X3 * estB[4,] #predicted y values
devianceB <- -2 * sum(dnorm(Y, mean = yB, sd = sqrt(estB[5,]), log = T)) #Dhat
mean_devianceB <- rep(0, nrow(rchains)) #storage
#calculate y, residuals and the deviance for each row in rchains (so all iterations in both chains combined)
for (i in 1:nrow(rchainsB)) {
  yB <- rchainsB[i,1] + X1 * rchainsB[i,2] + X2 * rchainsB[i,3] + X3 * rchainsB[i,4]
  mean_devianceB[i] <- -2 * sum(dnorm(Y, mean = yB, sd = sqrt(rchainsB[i,5]), log = T)) #mean is Dbar
}
penaltyB <- mean(mean_devianceB) - devianceB #Dbar - Dhat
DICB <- devianceB + 2*penaltyB #3024.837

### BAYES FACTOR  I - For hypothesis testing ###
#change into hstat to make the sign the same
BostonHousing$hstat <- 100 - BostonHousing$lstat #invert
BostonHousing$hstat.cent <- scale(BostonHousing$hstat, scale = F) #and center
m1 <- lm(medv ~ rm.cent + hstat.cent, data = BostonHousing) #run model
set.seed(10) #reproducibility
#test hypothesis, effect lstat bigger than rm
resultsh1 <- bain(m1, "rm.cent = hstat.cent; hstat.cent > rm.cent", standardize = T) #run bain
summary(resultsh1, ci = 0.95) #look at the results
resultsh1$BFmatrix #H1 more likely than H2

#again, to ensure stability 
set.seed(1) #reproducibility
resultsh1b <- bain(m1, "rm.cent = hstat.cent; hstat.cent > rm.cent", standardize = T) #test hypothesis, effect lstat bigger than rm
resultsh1b$BFmatrix #yes, the same, so that is good!

### BAYES FACTOR  II - For model comparison ###
M1 <- lmBF(medv ~ rm.cent + lstat.cent, data = BostonHousing) #model 1
M2 <- lmBF(medv ~ rm.cent + rm2 + lstat.cent, data = BostonHousing) #model 2
BF <- M2 / M1 #divide to compute BF

### COMPARISON OF FREQUENTIST AND BAYESIAN APPROACHES ###
## Small sensitivity analysis ##
#try some different prior means and variances for rm
set.seed(10) #reproducibility
chain1C <- Gibbs_withMH(Y, X1, X2, X3 = 0, n, #no need for the other default parameters
                        mu1= 2, s21 = 30,
                        beta0 = 5, beta1 = 0.1, beta2 = 0.1, beta3 = 0, s2 = 100)                       
chain2C <- Gibbs_withMH(Y, X1, X2, X3 = 0, n,
                        mu1 = 2, s21 = 30, 
                        beta0 = 50, beta1 = 2, beta2 = 2, beta3 = 0, s2 = 50) #different initial values
#look at the results 
est_ch1C <- data.frame(est = apply(chain1C[, c(1:3, 5)], 2, mean), 
                       cred.low = apply(chain1C[, c(1:3, 5)], 2, quantile, probs = 0.025), 
                       cred.high = apply(chain1C[, c(1:3, 5)], 2, quantile, probs = 0.975)) #did not affect results
est_ch2C <- data.frame(est = apply(chain2C[, c(1:3, 5)], 2, mean), 
                       cred.low = apply(chain2C[, c(1:3, 5)], 2, quantile, probs = 0.025), 
                       cred.high = apply(chain2C[, c(1:3, 5)], 2, quantile, probs = 0.975)) #did not affect results

#same, but a bigger value 
set.seed(10) #reproducibility
chain1D <- Gibbs_withMH(Y, X1, X2, X3 = 0, n, #no need for the other default parameters
                        mu1= 20, s21 = 30, 
                        beta0 = 5, beta1 = 0.1, beta2 = 0.1, beta3 = 0, s2 = 100)                       
chain2D <- Gibbs_withMH(Y, X1, X2, X3 = 0, n,
                        mu1 = 20, s21 = 30, 
                        beta0 = 50, beta1 = 2, beta2 = 2, beta3 = 0, s2 = 50) #different initial values
#look at the results 
est_ch1D <- data.frame(est = apply(chain1D[, c(1:3, 5)], 2, mean), 
                       cred.low = apply(chain1D[, c(1:3, 5)], 2, quantile, probs = 0.025), 
                       cred.high = apply(chain1D[, c(1:3, 5)], 2, quantile, probs = 0.975)) #did not affect results
est_ch2D <- data.frame(est = apply(chain2D[, c(1:3, 5)], 2, mean), 
                       cred.low = apply(chain2D[, c(1:3, 5)], 2, quantile, probs = 0.025), 
                       cred.high = apply(chain2D[, c(1:3, 5)], 2, quantile, probs = 0.975)) #did not affect results

#then try a larger variance
set.seed(10) #reproducibility
chain1E <- Gibbs_withMH(Y, X1, X2, X3 = 0, n, #no need for the other default parameters
                        mu1= 2, s21 = 50, 
                        beta0 = 5, beta1 = 0.1, beta2 = 0.1, beta3 = 0, s2 = 100)                       
chain2E <- Gibbs_withMH(Y, X1, X2, X3 = 0, n,
                        mu1 = 2, s21 = 50, 
                        beta0 = 50, beta1 = 2, beta2 = 2, beta3 = 0, s2 = 50) #different initial values
#look at the results 
est_ch1E <- data.frame(est = apply(chain1E[, c(1:3, 5)], 2, mean), 
                       cred.low = apply(chain1E[, c(1:3, 5)], 2, quantile, probs = 0.025), 
                       cred.high = apply(chain1E[, c(1:3, 5)], 2, quantile, probs = 0.975)) #did not affect results
est_ch2E <- data.frame(est = apply(chain2E[, c(1:3, 5)], 2, mean), 
                       cred.low = apply(chain2E[, c(1:3, 5)], 2, quantile, probs = 0.025), 
                       cred.high = apply(chain2E[, c(1:3, 5)], 2, quantile, probs = 0.975)) #did not affect results

#same, but a bigger value 
set.seed(10) #reproducibility
chain1F <- Gibbs_withMH(Y, X1, X2, X3 = 0, n, #no need for the other default parameters
                        mu1= 20, s21 = 30, 
                        beta0 = 5, beta1 = 0.1, beta2 = 0.1, beta3 = 0, s2 = 100)                       
chain2F <- Gibbs_withMH(Y, X1, X2, X3 = 0, n,
                        mu1 = 20, s21 = 30, 
                        beta0 = 50, beta1 = 2, beta2 = 2, beta3 = 0, s2 = 50) #different initial values
#look at the results 
est_ch1F <- data.frame(est = apply(chain1F[, c(1:3, 5)], 2, mean), 
                       cred.low = apply(chain1F[, c(1:3, 5)], 2, quantile, probs = 0.025), 
                       cred.high = apply(chain1F[, c(1:3, 5)], 2, quantile, probs = 0.975)) #did not affect results
est_ch2F <- data.frame(est = apply(chain2F[, c(1:3, 5)], 2, mean), 
                       cred.low = apply(chain2F[, c(1:3, 5)], 2, quantile, probs = 0.025), 
                       cred.high = apply(chain2F[, c(1:3, 5)], 2, quantile, probs = 0.975))

#thank you for making it this far