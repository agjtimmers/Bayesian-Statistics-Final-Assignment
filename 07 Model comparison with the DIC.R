###################################
## Model comparison with the DIC ##
###################################

## Load workspace ##
load("Workspaces/06 Posterior predtice check.RData")

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

## Save workspace ##
save.image("Workspaces/07 Model comparison with the DIC.RData")
