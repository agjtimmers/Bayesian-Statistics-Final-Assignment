################
## Estimation ##
################

## Load workspace ##
load("Workspaces/01 Introduction dataset.RData")

## Call function Gibbs sampler with Metropolis-Hastings step ##
source("Functions/Gibbs_withMH.R") 

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

## Save workspace ##
save.image("Workspaces/02 Estimation.RData")
