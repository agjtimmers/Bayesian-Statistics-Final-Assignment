###############################################
## Interpretation of estimates and intervals ##
###############################################

## Load workspace ##
load("Workspaces/03 Convergence.RData")

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

## Save workspace ##
save.image("Workspaces/04 Interpretation of estimates and intervals.RData")
