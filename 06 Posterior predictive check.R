################################
## Posterior predictive check ##
################################

## Load workspace ##
load("Workspaces/05 Bayes Factor I.RData")

## Call functions needed in this section ##
source("Functions/Posterior_Predictive_Check.R")
set.seed(10) #reproducibility
ppp(rchains, Y, X1, X2, n) #0.37

## Save workspace ##
save.image("Workspaces/06 Posterior predtice check.RData")
