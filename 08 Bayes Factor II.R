############################################
## Bayes Factor II - For model comparison ##
############################################

## Load workspace ##
load("Workspaces/07 Model comparison with the DIC.RData")

M1 <- lmBF(medv ~ rm.cent + lstat.cent, data = BostonHousing) #model 1
M2 <- lmBF(medv ~ rm.cent + rm2 + lstat.cent, data = BostonHousing) #model 2
BF <- M2 / M1 #divide to compute BF

## Save workspace ##
save.image("Workspaces/08 Bayes Factor II.RData")

