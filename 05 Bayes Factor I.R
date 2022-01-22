##############################################
## Bayes Factor  I - For hypothesis testing ##
##############################################

## Load workspace ##
load("Workspaces/04 Interpretation of estimates and intervals.RData")

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

## Save workspace ##
save.image("Workspaces/05 Bayes Factor I.RData")
