##########################
## Introduction dataset ##
##########################

## Preparation ##
#load the libraries
lapply(c("mlbench", "psych", "ggplot2", "cowplot", "reshape2", "bain", "BayesFactor"), library, character.only = TRUE)
#and the dataset
data("BostonHousing")

## Table 1 - Descriptive statistics ##
describe(BostonHousing[, c(14, 6, 13)])

## Table 2 - Bivariate associations ##
cor(BostonHousing[, c(14, 6, 13)])

## Save workspace ##
save.image("Workspaces/01 Introduction dataset.RData")
