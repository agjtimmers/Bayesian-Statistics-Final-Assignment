#######################################################
## Comparison of frequentist and Bayesian approaches ##
#######################################################

##Load workspace ##
load("Workspaces/08 Bayes Factor II.RData")

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

## Save workspace ##
save.image("Workspaces/09 Comparison of frequentist and Bayesian approaches.RData")
