#################
## Convergence ##
#################

## Load workspace ##
file("Workspaces/02 Estimations.RData")

## Call functions needed in this section ##
source("Functions/Convergence_Diagnostics.R") 

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
#this uses a function in the source("Functions/Convergence_Diagnostics.R") file
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
#this uses a function in the source("Functions/Convergence_Diagnostics.R") file
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

## Gelman-Rubin statistics ##
n.iters_burnin <- n.iters - burnin #prep, this is needed in all calculations
B <- n.iters_burnin * ((est_ch1[,1] - est[,1])^2 + (est_ch2[,1] - est[,1])^2) #between chains variance
s2j_ch1 <- (1/(n.iters_burnin - 1)) * sum((chain1 - est_ch1[,1])^2) #variance within chain 1
s2j_ch2 <- (1/(n.iters_burnin - 1)) * sum((chain2 - est_ch2[,1])^2) #variance within chain 2
W <- (s2j_ch1 + s2j_ch2) / 2 #pooled
R <- ((((n.iters_burnin - 1) / n.iters_burnin) * W) + ((1 / n.iters_burnin) * B)) / W #all just short of 1, good

## MC errors ##
#this uses a function in the source("Functions/Convergence_Diagnostics.R") file
mc_error(sd(rchains[, 1]), n.iters) #intercept, good below 5%
mc_error(sd(rchains[, 2]), n.iters) #beta1, good below 5%
mc_error(sd(rchains[, 3]), n.iters) #beta2, good below 5%
mc_error(sd(rchains[, 5]), n.iters) #variance, good below 5%

## Save workspace ##
save.image("Workspaces/03 Convergence.RData")
