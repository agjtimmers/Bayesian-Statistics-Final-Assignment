## Bayesian-Statistics-Final-Assignment ##

This repository contains the files for the final assignment of the course Bayesian Statistics at Utrecht Univesity. To get the same results as in the Final Assignment Report file, download this repository and run the R code in the Final Assignment Dataset and Code file (which makes use of the functions in the Functions folder). This repository also serves as a research archive for the course Markup Languages at Utrecht Univesity.

## Contents of the repository ##
   - README (.txt, you are here)
   - Final assignment report (.pdf)
   - Final Assignment Dataset and Code (.R, includes the dataset)
   - Gibbs_withMH (.R, in Functions folder)
   - Convergence_Diagnostics (.R, in Functions folder)
   - Posterior_Predictive_Check (.R, in Functions folder)

## Final assignment report ##
1. Research question 					page 1

2. Introduction dataset 				page 1
   - Table 1 - Descriptive statistics 			page 1
   - Table 2 - Bivariate associations 			page 1

3. Estimation						page 2
   - Gibbs sampler					page 2
   - The Metropolis-Hastings algorithm			page 2

4. Convergence						page 2 - 3

5. Interpretation of estimates and intervals		page 3
   - Table 3 - Results Model 1				page 3

6. Bayes Factor I					page 3 - 4

7. Posterior predictive check				page 4

8. Model comparison with the DIC			page 4 - 5
   - Table 4 - Results Model 2				page 5

9. Bayes Factor II					page 5

10. Comparison of frequentist and Bayesian approaches	page 5 - 6

11. References						page 6

## BostonHousing dataset ##
This dataset contains 506 cases and 14 variables, of which 3 were used:
   - medv: the median value of the housing prices in the census tracts in the Boston Standard Metropolitan Statistical Area in 1970.
   - rm: the average number of rooms per house in the tracts.
   - lstat: the percentage of the population of the tracts that is of a lower socioeconomic status.
The dataset can be found in R and is loaded in the attached R Code (see below). 

## R Code ##
1. Preparation
   - Load libraries 					line 8
   - Load dataset 					line 10

2. Introduction dataset 
   - Table 1 - Descriptive statistics 			line 14
   - Table 2 - Bivariate associations			line 17

3. Estimation						
   - Function with vague priors in default		Gibbs_withMH.R
   - Preparation input					lines 25 - 31
     - Center variables					lines 25 and 26
     - Assign input function				lines 28 - 31
   - Run two chains					lines 34 - 41
   - Calculate acceptance ratios 			lines 44 - 46
   - Preparation convergence plots and results		lines 50 - 58

4. Convergence
   - Trace plots					lines 78 - 103
   - Autocorrelation plots				lines 108 - 141
   - Trace and autocorrelation plots in one figure 	lines 144 and 145
   - Density plots					lines 149 - 174
   - Running mean plots					lines 179 - 208
   - Density and running mean plots in one figure	lines 211 and 212
   - Gelman-Rubin statistics				lines 215 - 220
   - MC errors						lines 224 - 227 

5. Interpretation of estimates and intervals		
   - Table 3 - Results Model 1				lines 62 - 70

6. Bayes Factor I					lines 290 - 302	

7. Posterior predictive check 				lines 231 - 233	

8. Model comparison with the DIC			
  - Run competing model					lines 237 - 251
  - Table 4 - Results Model 2				lines 254 - 261
  - DIC Model 1						lines 265 - 274
  - DIC Model 2						lines 277 - 286

9. Bayes Factor II					lines 305 - 307		

10. Comparison of frequentist and Bayesian approaches
   - Referenced ML estimates				line 292
   - Referenced small sensitivity analysis		lines 312 - 373

