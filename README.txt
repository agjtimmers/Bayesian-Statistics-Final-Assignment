## Bayesian-Statistics-Final-Assignment ##

This repository contains the files for the final assignment of the course Bayesian Statistics at Utrecht Univesity. To get the same results as in the Final Assignment Report file, download this repository and run the R code in the R files in the order in which they are numbered (which makes use of the functions in the Functions folder). This repository also serves as a research archive for the course Markup Languages at Utrecht Univesity.

## Contents of the repository ##
   - README (.txt, you are here)
   - Final assignment report (.pdf)
   - 01 Introduction dataset (.R, includes the dataset)
   - 02 Estimation (.R
   - 03 Convergence (.R)
   - 04 Interpretation of estimates and intervals (.R)
   - 05 Bayes Factor I (.R)
   - 06 Posterior predictive check (.R)
   - 07 Model comparison with the DIC (.R)
   - 08 Bayes Factor II (.R)
   - 09 Comparison of frequentist and Bayesian approaches (.R)
   - Workspaces folder with saved workspaces of all the above R files (01 - 09)
   - Functions folder:
      - Gibbs_withMH (.R)
      - Convergence_Diagnostics (.R)
      - Posterior_Predictive_Check (.R)

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
01 Introduction dataset 
   - Table 1 - Descriptive statistics 			line 12
   - Table 2 - Bivariate associations			line 15

02 Estimation						
   - Function with vague priors in default		Gibbs_withMH.R
   - Preparation input					lines 13 - 19
     - Center variables					lines 13 and 14
     - Assign input function				lines 16 - 19
   - Run two chains					lines 22 - 29
   - Calculate acceptance ratios 			lines 32 - 34
   - Preparation convergence plots and results		lines 38 - 46

03 Convergence
   - Trace plots					lines 13 - 38
   - Autocorrelation plots				lines 43 - 76
   - Density plots					lines 80 - 105
   - Running mean plots					lines 110 - 139
   - Gelman-Rubin statistics				lines 142 - 147
   - MC errors						lines 151 - 154 

04 Interpretation of estimates and intervals		
   - Table 3 - Results Model 1				

05 Bayes Factor I						

06 Posterior predictive check 					

07 Model comparison with the DIC			
  - Run competing model					lines 9 - 23
  - Table 4 - Results Model 2				lines 26 - 34
  - DIC Model 1						lines 37 - 46
  - DIC Model 2						lines 49 - 58

08 Bayes Factor II							

09 Comparison of frequentist and Bayesian approaches

