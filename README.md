## Repository containg my final assignment for the course Bayesian Statistics at the UU ##

## Contents of the repository ##
   - README (.md, you are here)
   - Final assignment report (.pdf)
   - R code (.R, includes the dataset)

## Final assignment report ##
1. Research question 					                  emsp*page 1*

2. Introduction dataset 				                  -*page 1*
   - Table 1 - Descriptive statistics 			         -*page 1*
   - Table 2 - Bivariate associations 			         -*page 1*

3. Estimation						                        -*page 2*
   - Gibbs sampler					                     -*page 2*
   - The Metropolis-Hastings algorithm			         -*page 2*

4. Convergence						                        page 2 - 3

5. Interpretation of estimates and intervals	         page 3
   - Table 3 - Results Model 1				            page 3

6. Bayes Factor I					                        page 3 - 4

7. Posterior predictive check				               page 4

8. Model comparison with the DIC			               page 4 - 5
   - Table 4 - Results Model 2				            page 5

9. Bayes Factor II					                     page 5

10. Comparison of frequentist and Bayesian approaches	page 5 - 6

11. References						                        page 6

## BostonHousing dataset ##
This dataset contains 506 cases and 14 variables, of which 3 were used:
   - medv: the median value of the housing prices in the census tracts in the Boston Standard Metropolitan Statistical Area in 1970.
   - rm: the average number of rooms per house in the tracts.
   - lstat: the percentage of the population of the tracts that is of a lower socioeconomic status
The dataset can be found in R and is loaded in the attached R Code (see below). 

## R Code ##
1. Preparation
   - Load libraries 					                     line 8
   - Load dataset 					                     line 10

2. Introduction dataset 
   - Table 1 - Descriptive statistics 			         line 14
   - Table 2 - Bivariate associations			         line 17

3. Estimation						
   - Function with vague priors in default		      lines 21 - 110
   - Preparation input					                  lines 114 - 120
     - Center variables					                  lines 114 and 115
     - Assign input function				               lines 117 - 120
   - Run two chains					                     lines 123 - 130
   - Calculate acceptance ratios 			            lines 133 - 135
   - Preparation convergence plots and results		   lines 139 - 147

4. Convergence						
   - Trace plots					                        lines 164 - 189
   - Autocorrelation plots				                  lines 193 - 236
   - Trace and autocorrelation plots in one figure 	lines 239 and 240
   - Density plots					                     lines 244 - 269
   - Running mean plots					                  lines 273 - 311
   - Density and running mean plots in one figure	   lines 314 and 315
   - Gelman-Rubin statistics				               lines 318 - 323
   - MC errors						                        lines 326 - 329 

5. Interpretation of estimates and intervals		
   - Table 3 - Results Model 1				            lines 151 - 159

6. Bayes Factor I					                        lines 407 - 419	

7. Posterior predictive check 				            lines 333 - 350	

8. Model comparison with the DIC			
  - Run competing model					                  lines 354 - 368
  - Table 4 - Results Model 2				               lines 371 - 379
  - DIC Model 1						                     lines 382 - 391
  - DIC Model 2						                     lines 394 - 403

9. Bayes Factor II					                     lines 422 - 424		

10. Comparison of frequentist and Bayesian approaches
   - Referenced ML estimates				               line 409 
   - Referenced small sensitivity analysis		      lines 429 - 490

