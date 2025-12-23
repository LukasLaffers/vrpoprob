# vrpoprob
**V**ariable **r**esponse **p**ropensity model with **o**rdered **prob**it
***************************************************************************************************************
*** README                                                                         			    ***
*** Replication file for: Correcting for Nonignorable Nonresponse Bias in Ordinal Observational Survey Data ***
***         							                    			    ***
*** Version: 0.6 (December 23, 2025)                                                        			    ***
***************************************************************************************************************


#### NOTES

Welcome to the replication package for "Correcting for Nonignorable Nonresponse Bias in Ordinal Observational Survey Data", encompassing the data, code, and documentation. This README provides detailed information on the  dataset, software environment, and reproducibility instructions. This replication repository generates the analysis and figures used in the paper and supplemental materials. The script `run_all.R` runs all scripts in the correct order. 

The raw data is from ANES (2025), which is available on the [www.electionstudies.org](https://electionstudies.org/data-center/2024-time-series-study/).


#### TABLE OF CONTENTS

* `run_all.R` - runs all the scripts in the replication, if `calculate_results` is set to `FALSE`, then instead of running `calc_results.R`, results are loaded from `results.RData`
* `README.md` - this file

* `data` folder:
  * `anes_timeseries_2024_csv_20250219.csv` dataset from ANES 2025
  * `anes_timeseries_2024_userguidecodebook_20250219.pdf` codebook from ANES 2025

* `code` folder:
  * `load_data.R` loads the data
  * `calc_weights.R` calculate sampling weights
  * `vrpoprob.R` contains the functions necessary for estimating the model
    * **Main functions**
      * `vrpoprob_estim(ydata, rdata, xdata, zdata, Nmiss, WXpop, Xpop, WZpop, Zpop)`  Main estimation routine. Returns estimated parameters (`alpha`, `beta`, `lambda`, `theta`, `rho`), standard errors, and estimated population proportions (`pphat`, `pphat_nonresp`, `pphat_resp`).
      * `vrpoprob_simdata(N, alpha, beta, lambda, theta, rho, WXpop, Xpop)`   Simulates dataset of size `N` including nonresponses. Returns `ydata`, `rdata`, `xdata`, `zdata`, `Nmiss`
    * **Internal / helper functions**
      * `vrpoprob_pack(alpha, beta, lambda, theta, rho)` – packs parameters into a single vector.  
      * `vrpoprob_unpack(xi, J, K, M, R)` – unpacks parameter vector into `alpha`, `beta`, `lambda`, `theta`, `rho`.  
      * `vrpoprob_eval_yrlogp(y, r, ystarhat, rstarhat, lambda, theta, rho)` – evaluates log-probability of `(y, r)` conditional on `(x, z)`.  
      * `vrpoprob_eval_npunc(beta, WZpop, Zpop)` – evaluates unconditional log-probability of nonresponse.  
      * `vrpoprob_loglik(xi, ydata, rdata, xdata, zdata, Nmiss, WZpop, Zpop)` – computes the full log-likelihood.  
      * `vrpoprob_sim1obs(alpha, beta, lambda, theta, rho, WXpop, Xpop)` – simulates a single observation.  
      * `vrpoprob_xi_to_pphat(xi, WXpop, Xpop, J, K, M, R)` – computes population outcome proportions from parameters.  
      * `vrpoprob_delta_se(f, x, V)` – delta-method computation of standard errors.  
      * `vrpoprob_xi_to_pphat_resp_nonresp(xi, WXpop, Xpop, Zpop, J, K, M, R)` – computes population proportions separately for respondents and nonrespondents.
  * `calc_results.R` perform all the calculations, replicates all the results used in the paper. It takes about ~6hrs on M1PRO 2021 16GB RAM laptop.
  * `plot_results.R` creates all the figures in the paper and saves them to `plots` folder

* `results` folder:
  * `results.RData` file with results of the replication

* `plots` folder:
  * `rating.pdf` Distribution of the response variable *Rating of interview* Figure 1
  * `serious.pdf` Distribution of the response variable *How often you take survey seriously* Figure 2
  * `response_life_all.pdf`, `life_all.pdf`, `response_comp_life_all.pdf` Figure 3
  * `response_economy_all.pdf`, `economy_all.pdf`, `response_comp_economy_all.pdf` Figure 4
  * `response_unemployment_all.pdf`, `unemployment_all.pdf`, `response_comp_unemployment_all.pdf` Figure 5 (Question 3)
  * `response_media_all.pdf`, `media_all.pdf`, `response_comp_media_all.pdf` Figure 6 (Question 4)
  * `response_votes_all.pdf`, `votes_all.pdf`, `response_comp_votes_all.pdf` Figure 7 (Question 5)
  * `response_religion_all.pdf`, `religion_all.pdf`, `response_comp_religion_all.pdf` Figure 8 (Question 6)
  * `response_abortions_all.pdf`, `abortions_all.pdf`, `response_comp_abortions_all.pdf` Figure 9 (Question 7)
  * `response_death_all.pdf`, `death_all.pdf`, `response_comp_death_all.pdf` Figure 10 (Question 8)


#### SOFTWARE DEPENDENCIES

All scripts were run under the following environment:

R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS 26.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Bratislava
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] colorspace_2.1-1  scales_1.3.0      compiler_4.4.2    R6_2.6.1          cli_3.6.4         tools_4.4.2       glue_1.8.0        rstudioapi_0.17.1
 [9] lifecycle_1.0.4   munsell_0.5.1     rlang_1.1.5      


#### REFERENCES

* American National Election Studies. 2025. ANES 2024 Time Series Study Preliminary Release: Combined Pre-Election and Post-Election Data [dataset and documentation]. April 30, 2025 
