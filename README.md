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

#######################################################################################


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
  * `XXX.pdf` Distribution of the response variable *Rating of interview* Figure 1
  * `response_life_all.pdf`, `life_all.pdf`, `response_comp_life_all.pdf` Figure 2
  * `response_economy_all.pdf`, `economy_all.pdf`, `response_comp_economy_all.pdf` Figure 3
  * `response_unemployment_all.pdf`, `unemployment_all.pdf`, `response_comp_unemployment_all.pdf` Figure 4 (Question 3)
  * `response_media_all.pdf`, `media_all.pdf`, `response_comp_media_all.pdf` Figure 5 (Question 4)
  * `response_votes_all.pdf`, `votes_all.pdf`, `response_comp_votes_all.pdf` Figure 6 (Question 5)
  * `response_religion_all.pdf`, `religion_all.pdf`, `response_comp_religion_all.pdf` Figure 7 (Question 6)
  * `response_abortions_all.pdf`, `abortions_all.pdf`, `response_comp_abortions_all.pdf` Figure 8 (Question 7)
  * `response_death_all.pdf`, `death_all.pdf`, `response_comp_death_all.pdf` Figure 9 (Question 8)


#######################################################################################


#### SOFTWARE DEPENDENCIES

All scripts were run under the following environment:

R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] latex2exp_0.9.6     RColorBrewer_1.1-3  scales_1.3.0        ggtext_0.1.2       
 [5] tictoc_1.2.1        lubridate_1.9.3     forcats_1.0.0       stringr_1.5.1      
 [9] dplyr_1.1.4         purrr_1.0.2         readr_2.1.5         tidyr_1.3.1        
[13] tibble_3.2.1        ggplot2_3.5.1       tidyverse_2.0.0     maxLik_1.5-2.1     
[17] miscTools_0.6-28    mvtnorm_1.3-3       numDeriv_2016.8-1.1

loaded via a namespace (and not attached):
 [1] sandwich_3.1-1    utf8_1.2.4        generics_0.1.3    xml2_1.3.6        stringi_1.8.4    
 [6] lattice_0.22-6    hms_1.1.3         digest_0.6.37     magrittr_2.0.3    grid_4.4.1       
[11] timechange_0.3.0  fansi_1.0.6       cli_3.6.3         rlang_1.1.4       crayon_1.5.3     
[16] bit64_4.5.2       munsell_0.5.1     withr_3.0.2       tools_4.4.1       parallel_4.4.1   
[21] tzdb_0.4.0        colorspace_2.1-1  vctrs_0.6.5       R6_2.5.1          zoo_1.8-14       
[26] lifecycle_1.0.4   bit_4.5.0         vroom_1.6.5       pkgconfig_2.0.3   pillar_1.9.0     
[31] gtable_0.3.6      glue_1.8.0        Rcpp_1.0.13-1     tidyselect_1.2.1  rstudioapi_0.17.1
[36] compiler_4.4.1    gridtext_0.1.5  


#######################################################################################


#### REFERENCES

* American National Election Studies. 2025. ANES 2024 Time Series Study Preliminary Release: Combined Pre-Election and Post-Election Data [dataset and documentation]. April 30, 2025 
