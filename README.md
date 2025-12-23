# vrpoprob
Correcting for Nonignorable Nonresponse Bias in Ordinal Observational Survey Data
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

`run_all.R` - runs all the scripts in the replication

* `data` folder:
  * `anes_timeseries_2024_csv_20250219.csv` dataset from ANES 2025
  * `anes_timeseries_2024_userguidecodebook_20250219.pdf` codebook from ANES 2025

* `code` folder:
  * `load_data.R` loads the data
  * `calc_weights.R` calculate sampling weights
  * `vrpoprob.R` contains the functions necessary for estimating the model
  * `calc_results.R` perform all the calculations, replicates all the
  * `plot_results.R` creates all the figures in the paper and saves them to `plots` folder

* `plots` folder:
  * `XXX.pdf` Figure 1
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
