# run_all.R
# Set working directory to REPLICATION PACKAGE/ before running.
# Toggle sections as needed.

run_main     <- TRUE
run_appC     <- TRUE
run_appD     <- TRUE
run_ols_diag <- TRUE


source("code/vrpoprob.R")
source("code/config.R")
source("code/load_data.R")
source("code/calc_weights.R")

if (run_main) {
  source("code/calc_results.R")     # --> results/results.RData, summary_table.csv
} else {
  load("results/results.RData")
}

source("code/plot_results.R")


# Appendix C: Monte Carlo simulation (~1.5 hrs on 9 cores)
if (run_appC) {
  source("code/simulation_appC.R")   # --> results/simulation.RData + tables
}

# Appendix D: web robustness with IW_ONLINE exclusion restriction
if (run_appD) {
  source("code/calc_results_appD.R") # --> results/results_web.RData
  source("code/plot_results_appD.R")
}

# Appendix D: OLS diagnostics for IW_ONLINE
if (run_ols_diag) {
  source("code/table_ols_diagnostics.R") # --> results/table_ols_diagnostics.csv
}


