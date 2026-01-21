calculate_results <- TRUE

source("code/load_data.R")
source("code/calc_weights.R")

if (calculate_results) {
  source("code/calc_results.R")
} else {
  load("results/results.RData")
}

source("code/plot_results.R")
