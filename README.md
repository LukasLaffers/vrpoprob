# vrpoprob

**V**ariable **r**esponse **p**ropensity model with **o**rdered **prob**it.

Replication package for *Correcting for Nonignorable Nonresponse Bias in Ordinal Observational Survey Data*.

Version 1.1 (May 2026).

The raw data is from the 2024 ANES Time Series Study, available at
[electionstudies.org](https://electionstudies.org/data-center/2024-time-series-study/).
`run_all.R` runs all scripts in the correct order.


## Layout

```
run_all.R                 entry point; toggles main / appendix C / appendix D
README.md                 this file
data/                     ANES 2024 csv + codebook
code/                     scripts (see below)
results/                  generated .RData and .csv (created on run)
plots/                    generated figures (created on run)
```

### `run_all.R` switches

| flag       | default | what it controls                                              |
|------------|---------|---------------------------------------------------------------|
| `run_main` | `TRUE`  | main ANES estimation. If `FALSE`, loads cached `results.RData` |
| `run_appC` | `TRUE`  | Monte Carlo simulation (Appendix C)                            |
| `run_appD` | `TRUE`  | web-only robustness with `IW_ONLINE` exclusion (Appendix D)    |

### `code/`

| file                    | purpose                                                                                  |
|-------------------------|------------------------------------------------------------------------------------------|
| `vrpoprob.R`            | the estimator (`vrpoprob_estim`) and its internal helpers                                |
| `config.R`              | shared outcome list, labels, and the nonresponse grid `{0.5, 0.65, 0.8}`                 |
| `load_data.R`           | reads the ANES csv and constructs `dff3` with renamed columns and basic filters          |
| `calc_weights.R`        | builds the `(X, W)` population grid from demographics                                    |
| `calc_results.R`        | runs `vrpoprob_estim` for all 8 outcomes × 3 nonresponse rates → `results/results.RData` and `results/summary_table.csv` |
| `plot_results.R`        | per-outcome bar plots and respondent/nonrespondent decomposition figures                 |
| `simulation_appC.R`     | Monte Carlo: 8 outcomes × 3 ρ × 3 p_miss × 500 reps; writes `simulation.RData`, `simulation_table1.csv`, `simulation_table2.csv` |
| `calc_weights_appD.R`   | population grids for the web subsample, with and without `online`                        |
| `calc_results_appD.R`   | web-only baseline (`Z = X`) and exclusion (`Z = X ∪ {online}`) → `results_web.RData`     |
| `plot_results_appD.R`   | web-robustness figures and `rho_web_comparison.csv`                                      |

### Top-level functions in `vrpoprob.R`

| function                                                                  | role                                                                                  |
|---------------------------------------------------------------------------|---------------------------------------------------------------------------------------|
| `vrpoprob_estim = function(ydata, rdata, xdata, zdata, Nmiss, Wpop, Xpop, Zpop)` | full estimation; returns `alpha`, `beta`, `lambda`, `theta`, `rho`, their SEs, `pphat`, `pphat_resp`, `pphat_nonresp`, and SEs |
| `vrpoprob_loglik(...)`                                                    | summed log-likelihood (used by `maxLik`)                                              |
| `vrpoprob_unpack(xi, J, K, M, R)`                                         | maps the optimizer's free vector to `(α, β, λ, θ, ρ)`                                  |
| `vrpoprob_xi_to_pphat(...)`                                               | population shares from parameters                                                     |
| `vrpoprob_xi_to_pphat_resp_nonresp(...)`                                  | population shares split by respondent status                                          |
| `vrpoprob_delta_se(f, x, V)`                                              | delta-method SEs via `numDeriv::jacobian`                                              |

### `results/` (generated)

| file                          | source                | content                                                  |
|-------------------------------|-----------------------|----------------------------------------------------------|
| `results.RData`               | `calc_results.R`      | `res_array[[p_miss]][[outcome]]` for all 8 × 3 fits      |
| `summary_table.csv`           | `calc_results.R`      | main-text Table 1 (raw vs corrected shares, ρ̂, 95% CI)   |
| `simulation.RData`            | `simulation_appC.R`   | per-replication Monte Carlo output                       |
| `simulation_summary.RData`    | `simulation_appC.R`   | cell-level summary statistics                            |
| `simulation_table1.csv`       | `simulation_appC.R`   | Appendix C Table 1 (naive vs ordinal VRP)                |
| `simulation_table2.csv`       | `simulation_appC.R`   | Appendix C Table 2 (binary VRP vs ordinal VRP)           |
| `results_web.RData`           | `calc_results_appD.R` | `res_array_web_noZ`, `res_array_web_withZ`               |
| `rho_web_comparison.csv`      | `plot_results_appD.R` | side-by-side ρ̂ across the two web specifications        |

### `plots/` (generated)

For each outcome `<outcome>`: `<outcome>_all.pdf` (corrected vs unadjusted), `response_<outcome>_all.pdf` (response distribution by interview rating), `response_comp_<outcome>_all.pdf` (respondents vs nonrespondents). Plus `rating.pdf`, `serious.pdf`. Web-robustness figures live in `plots/web_noZ/`, `plots/web_withZ/`, and `plots/web_rho_comparison.pdf`.


## Runtime

On an M1 Pro 2021 (16 GB RAM):

- main ANES (`calc_results.R`): ~5 minutes
- Monte Carlo (`simulation_appC.R`): ~1.5 hours, 9 cores via `mclapply`
- web robustness (`calc_results_appD.R`): ~6 minutes


## Software

```text
R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
```

Packages used (auto-installed by `load_data.R` if missing): `tictoc`, `readr`, `ggplot2`, `ggtext`, `scales`, `dplyr`, `tidyr`, `RColorBrewer`, `latex2exp`, `mnorm`, `parallel`, `maxLik`, `numDeriv`, `MASS`.


## References

American National Election Studies. 2025. *ANES 2024 Time Series Study Preliminary Release: Combined Pre-Election and Post-Election Data* (April 30, 2025).
