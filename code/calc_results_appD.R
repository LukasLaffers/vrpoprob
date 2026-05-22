# web-only robustness: V241620 (online) as exclusion restriction in Z
# produces two result objects
#   res_array_web_noZ    -- xdata = zdata = demographics
#   res_array_web_withZ  -- zdata = demographics + online

tic("Web robustness")

source("code/calc_weights_appD.R")  # dff3_web, Xpop_web, Wpop_web, Zpop_withZ, Wpop_withZ


res_array_web_noZ   <- list()
res_array_web_withZ <- list()

for (miss_p in miss_props) {

  res_array_web_noZ[[as.character(miss_p)]]   <- list()
  res_array_web_withZ[[as.character(miss_p)]] <- list()
  message("Running missing proportion: ", round(100 * miss_p), "%")

  for (iOutcome in seq_along(outcomes)) {

    message("  Calculating outcome: ", outcomes[iOutcome])

    dff4         <- dff3_web[dff3_web[[outcomes[iOutcome]]] > 0, ]
    outcome_here <- dff4[[outcomes[iOutcome]]]
    NN           <- nrow(dff4)
    Nmiss_here   <- round(NN * miss_p / (1 - miss_p))

    xdata <- data.matrix(dff4 %>% dplyr::select(married, black, gender_spouse, education))
    zdata <- data.matrix(dff4 %>% dplyr::select(married, black, gender_spouse, education, online))

    res_array_web_noZ[[as.character(miss_p)]][[iOutcome]] <- vrpoprob_estim(
      ydata = outcome_here,
      rdata = dff4$int_rating,
      xdata = xdata,
      zdata = xdata,
      Nmiss = Nmiss_here,
      Wpop = Wpop_web,
      Xpop  = data.matrix(Xpop_web),
      Zpop  = data.matrix(Xpop_web)
    )

    res_array_web_withZ[[as.character(miss_p)]][[iOutcome]] <- vrpoprob_estim(
      ydata = outcome_here,
      rdata = dff4$int_rating,
      xdata = xdata,
      zdata = zdata,
      Nmiss = Nmiss_here,
      Wpop = Wpop_withZ,
      Xpop  = data.matrix(Zpop_withZ)[, 1:4, drop = FALSE],
      Zpop  = data.matrix(Zpop_withZ)
    )
  }
}

save(res_array_web_noZ, res_array_web_withZ,
     dff3_web, Xpop_web, Wpop_web, Zpop_withZ, Wpop_withZ,
     file = "results/results_web.RData")

toc()
