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

    xvars <- c("married", "black", "gender_spouse", "education")
    zvars <- c("married", "black", "gender_spouse", "education", "online")

    # no-Z spec: x = z = demographics, evaluated on the demographics grid
    des_noZ <- vrpoprob_build_designs(dff4[, xvars], Xpop_web[, xvars], xvars)

    res_array_web_noZ[[as.character(miss_p)]][[iOutcome]] <- vrpoprob_estim(
      ydata = outcome_here,
      rdata = dff4$int_rating,
      xdata = des_noZ$ind,
      zdata = des_noZ$ind,
      Nmiss = Nmiss_here,
      Wpop = Wpop_web,
      Xpop  = des_noZ$pop,
      Zpop  = des_noZ$pop
    )

    # with-Z spec: z adds the 'online' exclusion restriction; both grids are
    # evaluated on the finer demographics x online grid (Zpop_withZ).
    des_x <- vrpoprob_build_designs(dff4[, xvars], Zpop_withZ[, xvars], xvars)
    des_z <- vrpoprob_build_designs(dff4[, zvars], Zpop_withZ[, zvars], zvars)

    res_array_web_withZ[[as.character(miss_p)]][[iOutcome]] <- vrpoprob_estim(
      ydata = outcome_here,
      rdata = dff4$int_rating,
      xdata = des_x$ind,
      zdata = des_z$ind,
      Nmiss = Nmiss_here,
      Wpop = Wpop_withZ,
      Xpop  = des_x$pop,
      Zpop  = des_z$pop
    )
  }
}

save(res_array_web_noZ, res_array_web_withZ,
     dff3_web, Xpop_web, Wpop_web, Zpop_withZ, Wpop_withZ,
     file = "results/results_web.RData")

toc()
