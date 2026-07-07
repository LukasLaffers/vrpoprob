tic("Full computation")

# Covariates entering the latent indices. Categorical covariates are expanded
# to dummies inside vrpoprob_build_designs(); see code/vrpoprob.R.
covars <- c("married", "black", "gender_spouse", "education")

res_array <- list()

for (miss_p in miss_props) {

  res_array[[as.character(miss_p)]] <- list()
  message("Running missing proportion: ", round(100 * miss_p), "%")

  for (iOutcome in seq_along(outcomes)) {

    message("  Calculating outcome: ", outcomes[iOutcome])

    dff4         <- dff3[dff3[[outcomes[iOutcome]]] > 0, ]
    outcome_here <- dff4[[outcomes[iOutcome]]]
    NN           <- nrow(dff4)
    Nmiss_here   <- round(NN * miss_p / (1 - miss_p))

    message("    Sample size NN = ", NN, "; Nmiss used = ", Nmiss_here)

    des <- vrpoprob_build_designs(dff4[, covars], Xpop[, covars], covars)

    res <- vrpoprob_estim(
      ydata = outcome_here,
      rdata = dff4$int_rating,
      xdata = des$ind,
      zdata = des$ind,
      Nmiss = Nmiss_here,
      Wpop = Wpop,
      Xpop  = des$pop,
      Zpop  = des$pop
    )

    res_array[[as.character(miss_p)]][[iOutcome]] <- res
  }
}

save(res_array, file = "results/results.RData")
message("All results saved to results/results.RData")


# Summary table: unadjusted vs corrected shares + rho 95% CI at p_miss = 0.65

labels <- c("Life satisfaction", "National economy", "Unemployment",
            "Media trust", "Votes counted accurately",
            "Religion important", "Abortion important", "Death penalty")

res65 <- res_array[["0.65"]]
rows  <- lapply(seq_along(outcomes), function(i) {
  r    <- res65[[i]]
  dff4 <- dff3[dff3[[outcomes[i]]] > 0, ]
  y    <- dff4[[outcomes[i]]]
  w    <- dff4$weights / sum(dff4$weights)
  M    <- max(y)
  raw_w <- sapply(1:M, function(j) sum(w[y == j]))
  data.frame(
    outcome   = labels[i],
    raw       = paste(round(100 * raw_w),   collapse = " / "),
    corrected = paste(round(100 * r$pphat),  collapse = " / "),
    rho       = round(r$rho,    3),
    rho_se    = round(r$rho_se, 3),
    ci_lo     = round(r$rho - 1.96 * r$rho_se, 3),
    ci_hi     = round(r$rho + 1.96 * r$rho_se, 3)
  )
})
tab <- do.call(rbind, rows)
write.csv(tab, "results/summary_table.csv", row.names = FALSE)
message("Summary table written to results/summary_table.csv")

toc()
