if (!exists("res_array_web_noZ")) load("results/results_web.RData")

dir.create("plots/web_noZ",   showWarnings = FALSE, recursive = TRUE)
dir.create("plots/web_withZ", showWarnings = FALSE, recursive = TRUE)

color_palette <- c(
  "Raw"      = "#56B4E9",
  "Adjusted" = "#999999",
  "Est50"    = "#E69F00",
  "Est65"    = "#009E73",
  "Est80"    = "#D55E00"
)


# per-spec, per-outcome bar plots

for (spec in c("noZ", "withZ")) {

  res_array <- get(paste0("res_array_web_", spec))
  res50  <- res_array[["0.5"]]
  res65  <- res_array[["0.65"]]
  res80  <- res_array[["0.8"]]
  outdir <- paste0("plots/web_", spec)
  tag    <- if (spec == "noZ") "no IW_ONLINE" else "with IW_ONLINE"

  pp <- function(r) r$pphat
  se <- function(r) r$pphat_se
  rh <- function(r) r$rho

  for (iOutcome in seq_along(outcomes)) {

    dff4         <- dff3_web[dff3_web[[outcomes[iOutcome]]] > 0, ]
    outcome_here <- dff4[[outcomes[iOutcome]]]
    NN           <- nrow(dff4)
    levels_x     <- outcome_levels[[iOutcome]]

    dff4$weights_adj <- dff4$weights / sum(dff4$weights)
    raw <- as.numeric(table(outcome_here) / NN)
    adj <- sapply(seq_along(raw), function(j) {
      value <- sort(unique(outcome_here))[j]
      sum(dff4$weights_adj[outcome_here == value])
    })

    est50 <- pp(res50[[iOutcome]]); se50 <- se(res50[[iOutcome]])
    est65 <- pp(res65[[iOutcome]]); se65 <- se(res65[[iOutcome]])
    est80 <- pp(res80[[iOutcome]]); se80 <- se(res80[[iOutcome]])

    df <- data.frame(
      Category = factor(levels_x, levels = levels_x, ordered = TRUE),
      Raw = raw, Adjusted = adj, Est50 = est50, Est65 = est65, Est80 = est80
    ) %>%
      pivot_longer(cols = c(Raw, Adjusted, Est50, Est65, Est80),
                   names_to = "Method", values_to = "Probability") %>%
      mutate(
        Method = factor(Method, levels = c("Raw", "Adjusted", "Est50", "Est65", "Est80")),
        error_lower = case_when(
          Method == "Raw"      ~ Probability - sqrt((Probability * (1 - Probability)) / NN),
          Method == "Adjusted" ~ Probability - sqrt((Probability * (1 - Probability)) / NN),
          Method == "Est50"    ~ Probability - se50,
          Method == "Est65"    ~ Probability - se65,
          Method == "Est80"    ~ Probability - se80,
          TRUE ~ NA_real_
        ),
        error_upper = case_when(
          Method == "Raw"      ~ Probability + sqrt((Probability * (1 - Probability)) / NN),
          Method == "Adjusted" ~ Probability + sqrt((Probability * (1 - Probability)) / NN),
          Method == "Est50"    ~ Probability + se50,
          Method == "Est65"    ~ Probability + se65,
          Method == "Est80"    ~ Probability + se80,
          TRUE ~ NA_real_
        )
      )

    pd <- position_dodge(width = 0.6)

    p <- ggplot(df, aes(x = Category, y = Probability, fill = Method)) +
      geom_bar(stat = "identity", position = pd, width = 0.6) +
      scale_fill_manual(values = color_palette,
                        labels = c("Raw" = "Unadjusted", "Adjusted" = "Adjusted",
                                   "Est50" = "Non-response 50%", "Est65" = "Non-response 65%",
                                   "Est80" = "Non-response 80%")) +
      labs(title = paste0(
             "<span style='font-size:25pt;'>", titles[[iOutcome]], "</span>",
             " <span style='font-size:14pt;'>[web, ", tag, "]</span>",
             "<br>(<span style='color:#E69F00;'>&#961; =", round(rh(res50[[iOutcome]]), 3), "</span>, ",
             "<span style='color:#009E73;'>&#961; =",     round(rh(res65[[iOutcome]]), 3), "</span>, ",
             "<span style='color:#D55E00;'>&#961; =",     round(rh(res80[[iOutcome]]), 3), "</span>)"),
           x = "", y = "Probability", fill = NULL) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position    = "bottom",
        axis.text.x        = element_text(angle = 30, hjust = 1),
        panel.grid.major.x = element_blank(),
        plot.title         = element_markdown()
      ) +
      scale_y_continuous(labels = label_percent(scale = 100)) +
      geom_errorbar(aes(ymin = error_lower, ymax = error_upper),
                    position = pd, width = 0.25, color = "black", na.rm = TRUE)

    ggsave(paste0(outdir, "/", pdf_titles_all[[iOutcome]]), p,
           device = cairo_pdf, width = 10, height = 7, units = "in")
  }
}


# rho comparison: CSV + summary figure

rho_rows <- list()
k <- 1
for (iOutcome in seq_along(outcomes)) {
  for (miss_p in miss_props) {
    for (spec in c("noZ", "withZ")) {
      r <- get(paste0("res_array_web_", spec))[[as.character(miss_p)]][[iOutcome]]
      rho_rows[[k]] <- data.frame(
        outcome = outcomes[iOutcome],
        miss_p  = miss_p,
        spec    = spec,
        rho     = round(r$rho, 3)
      )
      k <- k + 1
    }
  }
}
rho_tab <- do.call(rbind, rho_rows)
write.csv(rho_tab, "results/rho_web_comparison.csv", row.names = FALSE)

rho_tab$spec    <- factor(rho_tab$spec, levels = c("noZ", "withZ"),
                          labels = c("no IW_ONLINE", "with IW_ONLINE"))
rho_tab$miss_lab <- factor(paste0(round(100 * rho_tab$miss_p), "% missing"),
                           levels = c("50% missing", "65% missing", "80% missing"))

p <- ggplot(rho_tab, aes(x = outcome, y = rho, color = spec, shape = spec)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_point(size = 3, position = position_dodge(width = 0.45)) +
  facet_wrap(~miss_lab, nrow = 1) +
  scale_color_manual(values = c("#D55E00", "#0072B2")) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "", y = "ρ", color = "", shape = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x        = element_text(angle = 30, hjust = 1),
        legend.position    = "bottom",
        panel.grid.major.x = element_blank())

ggsave("plots/web_rho_comparison.pdf", p, device = cairo_pdf, width = 10, height = 4.2)
