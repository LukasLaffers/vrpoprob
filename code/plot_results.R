if (!exists("res_array")) load("results/results.RData")

dir.create("plots/", showWarnings = FALSE, recursive = TRUE)

res50 <- res_array[["0.5"]]
res65 <- res_array[["0.65"]]
res80 <- res_array[["0.8"]]

levels_r <- c("Liked a great deal", "Liked a moderate amount", "Liked a little",
              "Neither liked nor disliked", "Disliked a little",
              "Disliked a moderate amount", "Disliked a great deal")

levels_r2 <- c("Never serious", "Some of the time serious",
                "About half of the time serious", "Most of the time serious",
                "Always serious")

color_palette <- c(
  "Raw"     = "#56B4E9",
  "Adjusted"= "#999999",
  "Est50"   = "#E69F00",
  "Est65"   = "#009E73",
  "Est80"   = "#D55E00"
)


# rating and seriously distributions

df_hist <- dff3 %>%
  group_by(int_rating) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(prob = count / sum(count),
         int_rating = factor(int_rating))

rr <- ggplot(df_hist, aes(x = int_rating, y = prob, fill = int_rating)) +
  geom_bar(stat = "identity", color = "gray", alpha = 0.7) +
  labs(title = "Distribution of Rating of interview", x = "Rating", y = "") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = brewer.pal(n = 7, name = "Blues"), labels = levels_r) +
  guides(fill = guide_legend(title = "Rating"))

ggsave("plots/rating.pdf", rr, device = cairo_pdf, width = 10, height = 7)


df_hist <- dff3 %>%
  group_by(seriously) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(prob = count / sum(count),
         seriously = factor(seriously))

rr2 <- ggplot(df_hist, aes(x = seriously, y = prob, fill = seriously)) +
  geom_bar(stat = "identity", color = "gray", alpha = 0.7) +
  labs(title = "Distribution of How often you take survey seriously", x = "Serious", y = "") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = brewer.pal(n = 5, name = "Blues"), labels = levels_r2) +
  guides(fill = guide_legend(title = "How serious"))

ggsave("plots/serious.pdf", rr2, device = cairo_pdf, width = 10, height = 7)


# per-outcome plots

for (iOutcome in seq_along(outcomes)) {

  dff4         <- dff3[dff3[[outcomes[iOutcome]]] > 0, ]
  outcome_here <- dff4[[outcomes[iOutcome]]]
  NN           <- nrow(dff4)
  levels_x     <- outcome_levels[[iOutcome]]

  dff4$weights_adj <- dff4$weights / sum(dff4$weights)

  raw <- as.numeric(table(outcome_here) / NN)
  adj <- sapply(seq_along(raw), function(j) {
    value <- sort(unique(outcome_here))[j]
    sum(dff4$weights_adj[outcome_here == value])
  })

  est50 <- res50[[iOutcome]]$pphat
  est65 <- res65[[iOutcome]]$pphat
  est80 <- res80[[iOutcome]]$pphat

  df <- data.frame(
    Category = factor(levels_x, levels = levels_x, ordered = TRUE),
    Raw = raw, Adjusted = adj, Est50 = est50, Est65 = est65, Est80 = est80
  ) %>%
    pivot_longer(cols = c(Raw, Adjusted, Est50, Est65, Est80),
                 names_to = "Method", values_to = "Probability") %>%
    mutate(
      Method = factor(Method, levels = c("Raw", "Adjusted", "Est50", "Est65", "Est80")),
      error_lower = case_when(
        Method == "Raw"       ~ Probability - sqrt((Probability * (1 - Probability)) / NN),
        Method == "Adjusted"  ~ Probability - sqrt((Probability * (1 - Probability)) / NN),
        Method == "Est50"     ~ Probability - res50[[iOutcome]]$pphat_se,
        Method == "Est65"     ~ Probability - res65[[iOutcome]]$pphat_se,
        Method == "Est80"     ~ Probability - res80[[iOutcome]]$pphat_se,
        TRUE ~ NA_real_
      ),
      error_upper = case_when(
        Method == "Raw"       ~ Probability + sqrt((Probability * (1 - Probability)) / NN),
        Method == "Adjusted"  ~ Probability + sqrt((Probability * (1 - Probability)) / NN),
        Method == "Est50"     ~ Probability + res50[[iOutcome]]$pphat_se,
        Method == "Est65"     ~ Probability + res65[[iOutcome]]$pphat_se,
        Method == "Est80"     ~ Probability + res80[[iOutcome]]$pphat_se,
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
           "<br>(<span style='color:#E69F00;'>ρ = ", round(res50[[iOutcome]]$rho, 3), "</span>, ",
           "<span style='color:#009E73;'>ρ = ",      round(res65[[iOutcome]]$rho, 3), "</span>, ",
           "<span style='color:#D55E00;'>ρ = ",      round(res80[[iOutcome]]$rho, 3), "</span>)"),
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

  ggsave(paste0("plots/", pdf_titles_all[[iOutcome]]), p,
         device = cairo_pdf, width = 10, height = 7, units = "in")


  # response-by-rating plot

  dff4$outcome_var <- factor(outcome_here)
  dff4$int_rating  <- factor(dff4$int_rating)

  df_probs <- dff4 %>%
    group_by(int_rating, outcome_var) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(int_rating) %>%
    mutate(prob = count / sum(count))

  n_levels <- length(levels(df_probs$int_rating))

  q <- ggplot(df_probs, aes(x = outcome_var, y = prob, fill = int_rating)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = brewer.pal(n = n_levels, name = "Blues"), labels = levels_r) +
    scale_x_discrete(breaks = seq_along(levels_x), labels = levels_x) +
    labs(title = titles[[iOutcome]], x = "", y = "",
         fill = "Rating of the interview\n(response variable)") +
    theme_minimal(base_size = 14) +
    theme(plot.title  = element_text(size = 25),
          axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(paste0("plots/response_", pdf_titles_all[[iOutcome]]), q,
         device = cairo_pdf, width = 10, height = 7, units = "in")


  # respondents vs nonrespondents comparison (65% missing)

  res    <- res65[[iOutcome]]
  df_plot <- data.frame(
    category       = factor(seq_along(res$pphat)),
    Unconditional  = res$pphat,
    Respondents    = res$pphat_resp,
    Nonrespondents = res$pphat_nonresp
  )

  df_long <- pivot_longer(df_plot, cols = -category,
                          names_to = "Group", values_to = "Probability")

  df_faceted <- bind_rows(
    df_long |> mutate(Panel = "Overall"),
    df_long |> filter(Group == "Unconditional")  |> mutate(Panel = "Unconditional"),
    df_long |> filter(Group == "Nonrespondents") |> mutate(Panel = "Nonrespondents"),
    df_long |> filter(Group == "Respondents")    |> mutate(Panel = "Respondents")
  ) %>%
    mutate(Panel = factor(Panel,
                          levels = c("Overall", "Unconditional", "Nonrespondents", "Respondents")))

  df_se <- data.frame(
    category       = factor(seq_along(res$pphat_se)),
    Unconditional  = res$pphat_se,
    Respondents    = res$pphat_resp_se,
    Nonrespondents = res$pphat_nonresp_se
  ) %>%
    pivot_longer(cols = -category, names_to = "Group", values_to = "SE")

  df_faceted <- df_faceted |> left_join(df_se, by = c("category", "Group"))

  p <- ggplot(df_faceted, aes(x = category, y = Probability, fill = Group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = pmax(Probability - SE, 0), ymax = Probability + SE),
                  position = position_dodge(width = 0.8), width = 0.25) +
    scale_fill_manual(values = c("Unconditional" = "#B3B3B3",
                                 "Respondents"   = "steelblue",
                                 "Nonrespondents"= "firebrick")) +
    scale_x_discrete(labels = outcome_levels[[iOutcome]]) +
    scale_y_continuous(
      limits = c(0, max(df_faceted$Probability + df_faceted$SE, na.rm = TRUE) * 1.15),
      expand = expansion(mult = c(0, 0)),
      labels = label_percent(scale = 100)
    ) +
    facet_wrap(~Panel, nrow = 2, ncol = 2) +
    labs(title = paste0(
           "<span style='font-size:25pt;'>", titles[[iOutcome]], "</span>",
           "<br><span style='color:#000000;'>(ρ = ", round(res$rho, 3), ", 65% missing)</span>"),
         x = "Outcome category", y = "Population proportion", fill = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom",
      axis.text.x        = element_text(angle = 30, hjust = 1),
      panel.grid.major.x = element_blank(),
      plot.title         = element_markdown()
    )

  ggsave(paste0("plots/response_comp_", pdf_titles_all[[iOutcome]]), p,
         device = cairo_pdf, width = 10, height = 7, units = "in")
}
