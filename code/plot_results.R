
load("results/results.RData")

res20 <- res_array[["0.2"]]
res50 <- res_array[["0.5"]]
res70 <- res_array[["0.7"]]

outcome_levels <- list()
titles <- list()

outcomes <- c("life","economy","unemployment","media",
              "votes_accurate","religion","abortions","death")

outcome_levels[[1]] <- c("Extremely satisfied","Very satisfied","Moderately satisfied","Slightly satisfied","Not satisfied at all")
outcome_levels[[2]] <- c("Gotten much better","Gotten somewhat better","Stayed about the same","Gotten somewhat worse","Gotten much worse")
outcome_levels[[3]] <- c("Much better","Somewhat better","About the same","Somewhat worse","Much worse")
outcome_levels[[4]] <- c("None","A little","A moderate amount","A lot","A great deal")
outcome_levels[[5]] <- c("Not at all accurately","A little accurately","Moderately accurately","Very accurately","Completely accurately")
outcome_levels[[6]] <- c("Extremely important","Very important","Moderately important","A little important","Not important at all")
outcome_levels[[7]] <- c("Not at all important","Not too important","Somewhat important","Very important","Extremely important")
outcome_levels[[8]] <- c("Favor strongly","Favor not strongly","Oppose not strongly","Oppose strongly")

titles[[1]] <- "How satisfied are you with life?"
titles[[2]] <- "National economy has gotten better or worse?"
titles[[3]] <- "Unemployment is better or worse than last year?"
titles[[4]] <- "How much trust and confidence do you have in news?"
titles[[5]] <- "How accurately do you think the votes will be counted?"
titles[[6]] <- "Is religion an important part of your life?"
titles[[7]] <- "Importance of abortion issue"
titles[[8]] <- "Favor or oppose death penalty"

pdf_titles_all <- list()
pdf_titles_all[[1]] <- "life_all.pdf"
pdf_titles_all[[2]] <- "economy_all.pdf"
pdf_titles_all[[3]] <- "unemployment_all.pdf"
pdf_titles_all[[4]] <- "media_all.pdf"
pdf_titles_all[[5]] <- "votes_all.pdf"
pdf_titles_all[[6]] <- "religion_all.pdf"
pdf_titles_all[[7]] <- "abortions_all.pdf"
pdf_titles_all[[8]] <- "death_all.pdf"

levels_r <- c("Liked a great deal", "Liked a moderate amount", "Liked a little",
              "Neither liked nor disliked", "Disliked a little", "Disliked a moderate amount",
              "Disliked a great deal")

levels_r2 <- c("Never serious",
               "Some of the time serious",
               "About half of the time serious",
               "Most of the time serious",
               "Always serious")


df_hist <- dff3 %>%
  group_by(int_rating) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(prob = count / sum(count)) 

df_hist$int_rating <- factor(df_hist$int_rating)

rr <- ggplot(df_hist, aes(x = int_rating, y = prob, fill = int_rating)) + 
  geom_bar(stat = "identity", color = "gray", alpha = 0.7) +  # Bar plot with precomputed y-values
  labs(
    title = "Distribution of Rating of interview",
    x = "Rating",
    y = ""
  ) + 
  theme_minimal(base_size = 14) + 
  theme(
    axis.text.x = element_blank()
  )+
  scale_y_continuous(labels = scales::percent_format()) +  # Format as percentages
  scale_fill_manual(values = brewer.pal(n = 7, name = "Blues"),  # Lighter to darker blue-green palette
                    labels = levels_r)+
  guides(fill = guide_legend(title = "Rating"))

ggsave(paste0("plots/","rating.pdf"), rr, device = cairo_pdf,
       width = 10, height = 7)



df_hist <- dff3 %>%
  group_by(seriously) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(prob = count / sum(count)) 


df_hist$seriously <- factor(df_hist$seriously)

rr2 <- ggplot(df_hist, aes(x = seriously, y = prob, fill = seriously)) + 
  geom_bar(stat = "identity", color = "gray", alpha = 0.7) +  # Bar plot with precomputed y-values
  labs(
    title = "Distribution of How often you take survey seriously",
    x = "Serious",
    y = ""
  ) + 
  theme_minimal(base_size = 14) + 
  theme(
    axis.text.x = element_blank()
  )+
  scale_y_continuous(labels = scales::percent_format()) +  # Format as percentages
  scale_fill_manual(values = brewer.pal(n = 5, name = "Blues"),  # Lighter to darker blue-green palette
                    labels = levels_r2)+
  guides(fill = guide_legend(title = "How serious"))

ggsave(paste0("plots/","serious.pdf"), rr2, device = cairo_pdf,
       width = 10, height = 7)



for (iOutcome in c(1:8)) {
  
  dff4 <- dff3 %>% filter(eval(parse(text=paste0("dff3$", outcomes[iOutcome]))) > 0)
  outcome_here <- eval(parse(text=paste0("dff4$", outcomes[iOutcome])))
  NN <- dim(dff4)[1]
  
  # Adjusted weights
  dff4$weights_adj <- dff4$weights / sum(dff4$weights)
  
  # Calculate raw and adjusted proportions
  raw <- as.numeric(table(outcome_here) / NN)
  adj <- numeric(length(raw))
  for (iLevel in 1:length(adj)) {
    value <- sort(unique(outcome_here))[iLevel]
    adj[iLevel] <- sum(dff4$weights_adj[outcome_here == value]) / sum(dff4$weights_adj)
  }
  
  # Extract estimated proportions
  est20 <- res20[[iOutcome]]$pphat
  est50 <- res50[[iOutcome]]$pphat
  est70 <- res70[[iOutcome]]$pphat
  
  levels_x <- outcome_levels[[iOutcome]]
  
  # Create a data frame
  df <- data.frame(
    Category = factor(levels_x, levels = levels_x, ordered = TRUE),  # Ensure correct order
    Raw = raw,
    Adjusted = adj,
    Est20 = est20,
    Est50 = est50,
    Est70 = est70
  ) %>%
    pivot_longer(cols = c(Raw, Adjusted, Est20, Est50, Est70), names_to = "Method", values_to = "Probability")
  
  df$Method <- factor(df$Method, levels = c("Raw", "Adjusted", "Est20", "Est50", "Est70"))
  
  color_palette <- c(
    "Raw" = "#56B4E9",      # Light Blue
    "Adjusted" = "#999999", # Gray
    "Est20" = "#E69F00",    # Orange
    "Est50" = "#009E73",    # Green
    "Est70" = "#D55E00"     # Red
  )
  
  pd <- position_dodge(width = 0.6)
  
  df <- df %>%
    mutate(error_lower = case_when(
      Method == "Raw"      ~ Probability - sqrt((Probability * (1 - Probability)) / NN),
      Method == "Adjusted" ~ Probability - sqrt((Probability * (1 - Probability)) / NN),
      Method == "Est20"    ~ Probability - res20[[iOutcome]]$pphat_se,
      Method == "Est50"    ~ Probability - res50[[iOutcome]]$pphat_se,
      Method == "Est70"    ~ Probability - res70[[iOutcome]]$pphat_se,
      TRUE ~ NA_real_
    ),
    error_upper = case_when(
      Method == "Raw"      ~ Probability + sqrt((Probability * (1 - Probability)) / NN),
      Method == "Adjusted" ~ Probability + sqrt((Probability * (1 - Probability)) / NN),
      Method == "Est20"    ~ Probability + res20[[iOutcome]]$pphat_se,
      Method == "Est50"    ~ Probability + res50[[iOutcome]]$pphat_se,
      Method == "Est70"    ~ Probability + res70[[iOutcome]]$pphat_se,
      TRUE ~ NA_real_
    ))
  
  # Plot with properly filtered error bars
  p <- ggplot(df, aes(x = Category, y = Probability, fill = Method)) +
    geom_bar(stat = "identity", position = pd, width = 0.6) +  # Adjust bar width
    scale_fill_manual(values = color_palette,
                      labels = c("Raw" = "Unadjusted", "Adjusted" = "Adjusted",
                                 "Est20" = "Non-response 20%", "Est50" = "Non-response 50%",
                                 "Est70" = "Non-response 70%")) +  # Custom colors
    labs(title = paste0("<span style='font-size:25pt;'>", titles[[iOutcome]], "</span>",
                        "<br>(<span style='color:#E69F00;'>ρ = ",  
                        round(res20[[iOutcome]]$rho, 3), "</span>, ",  
                        "<span style='color:#009E73;'>ρ = ",  
                        round(res50[[iOutcome]]$rho, 3), "</span>, ",  
                        "<span style='color:#D55E00;'>ρ = ",  
                        round(res70[[iOutcome]]$rho, 3), "</span>)"),
         x = "",
         y = "Probability",
         fill = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 30, hjust = 1),  # Rotate x labels
      panel.grid.major.x = element_blank(),
      plot.title = element_markdown(),
      text = element_text(family = "Arial Unicode MS")  
    ) +
    scale_y_continuous(labels = label_percent(scale = 100)) +
    geom_errorbar(aes(ymin = error_lower, ymax = error_upper), 
                  position = pd, width = 0.25, color = "black", na.rm = TRUE)  
  
  ggsave(paste0("plots/",pdf_titles_all[[iOutcome]]), p, 
         device = cairo_pdf,
         width = 10, height = 7, units = "in")
  
  

  dff4_filtered <- dff3 %>% filter(eval(parse(text=paste0("dff3$", outcomes[iOutcome]))) > 0)
  outcome_here <- eval(parse(text=paste0("dff4$", outcomes[iOutcome])))
  
  
  # Ensure variables are factors
  dff4_filtered$outcome_var <- as.factor(outcome_here)
  dff4_filtered$int_rating <- as.factor(dff4_filtered$int_rating)
  
  # Compute proportions for each int_rating category
  df_probs <- dff4_filtered %>%
    group_by(int_rating, outcome_var) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(int_rating) %>%
    mutate(prob = count / sum(count))  # Normalize within each int_rating
  
  # Determine the number of levels in int_rating
  n_levels <- length(levels(df_probs$int_rating))
  
  # Plot with custom legend labels using a lighter-to-darker palette and rotated x-axis labels
  q <- ggplot(df_probs, aes(x = outcome_var, y = prob, fill = int_rating)) + 
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Probability bars
    scale_y_continuous(labels = scales::percent_format()) +  # Format as percentages
    scale_fill_manual(values = brewer.pal(n = n_levels, name = "Blues"),  # Lighter to darker blue-green palette
                      labels = levels_r) +  # Set custom legend labels from levels_r
    labs(
      title = titles[[iOutcome]],
      x = "",
      y = "",
      fill = "Rating of the interview\n(response variable)"
    ) + 
    scale_x_discrete(
      breaks = 1:length(levels_x),  # Ensure the ticks are placed at 1, 2, 3, 4, 5
      labels = levels_x  # Modify these labels as needed
    ) + 
    theme_minimal(base_size = 14) +  # Apply a minimal theme
    theme(
      plot.title  = element_text(size = 25),
      axis.text.x = element_text(angle = 30, hjust = 1)  # Rotate x-axis labels by 30 degrees
    )
  
  
  ggsave(paste0("plots/","response_",pdf_titles_all[[iOutcome]]), q, 
         device = cairo_pdf,
         width = 10, height = 7, units = "in")
  
  
  
  
  #comparison respondents/non-respondents (for 50% missing proportion)
  res <- res50[[iOutcome]]
  
  p_over <- res$pphat
  p_resp <- res$pphat_resp
  p_nonresp <- res$pphat_nonresp
  
  df_plot <- data.frame(
    category = factor(seq_along(p_resp)),
    Unconditional = p_over,
    Respondents = p_resp,
    Nonrespondents = p_nonresp
  )
  
  df_long <- pivot_longer(
    df_plot,
    cols = -category,
    names_to = "Group",
    values_to = "Probability"
  )
  
  df_all <- df_long |> mutate(Panel = "Overall")
  df_uncond <- df_long |> filter(Group == "Unconditional") |> mutate(Panel = "Unconditional")
  df_nonresp <- df_long |> filter(Group == "Nonrespondents") |> mutate(Panel = "Nonrespondents")
  df_resp <- df_long |> filter(Group == "Respondents") |> mutate(Panel = "Respondents")
  
  df_faceted <- bind_rows(df_all, df_uncond, df_nonresp, df_resp)
  
  df_faceted$Panel <- factor(
    df_faceted$Panel,
    levels = c("Overall", "Unconditional", "Nonrespondents", "Respondents")
  )
  
  df_se <- data.frame(
    category = factor(seq_along(res$pphat_se)),
    Unconditional = res$pphat_se,
    Respondents = res$pphat_resp_se,
    Nonrespondents = res$pphat_nonresp_se
  )
  
  df_se_long <- pivot_longer(
    df_se,
    cols = -category,
    names_to = "Group",
    values_to = "SE"
  )
  
  df_faceted <- df_faceted |>
    left_join(df_se_long, by = c("category", "Group"))
  
  p <- ggplot(df_faceted, aes(x = category, y = Probability, fill = Group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(
        ymin = pmax(Probability - SE, 0),
        ymax = Probability + SE
      ),
      position = position_dodge(width = 0.8),
      width = 0.25
    ) +
    scale_fill_manual(
      values = c(
        "Unconditional" = "#B3B3B3",
        "Respondents" = "steelblue",
        "Nonrespondents" = "firebrick"
      )
    ) +
    scale_x_discrete(labels = outcome_levels[[iOutcome]]) +
    scale_y_continuous(
      limits = c(0, max(df_faceted$Probability + df_faceted$SE, na.rm = TRUE) * 1.15),
      expand = expansion(mult = c(0, 0)),
      labels = label_percent(scale = 100)
    ) +
    facet_wrap(~Panel, nrow = 2, ncol = 2) +
    labs(
      title = paste0(
        "<span style='font-size:25pt;'>", titles[[iOutcome]], "</span>",
        "<br><span style='color:#000000;'>(ρ = ",
        round(res$rho, 3), ", 50% missing)","</span>"
      ),
      x = "Outcome category",
      y = "Population proportion",
      fill = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.major.x = element_blank(),
      plot.title = element_markdown(),
      text = element_text(family = "Arial Unicode MS")
    )
  
  ggsave(paste0("plots/","response_comp_",pdf_titles_all[[iOutcome]]), p, 
         device = cairo_pdf,
         width = 10, height = 7, units = "in")
  
  
}
