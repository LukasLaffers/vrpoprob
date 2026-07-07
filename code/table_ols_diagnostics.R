# Recreate Table 4: OLS instrument diagnostics with demographic controls.
#
# Run from the replication package root:
#   Rscript code/table_ols_diagnostics.R
#
# Output:
#   results/table_ols_diagnostics.csv

data_file <- "data/anes_timeseries_2024_csv_20250219.csv"
if (!file.exists(data_file)) {
  stop("Could not find ", data_file, ". Run this script from the replication package root.")
}

dir.create("results", showWarnings = FALSE, recursive = TRUE)

df <- read.csv(data_file)

dff <- df[, c(
  "V240002", "V240105a",
  "V241461x", "V241520", "V241501x", "V241465x",
  "V241621", "V241294x", "V241300x", "V241335",
  "V241314", "V241420", "V241303", "V241308x",
  "V241618", "V241619", "V241620"
)]

colnames(dff) <- c(
  "mode", "weights",
  "marital", "gender_spouse", "race", "education",
  "life", "economy", "unemployment", "media",
  "votes_accurate", "religion", "abortions", "death",
  "int_rating", "seriously", "online"
)

# Web-only Appendix D sample with valid demographic controls, response proxies,
# and IW_ONLINE. Outcome-specific missingness is handled below.
dff <- subset(
  dff,
  mode == 2 &
    marital > 0 &
    gender_spouse > -2 & gender_spouse < 3 &
    race > 0 &
    education > 0 &
    int_rating > 0 &
    seriously > 0 &
    online > 0
)

dff$married <- as.integer(dff$marital == 1)
dff$black <- as.integer(dff$race == 2)

# Table 4 uses OLS with demographic controls. Married and black are binary
# indicators; spouse/partner gender and education are categorical controls and
# therefore enter as factors/dummy indicators.
control_terms <- c(
  "married",
  "black",
  "factor(gender_spouse)",
  "factor(education)"
)

relevance <- c(
  int_rating = "int_rating (1 liked -- 7 disliked)",
  seriously = "seriously (1 never -- 5 always)"
)

exclusion <- c(
  life = "Life satisfaction",
  economy = "National economy",
  unemployment = "Unemployment",
  media = "Media trust",
  votes_accurate = "Votes counted accurately",
  religion = "Religion important",
  abortions = "Abortion important",
  death = "Death penalty"
)

run_online_ols <- function(dep_var, label, block) {
  d <- dff
  if (!dep_var %in% names(relevance)) {
    d <- d[d[[dep_var]] > 0, ]
  }

  fml <- as.formula(paste(dep_var, "~ online +", paste(control_terms, collapse = " + ")))
  fit <- lm(fml, data = d)
  s <- coef(summary(fit))["online", ]

  data.frame(
    block = block,
    dependent_variable = label,
    coef_online = unname(s["Estimate"]),
    se = unname(s["Std. Error"]),
    t = unname(s["t value"]),
    p = unname(s["Pr(>|t|)"]),
    n = nrow(d),
    controls = paste(control_terms, collapse = " + "),
    stringsAsFactors = FALSE
  )
}

tab4 <- rbind(
  do.call(rbind, Map(run_online_ols, names(relevance), relevance, "Relevance")),
  do.call(rbind, Map(run_online_ols, names(exclusion), exclusion, "Exclusion"))
)

write.csv(tab4, "results/table_ols_diagnostics.csv", row.names = FALSE)

print(tab4[, c("block", "dependent_variable", "coef_online", "se", "t", "p", "n")])
