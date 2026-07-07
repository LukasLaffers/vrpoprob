install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

packages <- c(
  "tictoc", "readr", "ggplot2", "ggtext", "scales",
  "dplyr", "tidyr", "RColorBrewer", "latex2exp",
  "mnorm", "parallel", "maxLik", "numDeriv", "MASS"
)

invisible(lapply(packages, install_if_missing))


df <- read_csv("data/anes_timeseries_2024_csv_20250219.csv")


# WEIGHTS
# V240105a - all weights - face-to-face, web, paper

# COVARIATES
# V241461x - marital status 2 categories - 1. married
# V241520  - gender of the spouse
# V241501x - race - 1. white, 2. black
# V241465x - education - 5 levels

# OUTCOME VARIABLES
# V241621  - how satisfied with life 1. extremely ... 5. not at all
# V241294x - economy got better 1. much better ... 5. much worse
# V241300x - unemployment 1. much better ... 5. much worse
# V241303  - importance of abortions 1. not at all ... 5. extremely important
# V241308x - death penalty 1. favor strongly ... 4. oppose strongly
# V241314  - votes counted accurately 1. not at all ... 5. completely accurately
# V241335  - how much trust in media 1. none ... 5. a great deal
# V241420  - religion is important 1. extremely ... 5. not at all

# RESPONSE VARIABLES
# V241618  - rate interview: 1. liked a great deal ... 7. disliked a great deal
# V241619  - how often you took survey seriously: 1. never ... 5. always

# WEB-ONLY
# V241620  - internet access ease (IW_ONLINE): 1-5


dff <- df %>% dplyr::select(
  "V240002", "V240105a",
  "V241461x", "V241520", "V241501x", "V241465x",
  "V241621", "V241294x", "V241300x", "V241335",
  "V241314", "V241420", "V241303", "V241308x",
  "V241618", "V241619",
  "V241620"
)

colnames(dff) <- c(
  "mode", "weights",
  "marital", "gender_spouse", "race", "education",
  "life", "economy", "unemployment", "media",
  "votes_accurate", "religion", "abortions", "death",
  "int_rating", "seriously",
  "online"
)

dff <- as.data.frame(dff)

dff3 <- dff %>% filter(
  marital       >  0,
  gender_spouse > -2,
  gender_spouse <  3,
  race          >  0,
  education     >  0,
  int_rating    >  0,
  seriously     >  0
)

dff3$married <- 1 * (dff3$marital == 1)
dff3$black   <- 1 * (dff3$race   == 2)
