outcomes <- c("life", "economy", "unemployment", "media",
              "votes_accurate", "religion", "abortions", "death")

miss_props <- c(0.5, 0.65, 0.8)

outcome_levels <- list(
  c("Extremely satisfied", "Very satisfied", "Moderately satisfied", "Slightly satisfied", "Not satisfied at all"),
  c("Gotten much better", "Gotten somewhat better", "Stayed about the same", "Gotten somewhat worse", "Gotten much worse"),
  c("Much better", "Somewhat better", "About the same", "Somewhat worse", "Much worse"),
  c("None", "A little", "A moderate amount", "A lot", "A great deal"),
  c("Not at all accurately", "A little accurately", "Moderately accurately", "Very accurately", "Completely accurately"),
  c("Extremely important", "Very important", "Moderately important", "A little important", "Not important at all"),
  c("Not at all important", "Not too important", "Somewhat important", "Very important", "Extremely important"),
  c("Favor strongly", "Favor not strongly", "Oppose not strongly", "Oppose strongly")
)

titles <- list(
  "How satisfied are you with life?",
  "National economy has gotten better or worse?",
  "Unemployment is better or worse than last year?",
  "How much trust and confidence do you have in news?",
  "How accurately do you think the votes will be counted?",
  "Is religion an important part of your life?",
  "Importance of abortion issue",
  "Favor or oppose death penalty"
)

pdf_titles_all <- list(
  "life_all.pdf", "economy_all.pdf", "unemployment_all.pdf", "media_all.pdf",
  "votes_all.pdf", "religion_all.pdf", "abortions_all.pdf", "death_all.pdf"
)
