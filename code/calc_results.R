tic("Full computation")

miss_props <- c(0.2, 0.5, 0.7)

outcome_levels <- list()
titles <- list()
res_array <- list()

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
titles[[7]] <- "Importance of abortion issue."
titles[[8]] <- "Favor or oppose death penalty"

for (miss_p in miss_props) {
  
  res_array[[as.character(miss_p)]] <- list()
  message("Running missing proportion: ", round(100*miss_p), "%")
  
  for (iOutcome in seq_along(outcomes)) {
    
    message("  Calculating outcome: ", outcomes[iOutcome])
    
    dff4 <- dff3 %>% filter(eval(parse(text = paste0("dff3$", outcomes[iOutcome]))) > 0)
    outcome_here <- eval(parse(text = paste0("dff4$", outcomes[iOutcome])))
    NN <- nrow(dff4)
    Nmiss_here <- round(NN * miss_p / (1 - miss_p))
    
    message("    Sample size NN = ", NN, "; Nmiss used = ", Nmiss_here)
    
    res <- vrpoprob_estim(
      ydata = outcome_here,
      rdata = dff4$int_rating,
      xdata = data.matrix(dff4 %>% select(married, black, gender_spouse, education)),
      zdata = data.matrix(dff4 %>% select(married, black, gender_spouse, education)),
      Nmiss = Nmiss_here,
      WXpop = WXpop,
      Xpop = data.matrix(Xpop),
      WZpop = WZpop,
      Zpop = data.matrix(Zpop)
    )
    
    res_array[[as.character(miss_p)]][[iOutcome]] <- res
  }
}

save(res_array, file = "results/results.RData")
message("All results saved to results/results.RData")


toc()
