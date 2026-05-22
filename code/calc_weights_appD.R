# calc_weights_appD.R
# Population grids for the web-only (Appendix D) robustness analysis.
# Requires: dff3 (from load_data.R)
# Produces:
#   dff3_web   -- web-mode subsample with online > 0
#   Xpop_web, WXpop_web  -- X grid: demographics (married x black x gender_spouse x education)
#   Zpop_withZ, WZpop_withZ -- Z grid: demographics x online

dff3_web <- dff3 %>% filter(mode == 2, online > 0)

var1 <- unique(dff3_web$married)
var2 <- unique(dff3_web$black)
var3 <- unique(dff3_web$gender_spouse)
var4 <- unique(dff3_web$education)
var5 <- sort(unique(dff3_web$online))

weightAll <- sum(dff3_web$weights)


# X grid: demographics only -----------------------------------------------

combinations <- list()
Wpop_web    <- numeric(length(var1)*length(var2)*length(var3)*length(var4))
index <- 1
for (v1 in var1) for (v2 in var2) for (v3 in var3) for (v4 in var4) {
  combinations[[index]] <- c(v1, v2, v3, v4)
  dff3W <- dff3_web %>% filter(married == v1, black == v2,
                                gender_spouse == v3, education == v4)
  Wpop_web[index] <- sum(dff3W$weights) / weightAll
  index <- index + 1
}
Xpop_web        <- do.call(rbind, combinations) |> as.data.frame()
colnames(Xpop_web) <- c("married", "black", "gender_spouse", "education")
Xpop_web        <- Xpop_web[Wpop_web > 0, ]
Wpop_web       <- Wpop_web[Wpop_web > 0]


# Z grid: demographics x online -------------------------------------------

combinations <- list()
Wpop_withZ  <- numeric(length(var1)*length(var2)*length(var3)*length(var4)*length(var5))
index <- 1
for (v1 in var1) for (v2 in var2) for (v3 in var3) for (v4 in var4) for (v5 in var5) {
  combinations[[index]] <- c(v1, v2, v3, v4, v5)
  dff3W <- dff3_web %>% filter(married == v1, black == v2,
                                gender_spouse == v3, education == v4,
                                online == v5)
  Wpop_withZ[index] <- sum(dff3W$weights) / weightAll
  index <- index + 1
}
Zpop_withZ        <- do.call(rbind, combinations) |> as.data.frame()
colnames(Zpop_withZ) <- c("married", "black", "gender_spouse", "education", "online")
Zpop_withZ        <- Zpop_withZ[Wpop_withZ > 0, ]
Wpop_withZ       <- Wpop_withZ[Wpop_withZ > 0]
