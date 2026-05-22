var1 <- unique(dff3$married)
var2 <- unique(dff3$black)
var3 <- unique(dff3$gender_spouse)
var4 <- unique(dff3$education)

weightAll    <- sum(dff3$weights)
combinations <- list()
Wpop        <- numeric(length(var1) * length(var2) * length(var3) * length(var4))
index <- 1
for (v1 in var1) for (v2 in var2) for (v3 in var3) for (v4 in var4) {
  combinations[[index]] <- c(v1, v2, v3, v4)
  dff3W <- dff3 %>% filter(married == v1, black == v2,
                            gender_spouse == v3, education == v4)
  Wpop[index] <- sum(dff3W$weights) / weightAll
  index <- index + 1
}
Xpop        <- do.call(rbind, combinations) |> as.data.frame()
colnames(Xpop) <- c("married", "black", "gender_spouse", "education")
Xpop        <- Xpop[Wpop > 0, ]
Wpop       <- Wpop[Wpop > 0]

Zpop  <- Xpop
