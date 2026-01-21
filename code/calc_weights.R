#calculate weights


#married x black x gender_spouse x education 
# 2 x 2 x 3 x 5


var1 <- unique(dff3$married)
var2 <- unique(dff3$black)
var3 <- unique(dff3$gender_spouse)
var4 <- unique(dff3$education)
combinations <- list()
WXpop <- numeric(60)
weightAll <- sum(dff3$weight)

index <- 1
for (v1 in var1) {
  for (v2 in var2) {
    for (v3 in var3) {
      for (v4 in var4) {
        combinations[[index]] <- c(v1, v2, v3, v4)
        dff3W <- dff3 %>% filter(married == v1,
                                 black == v2,
                                 gender_spouse == v3,
                                 education == v4)
        #recalculate the weights
        WXpop[index] <- sum(dff3W$weight)/weightAll
        index <- index + 1
      }
    }
  }
}
Xpop <- do.call(rbind, combinations) |> as.data.frame()
colnames(Xpop) <- c("married", "black", "gender_spouse", "education")

Xpop <- Xpop[WXpop>0,]
WXpop <- WXpop[WXpop>0]

Zpop <- Xpop
WZpop <- WXpop