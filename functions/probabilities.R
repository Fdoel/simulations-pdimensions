likert_sd <- function(x) {
  for(i in 1:length(x)) {
    x[i] <- x[i]*i
  }
  mu <- mean(x)
  for(i in 1:length(x)) {
    x[i] <- x[i]*i
  }
  mean(x) - mu
}

likert_sd(c(0.05, 0.25, 0.4, 0.25, 0.05))
