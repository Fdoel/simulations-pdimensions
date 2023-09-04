likert_sd <- function(x) {
  for(i in 1:length(x)) {
    x[i] <- x[i]*i
  }
  mu <- sum(x)
  for(i in 1:length(x)) {
    x[i] <- x[i]*i
  }
  sqrt(sum(x) - mu*mu)
}

likert_sd(c(0.05, 0.25, 0.4, 0.25, 0.05))
