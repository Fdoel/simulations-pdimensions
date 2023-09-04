likert_sd <- function(x) {
  for(i in 1:length(x)) {
    x[i] <- x[i]*i
  }
  print(x)
  mu <- sum(x)
  for(i in 1:length(x)) {
    x[i] <- x[i]*i
  }
  print(x)
  sqrt(sum(x) - mu*mu)
}
