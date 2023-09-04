# Function to get the standard deviation of a likert scale
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

# Function to get likert base probabilities
probabilty_likert <- function(bins, likert_mean = (1 + bins)/2, likert_sd = bins/4) {
  intervals = c(seq(1.5, bins - 0.5), Inf)
  x <- pnorm(intervals, mean = likert_mean, sd = likert_sd)
  y <- dplyr::lag(x, default = 0)
  x - y
}
