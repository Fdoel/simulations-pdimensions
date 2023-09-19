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
likert_prob <-
  function(bins,
           likert_mean = (1 + bins) / 2,
           likert_sd = bins / 4,
           method = "norm") {
    intervals = c(seq(1.5, bins - 0.5), Inf)
    x <- switch (
      method,
      "norm" =  pnorm(intervals, mean = likert_mean, sd = likert_sd),
      "unif" = punif(seq(1, bins), 0, bins),
      "cauchy" = pcauchy(intervals, likert_mean, likert_sd)
    )
    y <- dplyr::lag(x, default = 0)
    x - y
  }

mamd <- function(x) {
  mean_x <- mean(x)
  for(i in 1:length(x)) {
    x[i] <- abs(x[i] - mean_x)
  }
  median(x)/0.675
}

test1 <- likert_prob(5, method = "cauchy")

