# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
#         Floris van den Doel
#         Erasmus University Rotterdam
# ------------------------------------


# load packages, but functions should still be called with package::function()
library("parallel")
library("simstudy")
library("ccaPP")
source("functions/block_correlation_matrix_fixed.R")
source("functions/probabilities.R")
source("functions/Pn.R")


# control parameters for data generation
n <- 100                                # number of observations
num_scales <- 4                         # number of scales
p <- 5                                  # number of items per scale
prob = c(0.05, 0.25, 0.4, 0.25, 0.05)   # get base probabilities
L <- length(prob)                       # number of response categories
rho_b <- 0.4                            # target correlation between scales
rho_w_seq <- seq(0.4, 0.8, 0.1)         # target correlation within scales
R <- 50                                 # number of simulation runs
seed <- 20230111                        # seed of the random number generator
scale_estimator <- "Pn"                # type of scale estimator used
true_sd <- likert_sd(prob)              # get the true variance of each variable

# define matrix with probabilities of response categories per item
prob_mat <-  matrix(prob, nrow = num_scales * p, ncol = L, byrow = TRUE)

# control parameters for random respondents
epsilons <- seq(0, 0.3, by = 0.05)
epsilon_max <- max(epsilons)

# it is very easy to use parallel computing on Unix systems, but not on Windows
if (.Platform$OS.type == "windows") {
  n_cores <- 1              # use only one CPU core
} else {
  n_cores <- 2              # number of CPU cores to be used
  RNGkind("L'Ecuyer-CMRG")  # use parallel random number streams
}

# loop over different within correlation levels
results_list_rho <- parallel::mclapply(rho_w_seq, function(rho_w) {
  # define correlation matrix
  Rho <- cor_mat_block(num_scales, p, rho_w, rho_b)
  
  # define matrix with variance of each variable on the diagonal
  V <- matrix(0, ncol = p * num_scales, nrow = p * num_scales)
  diag(V) <- true_sd
  
  # Multiply it with standard deviation matrices to get Sigma
  Sigma <- V %*% Rho %*% V
  
  # run simulation
  cat(paste(Sys.time(), ": starting ...\n"))
  set.seed(seed)
  results_list <- parallel::mclapply(seq_len(R), function(r) {
    # print simulation run
    cat(paste(Sys.time(), sprintf(":   run = %d\n", r)))
    
    # initialize data table
    initial <- simstudy::genData(n)
    # generate correlated rating-scale items
    data <-
      simstudy::genOrdCat(
        initial,
        baseprobs = prob_mat,
        prefix = "item",
        corMatrix = Rho
      )
    # drop ID and convert to integer matrix
    data <- as.matrix(data[,-1])
    storage.mode(data) <- "integer"
    
    # generate probabilities of being a random respondents
    careless_probabilities <- runif(n)
    
    # order observations according to probabilities of being random respondents,
    # which makes it easier to keep previous careless respondents the same as the
    # contamination level increases (for maximum comparability)
    order <- order(careless_probabilities)
    data <- data[order,]
    careless_probabilities <- careless_probabilities[order]
    
    # generate random responses to be used for careless respondents
    n_careless_max <- sum(careless_probabilities < epsilon_max)
    data_careless <-
      replicate(p, sample.int(L, n_careless_max, replace = TRUE))
    
    # loop over contamination levels
    results_r <- lapply(epsilons, function(epsilon) {
      # turn selected observations into careless respondents: since the
      # observations are sorted according to the probability of being careless,
      # this keeps previous careless respondents the same as the contamination
      # level increases
      if (epsilon > 0) {
        careless <- which(careless_probabilities < epsilon)
        data[careless,] <- data_careless[careless,]
      }
      
      # Create a scale matrix with variance on the diagonal
      scale_estimate = matrix(0,
                              nrow = p * num_scales,
                              ncol = p * num_scales)
      for (i in 1:(p * num_scales)) {
        scale_estimate[i, i] <-
          #Change this line to preffered scale estimator
          Pn(data[, i])
      }
      
      # compute Frobenius norm of the estimation error matrix for Pearson
      df_pearson_norm <- tryCatch({
        if (det(scale_estimate) == 0) {
          NULL
        } else {
          # get the transformed Kendall matrix and multiply it with scale to get the scatter matrix
          pearson_matrix <-
            scale_estimate %*% cor(data, method = "pearson") %*% scale_estimate
          # get the Frobenius norm of the error matrix
          pearson_norm <- norm(Sigma - pearson_matrix, type = "F")
          data.frame(
            Run = r,
            rho_w = rho_w,
            epsilon = epsilon,
            Method = "Pearson",
            Norm = pearson_norm
          )
        }
      }, error = function(e)
        NULL, warning = function(w)
          NULL)
      
      # compute Frobenius norm of the estimation error matrix for Pearson
      df_trans_kendall_norm <- tryCatch({
        if (det(scale_estimate) == 0) {
          NULL
        } else {
          # get the transformed Kendall matrix and multiply it with scale to get the scatter matrix
          kendall_matrix <-
            scale_estimate %*% sin(pi * 0.5 * cor(data, method = "kendall")) %*% scale_estimate
          # get the Frobenius norm of the error matrix
          kendall_norm <- norm(Sigma - kendall_matrix, type = "F")
          data.frame(
            Run = r,
            rho_w = rho_w,
            epsilon = epsilon,
            Method = "Transformed Kendall",
            Norm = kendall_norm
          )
        }
      }, error = function(e)
        NULL, warning = function(w)
          NULL)
      
      # compute Frobenius norm of the estimation error matrix for Pearson
      df_covariance_norm <- tryCatch({
        # get the Frobenius norm of the error matrix
        cov_norm <- norm(Sigma - cov(data), type = "F")
        data.frame(
          Run = r,
          rho_w = rho_w,
          epsilon = epsilon,
          Method = "Covariance",
          Norm = cov_norm
        )
        
      }, error = function(e)
        NULL, warning = function(w)
          NULL)
      
      # combine results
      rbind(df_pearson_norm, df_trans_kendall_norm, df_covariance_norm)
      
    })
    
    # combine results from current simulation run into data frame
    do.call(rbind, results_r)
    
  }, mc.cores = n_cores)
  
  # combine results into data frame
  results <- do.call(rbind, results_list)
  
}, mc.cores = n_cores)

# combine results into dataframe
results_rho <- do.call(rbind, results_list_rho)

# Get the amount of missing observations

total_estimators_possible <- R * length(rho_w_seq) * length(epsilons)
total_estimators_actual <- nrow(results_rho) / 3
useless_count <- total_estimators_possible - total_estimators_actual

# save results to file
file_results <-
  "pearson_vs_kendall/results/Qn/Rhow/results_n=%d-scales=%d-p=%d-est=%s.RData"
save(
  results_rho,
  n,
  p,
  num_scales,
  prob,
  rho_w_seq,
  rho_b,
  useless_count,
  seed,
  file = sprintf(file_results, n, num_scales, p, scale_estimator)
)

# print message that simulation is done
cat(paste(Sys.time(), ": finished.\n"))


# aggregate results over the simulation runs
library("dplyr")
aggregated <- results_rho %>%
  group_by(epsilon, Method, rho_w) %>%
  summarize(Norm = mean(Norm),
            .groups = "drop")

# plot average results over the simulation runs
library("ggplot2")
p_line <-
  ggplot(aggregated, aes(x = rho_w, y = Norm, color = Method)) +
  geom_line() +
  facet_wrap( ~ factor(epsilon), nrow = 1) +
  
  # add some cosmetic changes
  scale_x_continuous(n.breaks = 3) +
  theme_minimal() +
  theme(panel.spacing = unit(2, "lines")) +
  xlab("Correlation between items in scale") +
  labs(
    title = sprintf(
      "Frobesius norm of variance error matrix using %s",
      scale_estimator
    ),
    subtitle = sprintf("n =% d correlation between scales = %f", n, rho_b),
    caption = sprintf(
      "number of scales = %d, number of items per scale = %d, possible runs = %d, missing runs = %f",
      num_scales,
      p,
      total_estimators_possible,
      useless_count
    )
  )

# save plot to file
file_plot <-
  "pearson_vs_kendall/figures/Qn/Rhow/results_n=%d-scales=%d-p=%d-estimator=%s.pdf"
pdf(
  file = sprintf(file_plot, n, num_scales, p, scale_estimator),
  width = 8,
  height = 5
)
print(p_line)
dev.off()
