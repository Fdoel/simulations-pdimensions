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

# Function to generate base probabilities with multiple categories
probabilty_likert <- function(bins, likert_mean = (1 + bins)/2, likert_sd = bins/4) {
  intervals = c(seq(1.5, bins - 0.5), Inf)
  x <- pnorm(intervals, mean = likert_mean, sd = likert_sd)
  y <- dplyr::lag(x, default = 0)
  x - y
}

# control parameters for data generation
n <- 100                               # number of observations
p <- 2                                  # number of items
rho <- 0.8                              # target correlation between items
R <- 100                                # number of simulation runs
seed <- 20230111                        # seed of the random number generator
l_increment= 30                          # increment for number of categories
l_max = 1000                              # max amount of categories
L <- seq(3, l_max, l_increment)

# control parameters for random respondents
epsilons <- seq(0, 0.3, by = 0.3)
epsilon_max <- max(epsilons)

# define correlation matrix
Rho <- matrix(rho, nrow = p, ncol = p)
diag(Rho) <- 1

# it is very easy to use parallel computing on Unix systems, but not on Windows
if (.Platform$OS.type == "windows") {
  n_cores <- 1              # use only one CPU core
} else {
  n_cores <- 2              # number of CPU cores to be used
  RNGkind("L'Ecuyer-CMRG")  # use parallel random number streams
}

# run simulation
cat(paste(Sys.time(), ": starting ...\n"))
set.seed(seed)

results_list_L <- parallel::mclapply(L, function(l) {
  # define matrix with probabilities of response categories per item
  prob <- probabilty_likert(l)
  prob_mat <- matrix(prob, nrow = p, ncol = l, byrow = TRUE)

  results_list <- parallel::mclapply(seq_len(R), function(r) {
    # print simulation run
    cat(paste(Sys.time(), sprintf(": L=%d\n", l), sprintf(":   run = %d\n", r)))

    # initialize data table
    initial <- simstudy::genData(n)
    # generate correlated rating-scale items
    data <- simstudy::genOrdCat(initial, baseprobs = prob_mat, prefix = "item",
                                corMatrix = Rho)
    # drop ID and convert to integer matrix
    data <- as.matrix(data[, -1])
    storage.mode(data) <- "integer"
    
    # generate probabilities of being a random respondents
    careless_probabilities <- runif(n)
    
    # order observations according to probabilities of being random respondents,
    # which makes it easier to keep previous careless respondents the same as the
    # contamination level increases (for maximum comparability)
    order <- order(careless_probabilities)
    data <- data[order, ]
    careless_probabilities <- careless_probabilities[order]
    
    # generate clustered responses to be used for careless respondents
    n_careless_max <- sum(careless_probabilities < epsilon_max)
    careless_prob_matrix <- probabilty_likert(l, l/1.1) %>%
    matrix(nrow = p, ncol = l, byrow = TRUE)
    data_careless <- simstudy::genData(n_careless_max) %>%
    simstudy::genOrdCat(baseprobs = careless_prob_matrix, prefix = "item", corMatrix = Rho)
    data_careless <- as.matrix(data_careless[, -1])
    storage.mode(data_careless) <- "integer"
    
    # In case you want to go for careless opposed to clustered data uncomment following 2 lines and comment previous block
    #n_careless_max <- sum(careless_probabilities < epsilon_max)
    #data_careless <- replicate(p, sample.int(l, n_careless_max, replace = TRUE))
    
    
    # loop over contamination levels
    results_r <- lapply(epsilons, function(epsilon) {
      
      # turn selected observations into careless respondents: since the
      # observations are sorted according to the probability of being careless,
      # this keeps previous careless respondents the same as the contamination
      # level increases
      if (epsilon > 0) {
        careless <- which(careless_probabilities < epsilon)
        data[careless, ] <- data_careless[careless, ]
      }
      
      # compute Pearson correlation
      df_pearson <- tryCatch({
        pearson <- ccaPP::corPearson(data[, 1], data[, 2])
        data.frame(Run = r, Epsilon = epsilon, Categories = l, Method = "Pearson",
                  Correlation = pearson)
      }, error = function(e) NULL, warning = function(w) NULL)
    
      # compute Kendall correlation
      df_kendall <- tryCatch({
        kendall <- ccaPP::corKendall(data[, 1], data[, 2])
        data.frame(Run = r, Epsilon = epsilon, Categories = l, Method = "Kendall",
                   Correlation = kendall)
      }, error = function(e) NULL, warning = function(w) NULL)
      
      df_transformed_kendall <- tryCatch({
        kendall <- ccaPP::corKendall(data[, 1], data[, 2])
        transformed_kendall <- sin(kendall*pi/2)
        data.frame(Run = r, Epsilon = epsilon, Categories = l, Method = "Transformed Kendall",
                   Correlation = transformed_kendall)
        
      }, error = function(e) NULL, warning = function(w) NULL)
      
      # combine results
      rbind(df_pearson, df_kendall, df_transformed_kendall)
    
    })
    
    # combine results from current simulation run into data frame
    do.call(rbind, results_r)
    
  }, mc.cores = n_cores)
  
  # combine results into data frame
  results <- do.call(rbind, results_list)
  
}, mc.cores = n_cores)

# combine results into data frame
results_L <- do.call(rbind, results_list_L)

# aggregate results over the simulation runs
library("dplyr")
aggregated <- results_L %>%
  group_by(Categories, Epsilon, Method) %>%
  summarize(Correlation = mean(Correlation),
            .groups = "drop")

# save results to file
file_path_results = "Pearson_vs_kendall/results/varied_L/clustered/results"
file_info = paste(
  sprintf("n=%d", n),
  sprintf("max=%d", l_max),
  sprintf("increment=%d", l_increment),
  sep = "_")

save(results_L, n, p, rho, l_max, l_increment, seed, file = paste(file_path_results, file_info, ".RData", sep = ""))

# print message that simulation is done
cat(paste(Sys.time(), ": finished.\n"))

# plot average results over the simulation runs

library("ggplot2")
plot_box <-  dplyr::filter(results_L, near(Epsilon, 0.0)) %>%
  ggplot(aes(x = factor(Categories), y = Correlation, color = Method)) +
    geom_boxplot()


plot_line <- ggplot(results_L, aes(x = Categories, y = Correlation, color = Method)) +
  geom_smooth(aes(linetype = factor(Epsilon)))

# save box plot to file
file_path_plot = "pearson_vs_kendall/figures/varied_L/clustered/results"
pdf(file = paste(file_path_plot, file_info, "box", ".pdf", sep = ""), width = 5, height = 3.5)
print(plot_box)
dev.off()

# save line plot to file
file_path_plot = "pearson_vs_kendall/figures/varied_L/clustered/results"
pdf(file = paste(file_path_plot, file_info, "line", ".pdf", sep = ""), width = 5, height = 3.5)
print(plot_line)
dev.off()
print(plot_line)
