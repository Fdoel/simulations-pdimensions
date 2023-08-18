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

# Function to generate base proabilities with multiple categories
probabilty_likert <- function(bins, likert_mean = (1 + bins)/2, likert_sd = bins/4) {
  intervals = c(seq(1.5, bins - 0.5), Inf)
  x <- pnorm(intervals, mean = likert_mean, sd = likert_sd)
  y <- dplyr::lag(x, default = 0)
  x - y
}

# control parameters for data generation
n <- 50                                 # number of observations
p <- 2                                  # number of items
rho <- 0.8                              # target correlation between items
R <- 100                                # number of simulation runs
seed <- 20230111                        # seed of the random number generator
l_increment= 2                          # increment for number of categories
l_max = 30                              # max amount of categories
L <- seq(3, l_max, l_increment)

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
    
    # compute Pearson correlation
    df_pearson <- tryCatch({
      pearson <- ccaPP::corPearson(data[, 1], data[, 2])
      data.frame(Run = r, Categories = l, Method = "Pearson",
                 Correlation = pearson)
    }, error = function(e) NULL, warning = function(w) NULL)
    
    # compute Kendall correlation
    df_kendall <- tryCatch({
      kendall <- ccaPP::corKendall(data[, 1], data[, 2])
      data.frame(Run = r, Categories = l, Method = "Kendall",
                 Correlation = kendall)
    }, error = function(e) NULL, warning = function(w) NULL)
    
    df_transformed_kendall <- tryCatch({
      kendall <- ccaPP::corKendall(data[, 1], data[, 2])
      transformed_kendall <- sin(kendall*pi/2)
      data.frame(Run = r, Categories = l, Method = "Transformed Kendall",
                 Correlation = transformed_kendall)
      
    }, error = function(e) NULL, warning = function(w) NULL)
    
    # combine results
    rbind(df_pearson, df_kendall, df_transformed_kendall)
    
  }, mc.cores = n_cores)
  
  # combine results into data frame
  results <- do.call(rbind, results_list)
  
}, mc.cores = n_cores)

# combine results into data frame
results_L <- do.call(rbind, results_list_L)

# aggregate results over the simulation runs
library("dplyr")
aggregated <- results_L %>%
  group_by(Categories, Method) %>%
  summarize(Correlation = mean(Correlation),
            .groups = "drop")



# save results to file
file_path_results = "Pearson_vs_kendall/results/varied_L/results"
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
p <- ggplot() +
  geom_line(aes(x = Categories, y = Correlation, color = Method),
            data = aggregated)

# save plot to file
file_path_plot = "pearson_vs_kendall/figures/varied_L/results"
pdf(file = paste(file_path_plot, file_info, ".pdf", sep = ""), width = 5, height = 3.5)
print(p)
dev.off()
