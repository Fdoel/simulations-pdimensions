# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------


# load packages, but functions should still be called with package::function()
library("parallel")
library("simstudy")
library("ccaPP")


# control parameters for data generation
n <- 10000                                # number of observations
p <- 2                                  # number of items
prob <- c(0.05, 0.075, 0.1, 0.125, 0.3, 0.125, 0.1, 0.075, 0.05)  # probabilities of response categories
L <- length(prob)                       # number of response categories
rho_seq <- seq(0.2, 0.8, 0.1)           # target correlation between items
R <- 100                                 # number of simulation runs
seed <- 20230111                        # seed of the random number generator

# define matrix with probabilities of response categories per item
prob_mat <- matrix(prob, nrow = p, ncol = L, byrow = TRUE)

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

#Loop over different target correlations

rho_results_list <- parallel::mclapply(rho_seq, function(rho) {

  # define correlation matrix
  Rho <- matrix(rho, nrow = p, ncol = p)
  diag(Rho) <- 1
  results_list <- parallel::mclapply(seq_len(R), function(r) {

    # print simulation run
    cat(paste(Sys.time(), sprintf(": rho=%g\n", rho), sprintf(":   run = %d\n", r)))

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
      data.frame(Run = r, Target = rho, Method = "Pearson",
                 Correlation = pearson)
    }, error = function(e) NULL, warning = function(w) NULL)

    # compute Kendall correlation
    df_kendall <- tryCatch({
      kendall <- ccaPP::corKendall(data[, 1], data[, 2])
      data.frame(Run = r, Target = rho, Method = "Kendall",
                 Correlation = kendall)
    }, error = function(e) NULL, warning = function(w) NULL)

    # combine results
    rbind(df_pearson, df_kendall)

}, mc.cores = n_cores)
  # combine results into data frame
  results <- do.call(rbind, results_list)

}, mc.cores = n_cores)
  
# combine results into data frame
rho_results <- do.call(rbind, rho_results_list)

# aggregate results over the simulation runs
library("dplyr")
aggregated <- rho_results %>%
  group_by(Target, Method) %>%
  summarize(Correlation = mean(Correlation),
            .groups = "drop")
# save results to file
file_results <- "C:/Users/flori/simulations-practice/pearson_vs_kendall/results/9cat/results_n=%d.RData"
save(rho_results, n, p, prob,rho_seq ,seed, file = sprintf(file_results, n))

# print message that simulation is done
cat(paste(Sys.time(), ": finished.\n"))


# plot average results over the simulation runs
library("ggplot2")
p <- ggplot() +
  geom_line(aes(x = Target, y = Correlation, color = Method),
            data = aggregated)

# save plot to file
file_plot <- "C:/Users/flori/simulations-practice/pearson_vs_kendall/figures/9cat/results_n=%d.pdf"
pdf(file = sprintf(file_plot, n), width = 5, height = 3.5)
print(p)
dev.off()
