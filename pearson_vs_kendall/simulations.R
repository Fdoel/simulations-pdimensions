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
library("ltm")

# control parameters for data generation
n <- 1000                               # number of observations
p_low <- 2                              # minimal number of items
p_high <- 2                             # maximal number of items
prob <- c(0.05, 0.25, 0.4, 0.25, 0.05)  # probabilities of response categories
L <- length(prob)                       # number of response categories
rho <- 0.7                              # target correlation between items
R <- 10                                 # number of simulation runs
seed <- 20230111                        # seed of the random number generator

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

p_results_list <- parallel::mclapply(seq(p_low, p_high), function(p) {
  # define matrix with probabilities of response categories per item
  prob_mat <- matrix(prob, nrow = p, ncol = L, byrow = TRUE)

  # define correlation matrix
  Rho <- matrix(rho, nrow = p, ncol = p)
  diag(Rho) <- 1
  
  results_list <- parallel::mclapply(seq_len(R), function(r) {
    # print simulation run
    cat(paste(Sys.time(), sprintf(":   run = %d\n", r)))

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
      pearson_matrix = cor(data, method = c("pearson"))
      data.frame(Items = p, Run = r, Method = "Pearson",
              Correlation = pearson_matrix)
    }, error = function(e) NULL, warning = function(w) NULL)
    # compute Kendall correlation
    df_kendall <- tryCatch({
      kendall_matrix <- cor(data, method = c("kendall"))
      data.frame(Items = p, Run = r, Method = "Kendall",
              Correlation = kendall_matrix)
    }, error = function(e) NULL, warning = function(w) NULL)

    # combine results
    rbind(df_pearson, df_kendall) 
    

}, mc.cores = n_cores)
# combine results into data frame
results <- do.call(rbind, results_list)

# aggregate results over the simulation runs

library("dplyr")

#First calculate the average correlation matrix under Pearson
pearson = data.matrix(filter(results, Method == "Pearson")[-1:-3])
print(pearson)
X <- matrix(0, ncol = p, nrow = p)
for(i in seq(1, (R*p), p)) {
  X <- X + pearson[i:(i+p-1), 1:p]
}

aggregated_pearson = data.frame(Items = p, Method = "Pearson", Covariance = X/R)

# Do the same for Kendall correlation

kendall = data.matrix(filter(results, Method == "Kendall")[-1:-3])
Y <<- matrix(0, ncol = p, nrow = p)

# p*p matrices therefore


for(j in seq(1, R*p, p)) {
  Y <- Y + kendall[j:(j+p-1), 1:p]
}

aggregated_kendall = data.frame(Items = p, Method = "Kendall", Covariance = Y/R)

#Combine the results into a data frame
rbind(aggregated_pearson, aggregated_kendall)


}, mc.cores = n_cores)


# print message that simulation is done
#cat(paste(Sys.time(), ": finished.\n"))


# plot average results over the simulation runs
#library("ggplot2")
#p <- ggplot() +
#  geom_line(aes(x = p, y = Correlation, color = Method),
#            data = aggregated)

# save plot to file
#file_plot <- "C:/Users/flori/simulations-practice/pearson_vs_kendall/results/results_n=%d.pdf"
#pdf(file = sprintf(file_plot, n), width = 5, height = 3.5)
#print(p)



# save results to file
file_results <- "C:/Users/flori/simulations-practice/pearson_vs_kendall/results/results_n=%d.RData"
save(results, n, p, prob, rho, seed, file = sprintf(file_results, n))

dev.off()
