source(file.path(getwd(), "/IHWpaper-master/R/ddhw_wrappers.R"))
source(file.path(getwd(), "/IHWpaper-master/R/benchmarking.R"))
source(file.path(getwd(), "/R/method_new_binning_admm.R"))


## NULL Simulation Benchmark
#----------------- General benchmark settings -------------------------------#
alphas <- seq(0.01, .1, length=10)
nreps <- 4000 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000

sim_funs <- list(null_sim_fun(ms))

fdr_methods <- list(ihw_binning_5fold, ihw_5fold)
fdr_methods <- lapply(fdr_methods, continuous_wrap)

results <- list()
for (i in 1:length(sim_funs)){
  res <- sim_fun_eval(sim_funs[[i]], fdr_methods, nreps = nreps, alphas, BiocParallel=F)
  results[[i]] <- res
}

final_output_sim2 <- do.call(rbind, results)

save(final_output_sim2, file = "final_results_binning_rep4000_sim_null.RData")