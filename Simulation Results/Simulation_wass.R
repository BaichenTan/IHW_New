source(file.path(getwd(), "/IHWpaper-master/R/ddhw_wrappers.R"))
source(file.path(getwd(), "/IHWpaper-master/R/benchmarking.R"))
source(file.path(getwd(), "/R/method_new_binning_admm.R"))

## Wasserman Normal Simulation Benchmark
#----------------- General benchmark settings -------------------------------#
alphas <- .1#c(0.1, 0.01)
nreps <- 500  #number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000

xi_maxs<- seq(3, 6, length=10)

sim_funs <- lapply(xi_maxs, function(x) wasserman_normal_sim_fun(20000,0.9,1, x))

fdr_methods <- list(ihw_binning_5fold, ihw_5fold)
fdr_methods <- lapply(fdr_methods, continuous_wrap)

results <- list()
for (i in 1:length(sim_funs)){
  res <- sim_fun_eval(sim_funs[[i]], fdr_methods, nreps = nreps, alphas, BiocParallel=F)
  results[[i]] <- res
}

final_output_sim3 <- do.call(rbind, results)

save(final_output_sim3, file = "final_results_sim_1_binning_rep500_sim_wass.RData")