source(file.path(getwd(), "/IHWpaper-master/R/ddhw_wrappers.R"))
source(file.path(getwd(), "/IHWpaper-master/R/benchmarking.R"))
source(file.path(getwd(), "/R/method_new_binning_admm.R"))

# T-test informative Simulation Benchmark
#----------------- General benchmark settings -------------------------------#
alphas <- 0.1
nreps <- 1000 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000

#sim_funs <- list(null_sim_fun(ms))
eff_sizes <- seq(1, 2.5, length=20)
sim_funs <- lapply(eff_sizes, function(x) du_ttest_sim_fun(ms,0.95,x, uninformative_filter=F))

fdr_methods <- list(ihw_binning_5fold, ihw_5fold)
fdr_methods <- lapply(fdr_methods, continuous_wrap)

results <- list()
for (i in 1:length(sim_funs)){
  res <- sim_fun_eval(sim_funs[[i]], fdr_methods, nreps = nreps, alphas, BiocParallel=F)
  results[[i]] <- res
  print(paste0("finished round: ", i))
}

final_output_sim1 <- do.call(rbind, results)

save(final_output_sim1, file = "final_results_sim_1_binning_rep1000_sim_ttest.RData")

