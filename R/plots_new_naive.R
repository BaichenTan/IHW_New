### This R file is the new plot weight function that is adapted to the new binning method
### I haven't fully figured out the adaptation without return_internal = T argument
### I use the bottomly dataset as an example
library("IHWpaper")
library(tidyverse)
library("DESeq2")
library("dplyr")
library(IHWold)
source(file.path(getwd(), "/R/method_new_binning_admm.R"))
source(file.path(getwd(), "/R/ihw_convex.R"))

bottomly <- analyze_dataset("bottomly")

ihw_res_binning <- ihw_binning_admm_new(bottomly$pvalue, bottomly$baseMean, 0.1, nbins = 50, nfolds_internal=4L, nfolds=5L, return_internal = T)

map_list <- ihw_res_binning$map_list
ws_mat <- ihw_res_binning$ws_mat
fold_lambda <- ihw_res_binning$fold_lambdas_idx
best_lambda_ave <- ihw_res_binning$lambda
best_lambda_weights <- ihw_res_binning$best_lambda_weights
rjs_mat_list <- ihw_res_binning$rjs_mat_list

# Prepare DataFrame for ggplot
df_plot <- data.frame()
for (fold in 1:(1+length(fold_lambda))) {
  if (fold != 1+length(fold_lambda)){
    lambda_idx <- fold_lambda[fold]
  } else{
    lambda_idx <- best_lambda_ave
  }
  
  map_mat <- map_list
  group_mapping <- map_mat[, lambda_idx]
  
  if (fold != 1+length(fold_lambda)){
    weights_vec <- ws_mat[[fold]][[lambda_idx]]
  } else{
    weights_vec <- best_lambda_weights
  }
  original_groups <- seq_along(group_mapping)
  new_weights <- sapply(group_mapping, function(g) weights_vec[g])
  
  df_plot <- rbind(df_plot, data.frame(
    fold = as.factor(fold),
    original_group = original_groups,
    weight = new_weights
  ))
}

# Plot
ggplot(df_plot, aes(x = original_group, y = weight, color = fold, group = fold)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Selected Weights by Original Group for Best Lambda per Fold",
    x = "Original Group ID",
    y = "Weight"
  ) +
  theme_minimal()