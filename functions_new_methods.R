##New Method
library(IHW)

find_lambda_max <- function(grenander_list, m_groups, alpha, t_0){
  t_inf <- IHW:::Rcpp_classicTautString_TV1(t_0, 1e10)
  density <- IHW:::GrenMix_pdf(grenander_list, m_groups, t_inf[1])
  mu_dual <- density / (1 - alpha * density)
  return(c(t_inf, IHW:::lambdaMax(grenander_list, m_groups, alpha, mu_dual, t_inf[1])))
}

optimal_ts_new <- function(grenander_list, m_groups, alpha, lambdas){
  #first obtain the unregularized rejection threholds
  t_0 <- IHW:::unregularized_thresholds_bh(grenander_list, m_groups, alpha = 0.1)
  
  #obtain the maximum lambda and the single t value
  t_lamda_inf <- find_lambda_max(grenander_list, m_groups, alpha, t_0)
  t_inf <- t_lamda_inf[1:length(t_lamda_inf)-1]
  lambda_max <- t_lamda_inf[length(t_lamda_inf)]
  
  #obtain the rest of the t values for various lamdbas
  lambdas <- lambda_max*lambdas
  print(lambdas)
  t_lam <- lapply(lambdas, function(x) IHW:::Rcpp_classicTautString_TV1(t_0, x))
  t_lam <- do.call(cbind, t_lam)
  return(cbind(t_inf, t_lam, t_0))
}
