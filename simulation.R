library(genefilter)

du_ttest_sim <- function(m, pi0, effect_size, n_samples=10, uninformative_filter=FALSE, seed=NULL){
  #if (!is.null(seed)) set.seed(seed) #don't do seed here
  m0 <- ceiling(m*pi0)
  false_nulls <- sample(1:m, m-m0)
  z_table <- matrix(rnorm(n_samples*m), ncol=n_samples)
  z_table[false_nulls,(n_samples/2+1):n_samples] <- matrix(rnorm(n_samples/2*(m-m0), effect_size), ncol=n_samples/2)
  H <- rep(0,m)
  H[false_nulls] <- 1
  ttests <- rowttests(z_table, factor(rep(1:2,each=n_samples/2)))
  sds <- rowSds(z_table)
  filter_pvalue <- 1-pchisq((n_samples-1)*sds^2,n_samples-1)
  if (uninformative_filter){
    filterstats <- runif(m)
  } else {
    filterstats <- sds
  }
  
  simDf <- data.frame(H=H, pvalue= ttests$p.value, sd= sds, filterstat=filterstats, 
                      filter_pvalue=filter_pvalue)
  return(simDf)
}

calculate_test_stats <- function(sim, fdr_method_result){
  rejected <- IHW::rejected_hypotheses(fdr_method_result)
  rjs <- sum(rejected)
  false_rjs <- sum(sim$H == 0 & rejected)
  rj_ratio <- rjs/nrow(sim)
  FDP <- ifelse(rjs == 0, 0, false_rjs/rjs)
  power <- ifelse(sum(sim$H) == 0, NA, sum(sim$H == 1 & rejected)/sum(sim$H ==1))
  # in principle I should take special care in case only alternatives, but we are not
  # interested in this scenario...
  FPR <-  sum(sim$H == 0 & rejected)/sum(sim$H == 0)
  FWER <- as.numeric(false_rjs > 0)
  df <- data.frame(rj_ratio=rj_ratio, FDP=FDP, power=power, FPR=FPR, FWER=FWER)
  df
}
