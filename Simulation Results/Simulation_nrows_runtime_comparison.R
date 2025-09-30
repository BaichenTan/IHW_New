library(IHW)
library(microbenchmark)
source(file.path(getwd(), "/R/method_new_binning_admm.R"))
source(file.path(getwd(), "/R/ihw_convex.R"))
library(IHWpaper)
library(genefilter)


## ============ General benchmark settings ============
alphas <- 0.1
nreps  <- 1            # Monte Carlo replicates per (m, nbins)
bench_times <- 10      # repeats inside microbenchmark

## ============ Simulation function ===================
ms <- round(exp(seq(log(20000), log(1e6), length.out = 20)))
effect_size <- 5.333333
sim_funs <- lapply(ms, function(x) wasserman_normal_sim_fun(x, 0.9, 1, effect_size))

## ============ Helpers ===============================
get_rejections <- function(fit, alpha) {
  ap_fun <- try(utils::getFromNamespace("adj_pvalues", ns = "IHW"), silent = TRUE)
  if (!inherits(ap_fun, "try-error")) return(ap_fun(fit) < alpha)
  if (exists("adj_pvalues")) return(adj_pvalues(fit) < alpha)
  rh_fun <- try(utils::getFromNamespace("rejected_hypotheses", ns = "IHW"), silent = TRUE)
  if (!inherits(rh_fun, "try-error")) return(rh_fun(fit))
  if (!is.null(fit$rejections) && is.logical(fit$rejections)) return(fit$rejections)
  stop("Could not extract rejections from the fit object.")
}

compute_metrics <- function(H, rejected) {
  m <- length(H)
  rjs <- sum(rejected)
  false_rjs <- sum(H == 0 & rejected)
  true_alt  <- sum(H == 1)
  true_null <- sum(H == 0)
  c(
    rj_ratio = rjs / m,
    FDP      = ifelse(rjs == 0, 0, false_rjs / rjs),
    power    = ifelse(true_alt == 0, NA_real_, sum(H == 1 & rejected) / true_alt),
    FPR      = ifelse(true_null == 0, NA_real_, sum(H == 0 & rejected) / true_null),
    FWER     = as.numeric(false_rjs > 0)
  )
}

## ============ Main loops ============================
alpha_level <- 0.1
nbins_grid  <- c(300L, 800L, 1500L)

results_list <- list()   # per-(m, nbins, method) summary rows
timing_rows  <- list()   # long microbenchmark rows across all settings

for (i in seq_along(ms)) {
  test_sim <- sim_funs[[i]]
  for (g in nbins_grid) {
    # storage for this (m, nbins)
    m_binning <- matrix(NA_real_, nreps, 5)
    m_ihwold  <- matrix(NA_real_, nreps, 5)
    colnames(m_binning) <- colnames(m_ihwold) <- c("rj_ratio","FDP","power","FPR","FWER")
    
    rt_binning_med <- rt_binning_q25 <- rt_binning_q75 <- numeric(nreps)
    rt_ihwold_med  <- rt_ihwold_q25  <- rt_ihwold_q75  <- numeric(nreps)
    
    for (j in seq_len(nreps)) {
      sim_data <- test_sim(j)  # seed = j
      pvals <- sim_data$pvalue; filt <- sim_data$filterstat; H <- sim_data$H
      
      ## ---- Metrics (fit once) ----
      rej_bin <- tryCatch({
        fit_bin <- ihw_binning_admm(pvals, filt, alpha = alpha_level, nbins = g, nfolds = 5L)
        get_rejections(fit_bin, alpha_level)
      }, error = function(e) rep(FALSE, length(pvals)))
      m_binning[j, ] <- compute_metrics(H, rej_bin)
      
      rej_old <- tryCatch({
        fit_old <- IHWold:::ihwOld(pvals, filt, alpha = alpha_level, nbins = g,
                                   nfolds = 5L, nsplits_internal = 1L)
        get_rejections(fit_old, alpha_level)
      }, error = function(e) rep(FALSE, length(pvals)))
      m_ihwold[j, ] <- compute_metrics(H, rej_old)
      
      ## ---- Microbenchmark runtimes ----
      mb <- microbenchmark(
        IHW_OLD = IHWold:::ihwOld(
          pvals, filt, alpha = alpha_level, nbins = g, nfolds = 5L, nsplits_internal = 1L
        ),
        IHW_Binning = ihw_binning_admm(
          pvals, filt, alpha = alpha_level, nbins = g, nfolds = 5L
        ),
        times = bench_times, unit = "ns"
      )
      
      s <- as_tibble(mb) %>%
        group_by(expr) %>%
        summarise(
          median_sec = median(time) / 1e9,
          q25_sec    = quantile(time, 0.25) / 1e9,
          q75_sec    = quantile(time, 0.75) / 1e9,
          .groups = "drop"
        ) %>%
        mutate(iter = i, rep = j, nrow = ms[i], nbins = g)
      
      timing_rows[[length(timing_rows) + 1L]] <- s
      
      rt_binning_med[j] <- s$median_sec[s$expr == "IHW_Binning"]
      rt_binning_q25[j] <- s$q25_sec[s$expr == "IHW_Binning"]
      rt_binning_q75[j] <- s$q75_sec[s$expr == "IHW_Binning"]
      rt_ihwold_med[j]  <- s$median_sec[s$expr == "IHW_OLD"]
      rt_ihwold_q25[j]  <- s$q25_sec[s$expr == "IHW_OLD"]
      rt_ihwold_q75[j]  <- s$q75_sec[s$expr == "IHW_OLD"]
    } # end reps
    
    # per-(m, nbins) averages
    avg_binning <- colMeans(m_binning, na.rm = TRUE)
    avg_ihwold  <- colMeans(m_ihwold,  na.rm = TRUE)
    sd_binning  <- apply(m_binning, 2, sd, na.rm = TRUE)
    sd_ihwold   <- apply(m_ihwold,  2, sd, na.rm = TRUE)
    
    df_i <- rbind(
      data.frame(iter = i, nrows = ms[i], nbins = g, method = "IHW Binning",
                 rj_ratio = avg_binning[1], FDP = avg_binning[2], FDP_sd = sd_binning[2],
                 power = avg_binning[3], FPR = avg_binning[4], FWER = avg_binning[5],
                 runtime_median_sec = mean(rt_binning_med, na.rm = TRUE)),
      data.frame(iter = i, nrows = ms[i], nbins = g, method = "IHW OLD",
                 rj_ratio = avg_ihwold[1], FDP = avg_ihwold[2], FDP_sd = sd_ihwold[2],
                 power = avg_ihwold[3], FPR = avg_ihwold[4], FWER = avg_ihwold[5],
                 runtime_median_sec = mean(rt_ihwold_med,  na.rm = TRUE))
    )
    
    results_list[[length(results_list) + 1L]] <- df_i
    message(sprintf("Finished m=%d, nbins=%d", ms[i], g))
  } # end nbins
}   # end m

final_results <- bind_rows(results_list)
rownames(final_results) <- NULL

# Long timing table across all (m, nbins, reps, methods)
timing_long <- bind_rows(timing_rows) %>%
  mutate(method = as.character(expr)) %>%
  select(iter, rep, nrow, nbins, method, median_sec, q25_sec, q75_sec)

# Compact per-setting timing summary
timing_by_setting <- timing_long %>%
  group_by(iter, nrow, nbins, method) %>%
  summarise(
    median_sec_mean = mean(median_sec, na.rm = TRUE),
    q25_sec_mean    = mean(q25_sec,    na.rm = TRUE),
    q75_sec_mean    = mean(q75_sec,    na.rm = TRUE),
    .groups = "drop"
  )

## ============ Save ======================
save(final_results, timing_long, timing_by_setting,
     file = "nrows_nbins_final_results_microbench.RData")

# (If you also want CSVs:)
# write.csv(final_results, "final_results_metrics.csv", row.names = FALSE)
# write.csv(timing_long, "timing_long_microbench.csv", row.names = FALSE)
# write.csv(timing_by_setting, "timing_by_setting_microbench.csv", row.names = FALSE)