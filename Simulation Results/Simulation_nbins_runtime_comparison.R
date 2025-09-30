library(IHW)
library(IHWpaper)
library(microbenchmark)
source(file.path(getwd(), "/R/method_new_binning_admm.R"))
source(file.path(getwd(), "/R/ihw_convex.R"))
library(genefilter)


## ============ General benchmark settings ============
alphas <- 0.1
nreps  <- 1           # Monte Carlo replicates per setting
bench_times <- 10     # number of repeats for microbenchmark per fit

## ============ Simulation function ===================
ms <- 20000
xi_maxs <- seq(3, 6, length = 10)
sim_funs <- lapply(xi_maxs, function(x) wasserman_normal_sim_fun(ms, 0.9, 1, x))
sim_fun <- sim_funs[[8]]           # choose one effect-size setting
nbins   <- seq(30, 120, length = 6)

## ============ Helpers ===============================

# Extract rejections (logical) from a fit at level alpha
get_rejections <- function(fit, alpha) {
  ap_fun <- try(utils::getFromNamespace("adj_pvalues", ns = "IHW"), silent = TRUE)
  if (!inherits(ap_fun, "try-error")) return(ap_fun(fit) < alpha)
  
  if (exists("adj_pvalues")) return(adj_pvalues(fit) < alpha)
  
  rh_fun <- try(utils::getFromNamespace("rejected_hypotheses", ns = "IHW"), silent = TRUE)
  if (!inherits(rh_fun, "try-error")) return(rh_fun(fit))
  
  if (!is.null(fit$rejections) && is.logical(fit$rejections)) return(fit$rejections)
  
  stop("Could not extract rejections from the fit object.")
}

# Compute metrics given truth H and rejected (logical)
compute_metrics <- function(H, rejected) {
  m <- length(H)
  rjs <- sum(rejected)
  false_rjs <- sum(H == 0 & rejected)
  true_alt  <- sum(H == 1)
  true_null <- sum(H == 0)
  
  rj_ratio <- rjs / m
  FDP   <- ifelse(rjs == 0, 0, false_rjs / rjs)
  power <- ifelse(true_alt == 0, NA_real_, sum(H == 1 & rejected) / true_alt)
  FPR   <- ifelse(true_null == 0, NA_real_, sum(H == 0 & rejected) / true_null)
  FWER  <- as.numeric(false_rjs > 0)
  
  c(rj_ratio = rj_ratio, FDP = FDP, power = power, FPR = FPR, FWER = FWER)
}

# Convenience: summarize microbenchmark (ns) to seconds
mb_summary_sec <- function(mb_obj) {
  tibble(
    median_sec = median(mb_obj$time) / 1e9,
    q25_sec    = quantile(mb_obj$time, 0.25) / 1e9,
    q75_sec    = quantile(mb_obj$time, 0.75) / 1e9
  )
}

## ============ Main loop =============================
alpha_level <- 0.1
n_outer <- length(nbins)   # 6
n_inner <- nreps

# per-setting tables
results_list <- vector("list", length = n_outer)  # your metrics table per nbins
timing_rows  <- list()                            # long timing rows per replicate

for (i in seq_len(n_outer)) {
  test_sim <- sim_fun
  nbin     <- round(nbins[i])
  
  m_binning <- matrix(NA_real_, n_inner, 5)
  m_ihwold  <- matrix(NA_real_, n_inner, 5)
  colnames(m_binning) <- colnames(m_ihwold) <-
    c("rj_ratio", "FDP", "power", "FPR", "FWER")
  
  # store per-replicate timing summaries (optional)
  rt_binning_med <- numeric(n_inner)
  rt_binning_q25 <- numeric(n_inner)
  rt_binning_q75 <- numeric(n_inner)
  rt_ihwold_med  <- numeric(n_inner)
  rt_ihwold_q25  <- numeric(n_inner)
  rt_ihwold_q75  <- numeric(n_inner)
  
  for (j in seq_len(n_inner)) {
    sim_data <- test_sim(j)  # uses j as the seed
    
    pvals <- sim_data$pvalue
    filt  <- sim_data$filterstat
    H     <- sim_data$H
    
    ## --- (A) Fit once for each method to get metrics ----------------
    # ADMM binning
    rej_bin <- tryCatch({
      fit_bin <- ihw_binning_admm(
        pvals, filt,
        alpha = alpha_level,
        nbins = nbin,
        nfolds = 5L
      )
      get_rejections(fit_bin, alpha_level)
    }, error = function(e) rep(FALSE, length(pvals)))
    m_binning[j, ] <- compute_metrics(H, rej_bin)
    
    # IHWold
    rej_old <- tryCatch({
      fit_old <- IHWold:::ihwOld(
        pvals, filt,
        alpha = alpha_level,
        nbins  = nbin,
        nfolds = 5L,
        nsplits_internal = 1L
      )
      get_rejections(fit_old, alpha_level)
    }, error = function(e) rep(FALSE, length(pvals)))
    m_ihwold[j, ] <- compute_metrics(H, rej_old)
    
    ## --- (B) Microbenchmark runtimes (independent of metrics) -------
    mb <- microbenchmark(
      IHW_OLD = IHWold:::ihwOld(
        pvals, filt, alpha = alpha_level, nbins = nbin,
        nfolds = 5L,
        nsplits_internal = 1L
      ),
      IHW_Binning = ihw_binning_admm(
        pvals, filt, alpha = alpha_level, nbins = nbin, nfolds = 5L
      ),
      times = bench_times, unit = "ns"
    )
    
    # Summaries in seconds
    s <- mb %>% as_tibble() %>% group_by(expr) %>%
      summarise(median_sec = median(time)/1e9,
                q25_sec    = quantile(time, 0.25)/1e9,
                q75_sec    = quantile(time, 0.75)/1e9,
                .groups = "drop") %>%
      mutate(iter = i, rep = j, nbins = nbin)
    
    # store for long table
    timing_rows[[length(timing_rows) + 1L]] <- s
    
    # also keep per-rep medians for quick per-setting averages
    rt_binning_med[j] <- s$median_sec[s$expr == "IHW_Binning"]
    rt_binning_q25[j] <- s$q25_sec[s$expr == "IHW_Binning"]
    rt_binning_q75[j] <- s$q75_sec[s$expr == "IHW_Binning"]
    rt_ihwold_med[j]  <- s$median_sec[s$expr == "IHW_OLD"]
    rt_ihwold_q25[j]  <- s$q25_sec[s$expr == "IHW_OLD"]
    rt_ihwold_q75[j]  <- s$q75_sec[s$expr == "IHW_OLD"]
  }
  
  # per-setting averages (metrics)
  avg_binning <- colMeans(m_binning, na.rm = TRUE)
  avg_ihwold  <- colMeans(m_ihwold,  na.rm = TRUE)
  
  sd_binning <- apply(m_binning, 2, sd, na.rm = TRUE)
  sd_ihwold  <- apply(m_ihwold,  2, sd, na.rm = TRUE)
  
  # per-setting average timings (median-of-medians, etc.)
  rt_binning_mean_med <- mean(rt_binning_med, na.rm = TRUE)
  rt_ihwold_mean_med  <- mean(rt_ihwold_med,  na.rm = TRUE)
  
  df_i <- rbind(
    data.frame(iter = i, nbins = nbin, method = "IHW Binning",
               rj_ratio = avg_binning[1], FDP = avg_binning[2],
               FDP_sd = sd_binning[2], power = avg_binning[3],
               FPR = avg_binning[4], FWER = avg_binning[5],
               runtime_median_sec = rt_binning_mean_med),
    data.frame(iter = i, nbins = nbin, method = "IHW OLD",
               rj_ratio = avg_ihwold[1], FDP = avg_ihwold[2],
               FDP_sd = sd_ihwold[2], power = avg_ihwold[3],
               FPR = avg_ihwold[4], FWER = avg_ihwold[5],
               runtime_median_sec = rt_ihwold_mean_med)
  )
  
  results_list[[i]] <- df_i
  message(sprintf("Finished setting %d / %d (nbins = %d)", i, n_outer, nbin))
}

final_results <- do.call(rbind, results_list)
rownames(final_results) <- NULL

# Long timing table across all iters/reps/methods (seconds)
timing_long <- bind_rows(timing_rows) %>%
  mutate(expr = as.character(expr)) %>%
  dplyr::rename(method = expr)

# Optional: a compact per-setting timing summary (median of medians)
timing_by_setting <- timing_long %>%
  group_by(iter, nbins, method) %>%
  summarise(
    median_sec_mean = mean(median_sec, na.rm = TRUE),
    q25_sec_mean    = mean(q25_sec, na.rm = TRUE),
    q75_sec_mean    = mean(q75_sec, na.rm = TRUE),
    .groups = "drop"
  )

## ============ Save everything ======================
save(final_results, timing_long, timing_by_setting,
     file = "nbins_final_results_computational_time_microbench.RData")

# (If you also want CSVs:)
# write.csv(final_results, "final_results_metrics.csv", row.names = FALSE)
# write.csv(timing_long, "timing_long_microbench.csv", row.names = FALSE)
# write.csv(timing_by_setting, "timing_by_setting_microbench.csv", row.names = FALSE)