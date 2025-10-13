#' @export
ihw <- function(...)
{
  UseMethod("ihw")
}


#' ihw: Main function for Independent Hypothesis Weighting
#'
#' Given a vector of p-values, a vector of covariates which are independent of the p-values under the null hypothesis and
#' a nominal significance level alpha, IHW learns multiple testing weights and then applies the weighted Benjamini Hochberg
#' (or Bonferroni) procedure.
#'
#'
#' @param pvalues  Numeric vector of unadjusted p-values.
#' @param covariates  Vector which contains the one-dimensional covariates (independent under the H0 of the p-value)
#'                for each test. Can be numeric or a factor. (If numeric it will be converted into factor by binning.)
#' @param alpha   Numeric, sets the nominal level for FDR control.
#' @param covariate_type  "ordinal" or "nominal" (i.e. whether covariates can be sorted in increasing order or not)
#' @param nbins  Integer, number of groups into which p-values will be split based on covariate. Use "auto" for
#'             automatic selection of the number of bins. Only applicable when covariates is not a factor.
#' @param m_groups Integer vector of length equal to the number of levels of the covariates (only to be specified
#'               when the latter is a factor/categorical). Each entry corresponds to the number of hypotheses to be tested in
#'               each group (stratum). This argument needs to be given when the complete vector of p-values is
#'               not available, but only p-values below a given threshold, for example because of memory reasons.
#'               See the vignette for additional details and an example of how this principle can be applied with
#'               numerical covariates.
#' @param folds  Integer vector or NULL. Pre-specify assignment of hypotheses into folds.
#' @param quiet  Boolean, if False a lot of messages are printed during the fitting stages.
#' @param nfolds Number of folds into which the p-values will be split for the pre-validation procedure
#' @param nfolds_internal  Within each fold, a second  (nested) layer of cross-validation can be conducted to choose a good
#'              regularization parameter. This parameter controls the number of nested folds.
#' @param lambdas  Numeric vector which defines the grid of possible regularization parameters.
#'				Use "auto" for automatic selection.
#' @param seed Integer or NULL. Split of hypotheses into folds is done randomly. To have the output of the function be reproducible,
#'	the seed of the random number generator is set to this value at the start of the function. Use NULL if you don't want to set the seed.
#' @param adjustment_type Character ("BH" or "bonferroni") depending on whether you want to control FDR or FWER.
#' @param null_proportion Boolean, if True (default is False), a modified version of Storey's estimator
#'        is used within each bin to estimate the proportion of null hypotheses.
#' @param null_proportion_level Numeric, threshold for Storey's pi0 estimation procedure, defaults to 0.5
#' @param return_internal Returns a lower level representation of the output (only useful for debugging purposes).
#' @param ... Arguments passed to internal functions.
#'
#' @return A ihwResult object.
#' @seealso ihwResult, plot,ihwResult-method, ihw.DESeqResults
#'
#' @examples
#'
#' save.seed <- .Random.seed; set.seed(1)
#' X   <- runif(20000, min=0, max=2.5)   # covariate
#' H   <- rbinom(20000,1,0.1)            # hypothesis true or false
#' Z   <- rnorm(20000, H*X)              # Z-score
#' .Random.seed <- save.seed
#' pvalue <- 1-pnorm(Z)                  # pvalue
#'
#' ihw_fdr <- ihw(pvalue, X, .1)        # Standard IHW for FDR control
#' ihw_fwer <- ihw(pvalue, X, .1, adjustment_type = "bonferroni")    # FWER control
#
#' table(H[adj_pvalues(ihw_fdr) <= 0.1] == 0) #how many false rejections?
#' table(H[adj_pvalues(ihw_fwer) <= 0.1] == 0)
#'
#'
#' @export
#' @aliases ihw
#'
ihw.default <- function(pvalues,
                        covariates,
                        alpha,
                        covariate_type = "ordinal",
                        nbins = "auto",
                        m_groups = NULL,
                        folds = NULL,
                        quiet = TRUE ,
                        nfolds = 5L,
                        admm = T,
                        order = F,
                        nfolds_internal = nfolds - 1L,
                        lambdas = "auto",
                        seed = 1L,
                        adjustment_type = "BH",
                        null_proportion = FALSE,
                        null_proportion_level = 0.5,
                        return_internal = FALSE,
                        ...) {
  # This function essentially wraps the lower level function ihw_internal
  # e.g. takes care of NAs, sorts pvalues and then
  # returns a nice ihw object
  # covariate_type = "ordinal" hasn't been implemented yet
  
  
  
  if (!adjustment_type %in% c("BH", "bonferroni")) {
    stop("IHW currently only works with BH or bonferroni types of multiple testing corrections")
  }
  
  
  if (is.null(folds)) {
    nfolds <- as.integer(nfolds)
  } else {
    nfolds <- length(unique(folds))
  }
  
  nfolds_internal <- as.integer(nfolds_internal)
  
  
  if (nfolds == 1) {
    message("Using only 1 fold! Only use this if you want to learn the weights, but NEVER for testing!")
  }
  
  
  # gracefully handle NAs
  nna <- !is.na(pvalues)
  
  if (all(nna)) {
    nna <- TRUE
  }
  
  # filter our p-values
  pvalues <- pvalues[nna]
  covariates <- covariates[nna]
  weights <- rep(NA, length(pvalues))
  
  if (any(is.na(covariates))) {
    stop("Covariates corresponding to non-NA p-values should never be NA. Aborting.")
  }
  
  
  if (covariate_type == "ordinal" & is.numeric(covariates)) {
    if (!is.null(m_groups)) {
      stop("m_groups should only be specified when the covariates are a factor")
    }
    
    if (nbins == "auto") {
      nbins <-
        max(1, floor(length(pvalues) /500)) # rule of thumb..
    }
    groups <-
      as.factor(groups_by_filter(covariates, nbins, seed = seed))
    penalty <- "total variation"
    
  } else if (is.factor(covariates)) {
    groups <- droplevels(covariates)
    if (nbins != "auto" & nbins != nlevels(groups)) {
      warning(
        "Overwriting manually specified nbins, since it has to equal the number
				of levels for categorical covariates"
      )
    }
    
    nbins <- nlevels(groups)
    
    if (covariate_type == "nominal") {
      penalty <- "uniform deviation"
    } else if (covariate_type == "ordinal") {
      penalty <- "total variation"
    }
  } else {
    stop("Covariates are not of the appropriate type")
  }
  
  if ((length(lambdas) == 1) & (lambdas[1] == "auto")) {
    lambdas <- c(Inf, 1 / 2^(1:7), 0)
  }
  
  if (nbins < 1) {
    stop("Cannot have less than one bin.")
  }
  
  
  group_levels <- levels(droplevels(groups))
  # sort pvalues globally
  order_pvalues <-
    order(pvalues) #TODO BREAK ties by filter statistic rank
  reorder_pvalues <- order(order_pvalues)
  
  sorted_groups <- groups[order_pvalues]
  sorted_pvalues <- pvalues[order_pvalues]
  
  if (!is.null(folds)) {
    sorted_folds <- as.factor(folds[order_pvalues])
  } else {
    sorted_folds <- NULL
  }
  
  if (is.null(m_groups)) {
    if (is.null(folds)) {
      m_groups <- table(sorted_groups)
    } else {
      m_groups <- table(sorted_groups, sorted_folds)
    }
  } else {
    # TODO: also check if m_groups entered by user is plausible based on observed m_groups and max p_value
  }
  
  #--- check if m_groups is correctly specified----------------------
  if (is.null(folds)) {
    if (length(group_levels) != length(m_groups)) {
      stop("Number of unique levels should be equal to length of m_groups.")
    }
  } else {
    if (any(dim(m_groups) != c(length(group_levels), nfolds))) {
      stop("m_groups should be of dimension nbins*nfolds")
    }
  }
  # once we have groups, check whether they include enough p-values
  if (nbins > 1 & any(m_groups < 200)) {
    message(
      "We recommend that you supply (many) more than 200 p-values for meaningful data-driven hypothesis weighting results."
    )
  }
  
  # start by making sure seed is correctly specified
  if (!is.null(seed)) {
    #http://stackoverflow.com/questions/14324096/setting-seed-locally-not-globally-in-r?rq=1
    tmp <- runif(1)
    old <- .Random.seed
    on.exit({
      .Random.seed <<- old
    })
    set.seed(as.integer(seed)) #seed
  }
  
  if (nbins == 1) {
    nfolds <- 1L
    nfolds_internal <- 1L
    message("Only 1 bin; IHW reduces to Benjamini Hochberg (uniform weights)")
    sorted_adj_p <-
      p.adjust(sorted_pvalues,
               method = adjustment_type,
               n = sum(m_groups))
    rjs <- sum(sorted_adj_p <= alpha)
    res <-
      list(
        lambda = 0,
        fold_lambdas = 0,
        rjs = rjs,
        sorted_pvalues = sorted_pvalues,
        sorted_weighted_pvalues = sorted_pvalues,
        sorted_adj_p = sorted_adj_p,
        sorted_weights = rep(1, length(sorted_pvalues)),
        sorted_groups = as.factor(rep(1, length(sorted_pvalues))),
        sorted_folds = as.factor(rep(1, length(sorted_pvalues))),
        weight_matrix = matrix(1)
      )
  } else if (nbins > 1) {
    res <- ihw_internal(
      sorted_groups,
      sorted_pvalues,
      alpha,
      lambdas,
      m_groups,
      order,
      penalty = penalty,
      quiet = quiet,
      sorted_folds = sorted_folds,
      nfolds = nfolds,
      admm = admm,
      nfolds_internal = nfolds_internal,
      seed = NULL,
      adjustment_type = adjustment_type,
      null_proportion = null_proportion,
      null_proportion_level = null_proportion_level,
      ...
    )
  }
  
  if (return_internal) {
    return(res)
  }
  
  # resort back to original form
  weights <-
    fill_nas_reorder(res$sorted_weights, nna, reorder_pvalues)
  pvalues <-
    fill_nas_reorder(res$sorted_pvalues, nna, reorder_pvalues)
  weighted_pvalues <-
    fill_nas_reorder(res$sorted_weighted_pvalues, nna, reorder_pvalues)
  adj_pvalues <-
    fill_nas_reorder(res$sorted_adj_p, nna, reorder_pvalues)
  groups     <-
    factor(fill_nas_reorder(res$sorted_groups, nna, reorder_pvalues),
           levels = group_levels)
  fold_levels <- levels(res$sorted_folds)
  folds      <-
    factor(fill_nas_reorder(res$sorted_folds, nna, reorder_pvalues),
           levels = fold_levels)
  covariates <-
    fill_nas_reorder(covariates, nna, 1:length(covariates))
  
  df <- data.frame(
    pvalue = pvalues,
    adj_pvalue = adj_pvalues,
    weight = weights,
    weighted_pvalue = weighted_pvalues,
    group = groups,
    covariate = covariates,
    fold = as.factor(folds)
  )
  
  ihw_obj <- new(
    "ihwResult",
    df = df,
    weights = res$weight_matrix,
    alpha = alpha,
    nbins = as.integer(nbins),
    nfolds = nfolds,
    regularization_term = res$fold_lambdas,
    m_groups = as.integer(m_groups),
    penalty = penalty,
    covariate_type = covariate_type,
    adjustment_type = adjustment_type,
    reg_path_information = data.frame(),
    solver_information = list()
  )
  
  ihw_obj
  
}

# operate on pvalues which have already been sorted and without any NAs
ihw_internal <-
  function(sorted_groups,
           sorted_pvalues,
           alpha,
           lambdas,
           m_groups,
           penalty = "total variation",
           return_num_rejection_list = FALSE,
           inner_cross = FALSE,
           seed = NULL,
           quiet = TRUE,
           sorted_folds = NULL,
           nfolds = 10L,
           admm = T,
           new_method = F,
           nfolds_internal = nfolds,
           adjustment_type = "BH",
           null_proportion = FALSE,
           null_proportion_level = 0.5,
           debug_flag = FALSE,
           ...) {
  
    folds_prespecified <- !is.null(sorted_folds)
    n_lambdas <- length(lambdas)
    m <- length(sorted_pvalues)
    split_sorted_pvalues <- split(sorted_pvalues, sorted_groups)
    m_groups_available <- sapply(split_sorted_pvalues, length)
    
    
    #  do the k-fold strategy
    if (!is.null(seed)) {
      #http://stackoverflow.com/questions/14324096/setting-seed-locally-not-globally-in-r?rq=1
      old <- .Random.seed
      on.exit({
        .Random.seed <<- old
      })
      set.seed(as.integer(seed)) #seed
    }
    
    if (!folds_prespecified) {
      sorted_folds <-
        factor(sample(1:nfolds, m, replace = TRUE), levels = 1:nfolds)
    }
    
    sorted_weights <- rep(NA, m)
    fold_lambdas <- rep(NA, nfolds)
    
    weight_matrix <- matrix(NA, nlevels(sorted_groups), nfolds)
    
    ws_mat_list <- list()
    rjs_mat_list <- list()
    
    original_sorted_groups <- list()
    map_list <- list()
    initial_weight_list <- list()
    regroup_list <- list()
    

    for (i in seq_along(levels(sorted_folds))) {
      if (!quiet)
        message(paste("Estimating weights for fold", i))
      if (nfolds == 1) {
        filtered_sorted_groups <- sorted_groups
        filtered_sorted_pvalues <- sorted_pvalues
        filtered_split_sorted_pvalues <-
          split(filtered_sorted_pvalues, filtered_sorted_groups)
        
        m_groups_holdout_fold <- m_groups
        m_groups_other_folds <- m_groups
        m_groups_other_folds_uncollapsed <- m_groups_other_folds
      } else {
        fold_i = levels(sorted_folds)[i]
        filtered_sorted_groups <- sorted_groups[sorted_folds != fold_i]
        filtered_sorted_pvalues <- sorted_pvalues[sorted_folds != fold_i]
        filtered_split_sorted_pvalues <-
          split(filtered_sorted_pvalues, filtered_sorted_groups)
        
        #record the original groups for each fold
        original_sorted_groups[[i]] <- sorted_groups[sorted_folds == fold_i]
        
        if (!folds_prespecified || length(dim(m_groups)) == 1) {
          # TODO: add check to see if for common use case the first term is 0
          m_groups_holdout_fold <- (m_groups - m_groups_available) / nfolds +
            m_groups_available - sapply(filtered_split_sorted_pvalues, length)
          m_groups_other_folds <- m_groups - m_groups_holdout_fold
          m_groups_other_folds_uncollapsed <- m_groups_other_folds
        } else {
          m_groups_holdout_fold <- m_groups[, colnames(m_groups) == fold_i]
          m_groups_other_folds_uncollapsed <-
            m_groups[, colnames(m_groups) != fold_i, drop = FALSE]
          m_groups_other_folds <-
            rowSums(m_groups_other_folds_uncollapsed)
        }
      }
    

      #so far we have successfully record the group information for the hold-out fold and the other folds.
      #Now for each fold we calculate the admm path t_g(\lambda) for each lambda
        
      if (!return_num_rejection_list){
        initial_weights <-
        ihw_convex_weights(
          filtered_split_sorted_pvalues,
          alpha,
          m_groups_holdout_fold,
          m_groups_other_folds,
          lambdas,
          admm,
          penalty = penalty,
          adjustment_type = adjustment_type,
          quiet = quiet,
          ...
        )

        initial_weight_list[[i]] <- initial_weights
        
        #Now for each lambda, we start the groupping
        regrouping_res <- iter_fused_regroup_admm(filtered_sorted_pvalues, filtered_sorted_groups, lambdas, initial_weights)
        regroup_pvalue_groups <- regrouping_res$groups
        map_new_old_group <- regrouping_res$group_id_matrix
        
        regroup_list[[i]] <- regroup_pvalue_groups
        map_list[[i]] <- map_new_old_group
        
        #build, for each lambda, the TRUE counts per new bin by aggregating the
        #original counts of the "other folds"
        m_groups_regrouped_list_other <- vector("list", length(lambdas))
        m_groups_regrouped_list_holdout <- vector("list", length(lambdas))
        
        for (l in seq_along(lambdas)) {
          map_l <- as.integer(map_new_old_group[, l])   # old_bin_id -> new_bin_id
          # aggregate original counts to regrouped bins
          # (tapply returns a named vector with names = new_bin_id)
          mg_holdout_new <- tapply(m_groups_holdout_fold, INDEX = map_l, FUN = sum)
          mg_other_new <- tapply(m_groups_other_folds, INDEX = map_l, FUN = sum)
          # be explicit about names so compute_BH_* can align by factor levels
          names(mg_other_new) <- as.character(names(mg_other_new))
          names(mg_holdout_new) <- as.character(names(mg_holdout_new))
          m_groups_regrouped_list_holdout[[l]] <- mg_holdout_new
          m_groups_regrouped_list_other[[l]] <- mg_other_new
        }
        
        # call the helper WITH the regrouped counts
        thres_weight_res <- compute_BH_thresholds_and_weights(
          filtered_sorted_pvalues,
          regroup_pvalue_groups,
          lambdas,
          alpha = alpha,
          grenander_binsize = 1,
          quiet = TRUE,
          m_groups_regrouped_list_holdout <- m_groups_regrouped_list_holdout,
          m_groups_regrouped_list_other = m_groups_regrouped_list_other  # <<< pass TRUE counts here
        )
        
        ws_mat_list[[i]] <- lapply(thres_weight_res, `[[`, "weights")


      }
      
      #only for inner cross validation, we apply cheap cross validation through unregularized BH 
      #and using that to construct "cheap weights" in the second layer of cross validation
      if (inner_cross){
        bh_result <- unregularized_BH(
          split_sorted_pvalues = filtered_split_sorted_pvalues,
          m_groups = m_groups_other_folds_uncollapsed,
          alpha = alpha,
          grenander_binsize = 1,
          quiet = quiet
        )
        weights <- thresholds_to_weights(as.vector(bh_result$t), m_groups_holdout_fold)
        ws_mat_list[[i]] <- weights
        
      }
      
      
      if (!return_num_rejection_list && length(lambdas) > 1) {
          
          nsplits_internal <- 1L
          rjs_mat  <- matrix(NA, nrow = n_lambdas, ncol = nsplits_internal)
          
          for (s in seq_len(nsplits_internal)) {
            if (nfolds > 2) {
              sorted_folds_internal = droplevels(sorted_folds[sorted_folds != fold_i])
            } else {
              sorted_folds_internal = NULL
            }
            
            rjs_lam_list <- rep(NA, n_lambdas)
            for (l in 1:length(lambdas)){
      
              filtered_sorted_groups_internal <- regroup_pvalue_groups[[l]]
              
              map_l <- as.integer(map_new_old_group[, l])   # old_bin_id -> new_bin_id
              # aggregate original counts to regrouped bins
              # (tapply returns a named vector with names = new_bin_id)
              mg_other_new <- tapply(m_groups_other_folds, INDEX = map_l, FUN = sum)
              # be explicit about names so compute_BH_* can align by factor levels
              names(mg_other_new) <- as.character(names(mg_other_new))
              m_groups_internal  <- mg_other_new
              
              if (nfolds_internal == 0){
                nfolds_internal <- 1
              }
              
              
              rjs_int <-
                ihw_internal(
                  filtered_sorted_groups_internal,
                  filtered_sorted_pvalues,
                  alpha,
                  lambdas,
                  m_groups_internal,
                  return_num_rejection_list = TRUE,
                  inner_cross = TRUE,
                  seed = NULL,
                  quiet = quiet,
                  sorted_folds = sorted_folds_internal,
                  nfolds = nfolds_internal,
                  admm = admm,
                  new_method = new_method,
                  adjustment_type = adjustment_type,
                  ...
                )

              rjs_lam_list[l] <- rjs_int
            }
          rjs_mat[, s] <- rjs_lam_list
        }
        rjs_mat_list[[i]] <- rjs_mat
      }
    }

    #if just return the rejection thresholds without inner cross
    if (return_num_rejection_list && !inner_cross){
      rjs_lam_list <- rep(NA, n_lambdas)
      for (l in 1:length(lambdas)){
        for (i in seq_along(levels(sorted_folds))) {
          fold_i <- levels(sorted_folds)[i]
          group_ids <- sorted_groups[sorted_folds == fold_i]
          lambda_weights <- ws_mat_list[[i]][[l]]
          sorted_weights[sorted_folds == fold_i] <- lambda_weights[group_ids]
        }
        
        sorted_weighted_pvalues <-
          mydiv(sorted_pvalues, sorted_weights)
        sorted_adj_p <-
          p.adjust(sorted_weighted_pvalues,
                   method = adjustment_type,
                   n = sum(m_groups))
        rjs  <- sum(sorted_adj_p <= alpha)
        rjs_lam_list[l] <- rjs_int
      }
      
      return(rjs_lam_list)
    }
    
    
    #if returning rejection threshold with inner cross
    if (return_num_rejection_list && inner_cross) {
      for (i in seq_along(levels(sorted_folds))) {
        fold_i <- levels(sorted_folds)[i]
        group_ids <- sorted_groups[sorted_folds == fold_i]
        lambda_weights <- ws_mat_list[[i]]
        sorted_weights[sorted_folds == fold_i] <- lambda_weights[group_ids]
      }
      
      sorted_weighted_pvalues <-
        mydiv(sorted_pvalues, sorted_weights)
      sorted_adj_p <-
        p.adjust(sorted_weighted_pvalues,
                 method = adjustment_type,
                 n = sum(m_groups))
      rjs  <- sum(sorted_adj_p <= alpha)
    
      return(rjs)
    }
    
    #create a dataset that saves the groupping result for each fold
    new_sorted_groups <- list()
    old_weight_list <- list()
    fold_lambdas_idx <- c()
    
    
    if (!return_num_rejection_list) {
      for (i in seq_along(levels(sorted_folds))) {
        fold_i = levels(sorted_folds)[i]
        if (length(lambdas) > 1) {
          #threshold_per_fold <- ceiling(0.98*max(rjs_mat_list[[i]]))
          #lambda_idx_candidates_per_fold <- which(rjs_mat_list[[i]]>= threshold_per_fold)
          
          ## set rejection threshold and select lambda max
          threshold_per_fold <- ceiling(0.98*(max(rjs_mat_list[[i]]) - rjs_mat_list[[i]][1, 1]))
          lambda_idx_candidates_per_fold <- which(rjs_mat_list[[i]] - rjs_mat_list[[i]][1, 1] >= threshold_per_fold)
          
          lambda_idx_per_fold <- min(lambda_idx_candidates_per_fold)
          fold_lambdas[i] <- lambdas[lambda_idx_per_fold]
          
          
          #now assign group ids
          group_ids_original <- sorted_groups[sorted_folds == fold_i]
          map_list_fold_i <- map_list[[i]][, lambda_idx_per_fold]
          group_ids <- map_list_fold_i[group_ids_original]
          
          lambda_weights <- ws_mat_list[[i]][[lambda_idx_per_fold]]
          weight_matrix[, i] <- lambda_weights[map_list_fold_i]
          
          sorted_weights[sorted_folds == fold_i] <- lambda_weights[group_ids]
          
          new_sorted_groups[[i]] <- group_ids
          
          old_weight_list[[i]] <- initial_weight_list[[i]][, lambda_idx_per_fold]
          
          fold_lambdas_idx[i] <- lambda_idx_per_fold
          
          
        } else {
          lambda_idx <- 1
        }
        
      }
      

      sorted_weighted_pvalues <-
        mydiv(sorted_pvalues, sorted_weights)
      sorted_adj_p <-
        p.adjust(sorted_weighted_pvalues,
                 method = adjustment_type,
                 n = sum(m_groups))
      
      rjs   <- sum(sorted_adj_p <= alpha)
      
      

      lst <-
        list(
          lambda = 0.0,
          fold_lambdas = fold_lambdas,
          fold_lambdas_idx = fold_lambdas_idx,
          rjs = rjs,
          sorted_pvalues = sorted_pvalues,
          sorted_weighted_pvalues = sorted_weighted_pvalues,
          sorted_adj_p = sorted_adj_p,
          sorted_weights = sorted_weights,
          sorted_groups = sorted_groups,
          sorted_folds = sorted_folds,
          weight_matrix = weight_matrix,
          original_sorted_groups = original_sorted_groups,
          new_sorted_groups = new_sorted_groups,
          old_weight_list = old_weight_list,
          map_list = map_list,
          ws_mat = ws_mat_list,
          rjs_mat_list = rjs_mat_list
        )
      return(lst)
    }
    
  }

#' @rdname ihw.default
#' @param formula \code{\link{formula}}, specified in the form pvalue~covariate (only 1D covariate supported)
#' @param data data.frame from which the variables in formula should be taken
#' @export
ihw.formula <- function(formula, data = parent.frame(), ...) {
  if (length(formula) != 3) {
    stop("expecting formula of the form pvalue~covariate")
  }
  pv_name <- formula[[2]]
  cov_name <- formula[[3]]
  pvalues <- eval(pv_name, data, parent.frame())
  covariates <- eval(cov_name, data, parent.frame())
  ihw(pvalues, covariates, ...)
}


#' ihw.DESeqResults: IHW method dispatching on DESeqResults objects
#'
#' @param deseq_res "DESeqResults" object
#' @param filter Vector of length equal to number of rows of deseq_res object. This is used
#'         for the covariates in the call to ihw. Can also be a character,
#'         in which case deseq_res[[filter]] is used as the covariate
#' @param alpha   Numeric, sets the nominal level for FDR control.
#' @param adjustment_type Character ("BH" or "bonferroni") depending on whether you want to control FDR or FWER.
#' @param ... Other optional keyword arguments passed to ihw.
#'
#' @return A "DESeqResults" object, which includes weights and adjusted p-values returned
#'         	by IHW. In addition, includes a metadata slot with an "ihwResult" object.
#' @seealso ihw, ihwResult
#'
#' @examples \dontrun{
#'    library("DESeq2")
#'    library("airway")
#'    data("airway")
#'    dds <- DESeqDataSet(se = airway, design = ~ cell + dex)
#'    dds <- DESeq(dds)
#'    deseq_res <- results(dds)
#'    deseq_res <- ihw(deseq_res, alpha=0.1)
#'    #equivalent: deseq_res2 <- results(dds, filterFun = ihw)
#' }
#'
#' @export
#'
ihw.DESeqResults <- function(deseq_res,
                             filter = "baseMean",
                             alpha = 0.1,
                             adjustment_type = "BH",
                             ...) {
  if (missing(filter)) {
    filter <- deseq_res$baseMean
  } else if (is.character(filter)) {
    stopifnot(filter %in% colnames(deseq_res))
    filter <- deseq_res[[filter]]
  }
  
  stopifnot(length(filter) == nrow(deseq_res))
  
  ihw_res <- ihw(deseq_res$pvalue,
                 filter,
                 alpha = alpha,
                 adjustment_type = adjustment_type,
                 ...)
  
  deseq_res$padj <- adj_pvalues(ihw_res)
  deseq_res$weight <- weights(ihw_res)
  
  mcols(deseq_res)$type[names(deseq_res) == "padj"] <- "results"
  mcols(deseq_res)$description[names(deseq_res) == "padj"] <-
    paste("Weighted", adjustment_type, "adjusted p-values")
  
  mcols(deseq_res)$type[names(deseq_res) == "weight"] <- "results"
  mcols(deseq_res)$description[names(deseq_res) == "weight"] <-
    "IHW weights"
  
  metadata(deseq_res)[["alpha"]] <- alpha
  metadata(deseq_res)[["ihwResult"]] <- ihw_res
  deseq_res
}


iter_fused_regroup_admm <- function(p,                # vector of p‑values
                                    g,                # initial group labels
                                    lambdas,          # numeric vector of λ’s
                                    weight_matrix,    # matrix of group weights: (original groups) x (lambdas)
                                    eps = 0.001) {
  # --- Sanity checks ---
  if (length(p) != length(g)) {
    stop("Length of p-values and group labels (g) must be equal.")
  }
  if (ncol(weight_matrix) != length(lambdas)) {
    stop("Number of columns in weight_matrix must match number of lambdas.")
  }
  
  m <- nrow(weight_matrix)     # number of original groups
  num_lambda <- length(lambdas)
  g <- as.integer(g)
  
  # List of regrouped group labels per lambda (length = number of lambdas)
  regrouped_pvalue_groups <- vector("list", length = num_lambda)
  
  # Matrix of new group ids for original groups (m x p)
  group_id_matrix <- matrix(NA_integer_, nrow = m, ncol = num_lambda)
  
  for (j in seq_len(num_lambda)) {
    wj <- weight_matrix[, j]  # weights for λ_j
    
    # Start with group_labels 1 through m
    group_labels <- seq_len(m)
    
    # Merge adjacent groups if weights are within eps
    k <- 1
    while (k < m) {
      if (abs(wj[k + 1] - wj[k]) <= eps) {
        group_labels[group_labels == group_labels[k + 1]] <- group_labels[k]
      }
      k <- k + 1
    }
    
    # Reindex group labels to consecutive integers: 1, 2, ...
    group_labels <- as.integer(factor(group_labels))
    
    # Apply to p-values: assign each p-value to its new group
    regrouped_pvalue_groups[[j]] <- group_labels[g]
    
    # Save mapping for original groups
    group_id_matrix[, j] <- group_labels
  }
  
  return(list(
    groups = regrouped_pvalue_groups,     # each p-value’s new group ID per lambda
    map = do.call(cbind, regrouped_pvalue_groups),  # n x p matrix (n = number of p-values)
    group_id_matrix = group_id_matrix     # m x p matrix (m = original groups)
  ))
}


compute_BH_thresholds_and_weights <- function(
    p,
    regroup_pvalue_groups,          # list of length lambdas; each is a vector of new bin labels
    lambdas,
    alpha = 0.1,
    grenander_binsize = 1,
    quiet = TRUE,
    m_groups_regrouped_list_holdout, #list of vectors with TRUE counts per new bin in the holdout sets
    m_groups_regrouped_list_other  #list of vectors with TRUE counts per new bin in the other sets
) {
  stopifnot(is.numeric(p), is.list(regroup_pvalue_groups))
  stopifnot(length(regroup_pvalue_groups) == length(lambdas))
  
  results <- vector("list", length(lambdas))
  
  for (j in seq_along(lambdas)) {
    # ensure we split with a stable bin ordering
    g_lab  <- regroup_pvalue_groups[[j]]
    g_fac  <- factor(g_lab, levels = sort(unique(g_lab)))
    split_p <- split(p, g_fac)
    lvls    <- levels(g_fac)
    
    # --- crucial: use provided TRUE counts after regrouping, else fall back to observed counts
    mg_in <- m_groups_regrouped_list_other[[j]]
    mg_in_holdout <- m_groups_regrouped_list_holdout[[j]]
    # allow both named/unnamed input; reorder by current levels
    if (is.null(names(mg_in))) {
      # assume bins are 1..K; coerce to character to index by 'lvls'
      names(mg_in) <- as.character(seq_along(mg_in))
      names(mg_in_holdout) <- as.character(seq_along(mg_in_holdout))
    }
    m_groups <- as.integer(mg_in[lvls])
    m_groups_holdout <- as.integer(mg_in_holdout[lvls])
    if (any(is.na(m_groups))) {
      stop("m_groups_regrouped_list[[", j, "]] does not match regrouped bin labels.")
    }
    
    
    # unregularized BH using TRUE counts for these new bins
    #print(m_groups)
    #print("hehehehehe")
    bh_result <- unregularized_BH(
      split_sorted_pvalues = split_p,
      m_groups = m_groups,
      alpha = alpha,
      quiet = quiet,
      grenander_binsize = grenander_binsize
    )
    
    # thresholds -> weights (same semantics as IHW)
    weights <- IHW:::thresholds_to_weights(as.vector(bh_result$t), m_groups_holdout)
    
    results[[j]] <- list(t = as.vector(bh_result$t),
                         weights = weights,
                         m_groups = m_groups_holdout)
  }
  results
}


##helper function that obtains the unregularized BH to start the fused regrouping 
unregularized_BH <- function(split_sorted_pvalues, 
                             m_groups,
                             alpha = 0.1,
                             quiet = T,
                             grenander_binsize = 1) {
  nbins <- length(split_sorted_pvalues)
  
  # also there temporarily to help with numerical stability of downstream optimization
  max_pval <- 0.9 * max(sapply(split_sorted_pvalues, max))
  if (max_pval < 0.1) {
    split_sorted_pvalues <-
      lapply(split_sorted_pvalues, function(x)
        x / max_pval)
    alpha <- alpha / max_pval
  }
  
  # preprocessing:  Set very low p-values to 10^(-20) to help with numerical stability
  # Note: This only affects the internal optimization, the higher level functions will
  # still return adjusted pvalues based on the original p-values
  split_sorted_pvalues <-
    lapply(split_sorted_pvalues, function(x)
      ifelse(x > 10 ^ (-20), x, 10 ^ (-20)))
  m <- sum(m_groups)
  
  if (nbins != length(m_groups)) {
    stop("length of m_groups should be equal to number of bins")
  }
  
  #lapply grenander...
  if (!quiet)
    message("Applying Grenander estimator within each bin.")
  
  grenander_list <- mapply(
    presorted_grenander,
    split_sorted_pvalues,
    m_groups,
    grenander_binsize = grenander_binsize,
    quiet = quiet,
    SIMPLIFY = FALSE
  )
  
  
  ts <-  unregularized_thresholds_bh(grenander_list, m_groups, alpha)
  ts <- matrix(ts, ncol = 1)
  
  return(list(t = ts))
  
}

#from pvalues to weights 
ihw_convex_weights_unregularized <-
  function(split_sorted_pvalues,
           alpha,
           m_groups,
           m_groups_grenander,
           penalty = "total variation",
           adjustment_type = "BH",
           grenander_binsize = 1,
           quiet = quiet) {
    # m_groups used for balancing (weight budget)
    # m_groups_grenander (used for grenander estimator)
    # Practically always it will hold: m_groups \approx m_groups_grenander
    nbins <- length(split_sorted_pvalues)
    
    # also there temporarily to help with numerical stability of downstream optimization
    max_pval <- 0.9 * max(sapply(split_sorted_pvalues, max))
    if (max_pval < 0.1) {
      split_sorted_pvalues <-
        lapply(split_sorted_pvalues, function(x)
          x / max_pval)
      alpha <- alpha / max_pval
    }
    # preprocessing:  Set very low p-values to 10^(-20) to help with numerical stability
    # Note: This only affects the internal optimization, the higher level functions will
    # still return adjusted pvalues based on the original p-values
    split_sorted_pvalues <-
      lapply(split_sorted_pvalues, function(x)
        ifelse(x > 10 ^ (-20), x, 10 ^ (-20)))
    m <- sum(m_groups)
    
    if (nbins != length(m_groups)) {
      stop("length of m_groups should be equal to number of bins")
    }
    
    #lapply grenander...
    if (!quiet)
      message("Applying Grenander estimator within each bin.")
    
    grenander_list <- mapply(
      presorted_grenander,
      split_sorted_pvalues,
      m_groups_grenander,
      grenander_binsize = grenander_binsize,
      quiet = quiet,
      SIMPLIFY = FALSE
    )
    
    ts <-  unregularized_thresholds_bh(grenander_list, m_groups, alpha)
    ts <- matrix(ts, ncol = 1)
    
    ws <- apply(ts, 2, function(x)
      thresholds_to_weights(x, m_groups))
    
    # return matrix as follows: nbins x nlambdas)
    return(ws)
  }

ihw_convex_weights <-
  function(split_sorted_pvalues,
           alpha,
           m_groups,
           m_groups_grenander,
           penalty = "total variation",
           lambdas,
           admm = T,
           adjustment_type = "BH",
           grenander_binsize = 1,
           quiet = quiet) {
    # m_groups used for balancing (weight budget)
    # m_groups_grenander (used for grenander estimator)
    # Practically always it will hold: m_groups \approx m_groups_grenander
    
    nbins <- length(split_sorted_pvalues)
    
    if (length(lambdas) == 1 && lambdas[1] == Inf) {
      ws <- matrix(1, nbins, 1)
      return(ws)
    }
    
    zero_included <- any(lambdas == 0)
    infty_included <- any(is.infinite(lambdas))
    
    lambdas_filt <- lambdas[!(lambdas %in% c(0, Inf))]
    # also special case if lambdas == 0
    
    # also there temporarily to help with numerical stability of downstream optimization
    max_pval <- 0.9 * max(sapply(split_sorted_pvalues, max))
    if (max_pval < 0.1) {
      split_sorted_pvalues <-
        lapply(split_sorted_pvalues, function(x)
          x / max_pval)
      alpha <- alpha / max_pval
    }
    # preprocessing:  Set very low p-values to 10^(-20) to help with numerical stability
    # Note: This only affects the internal optimization, the higher level functions will
    # still return adjusted pvalues based on the original p-values
    split_sorted_pvalues <-
      lapply(split_sorted_pvalues, function(x)
        ifelse(x > 10 ^ (-20), x, 10 ^ (-20)))
    m <- sum(m_groups)
    
    if (nbins != length(m_groups)) {
      stop("length of m_groups should be equal to number of bins")
    }
    
    #lapply grenander...
    if (!quiet)
      message("Applying Grenander estimator within each bin.")
    
    grenander_list <- mapply(
      presorted_grenander,
      split_sorted_pvalues,
      m_groups_grenander,
      grenander_binsize = grenander_binsize,
      quiet = quiet,
      SIMPLIFY = FALSE
    )
    
    if (length(lambdas_filt) == 0 && zero_included) {
      ts <-  unregularized_thresholds_bh(grenander_list, m_groups, alpha)
      if (infty_included) {
        ts <- cbind(1, ts)
      } else {
        ts <- matrix(ts, ncol = 1)
      }
    } else {
      if (admm){
        ts <-  optimal_ts(grenander_list, m_groups, alpha, lambdas_filt)
      }else{
        #add function cvxr here
        ts <- optimal_cvxr(grenander_list, m_groups, alpha, lambdas_filt)
      }
      
      if (!zero_included) {
        ts <- ts[,-ncol(ts), drop = FALSE]
      }
      if (!infty_included) {
        ts <- ts[,-1, drop = FALSE]
      }
      ts <- pmax(ts, 0)
    }
    
    ws <- apply(ts, 2, function(x)
      thresholds_to_weights(x, m_groups))
    
    # return matrix as follows: nbins x nlambdas)
    return(ws)
  }

presorted_grenander <-
  function(sorted_pvalues,
           m_total = length(sorted_pvalues),
           grenander_binsize = 1,
           quiet = TRUE) {
    unique_pvalues <- unique(sorted_pvalues)
    ecdf_values <-
      cumsum(tabulate(match(sorted_pvalues, unique_pvalues))) / m_total
    
    # I think fdrtool neglects this borderline case and this causes returned object
    # to be discontinuous hence also not concave
    unique_pvalues <- c(0, unique_pvalues)
    ecdf_values   <- c(0, ecdf_values)
    
    
    if (grenander_binsize != 1) {
      nmax = length(unique_pvalues)
      idx_thin <- c(seq(1, nmax - 1, by = grenander_binsize), nmax)
      unique_pvalues <- unique_pvalues[idx_thin]
      ecdf_values <- ecdf_values[idx_thin]
    }
    
    ll <- fdrtool::gcmlcm(unique_pvalues, ecdf_values, type = "lcm")
    
    gren <-
      list(
        locs = ll$x.knots,
        Fs = ll$y.knots,
        fs = ll$slope.knots
      )
    if (!quiet)
      message(paste("Grenander fit with", length(ll$slope.knots), "knots."))
    gren
  }


