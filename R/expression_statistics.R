# Compute per-feature statistics between two survival-defined expression matrices.
# Input: expr_high, expr_low
#   - standardized expression objects: rows = sample_id, cols = feature_coord
# Output: data.frame with the required "intermediate/significants/*" header schema.

compute_significants_from_split <- function(
    expr_low,
    expr_high,
    clinical_high = NULL,
    clinical_low  = NULL,
    surv_time_col   = NULL,
    surv_status_col = NULL
) {
  # ---- helper: coerce to numeric matrix ----
  to_num_matrix <- function(x, x_name) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (!is.matrix(x)) stop(sprintf("%s must be a matrix or data.frame", x_name), call. = FALSE)
    storage.mode(x) <- "numeric"
    x
  }
  
  # ---- helper: safe p-value wrappers ----
  safe_t_p <- function(a, b) {
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 2L || length(b) < 2L) return(NA_real_)
    if (sd(a) == 0 && sd(b) == 0) return(NA_real_)
    out <- try(stats::t.test(a, b), silent = TRUE)
    if (inherits(out, "try-error")) NA_real_ else as.numeric(out$p.value)
  }
  
  safe_ks <- function(a, b) {
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 1L || length(b) < 1L) return(c(NA_real_, NA_real_))
    out <- try(stats::ks.test(a, b), silent = TRUE)
    if (inherits(out, "try-error")) return(c(NA_real_, NA_real_))
    c(as.numeric(out$p.value), as.numeric(out$statistic))
  }
  
  safe_wilcox_p <- function(a, b) {
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 1L || length(b) < 1L) return(NA_real_)
    out <- try(stats::wilcox.test(a, b, exact = FALSE), silent = TRUE)
    if (inherits(out, "try-error")) NA_real_ else as.numeric(out$p.value)
  }
  
  safe_logrank_p <- function(time_high, status_high, time_low, status_low) {
    if (is.null(time_high) || is.null(status_high) || is.null(time_low) || is.null(status_low)) return(NA_real_)
    if (!requireNamespace("survival", quietly = TRUE)) return(NA_real_)
    
    time_high   <- suppressWarnings(as.numeric(time_high))
    time_low    <- suppressWarnings(as.numeric(time_low))
    status_high_num <- suppressWarnings(as.numeric(status_high))
    status_low_num  <- suppressWarnings(as.numeric(status_low))
    
    if (all(!is.finite(status_high_num)) && any(!is.na(status_high))) {
      status_high_num <- vital_status_to_event(status_high)
    }
    if (all(!is.finite(status_low_num)) && any(!is.na(status_low))) {
      status_low_num <- vital_status_to_event(status_low)
    }
    
    keep_h <- is.finite(time_high) & is.finite(status_high_num)
    keep_l <- is.finite(time_low)  & is.finite(status_low_num)
    
    time_high   <- time_high[keep_h];   status_high_num <- status_high_num[keep_h]
    time_low    <- time_low[keep_l];    status_low_num  <- status_low_num[keep_l]
    
    if (length(time_high) < 2L || length(time_low) < 2L) return(NA_real_)
    
    time   <- c(time_high, time_low)
    status <- c(status_high_num, status_low_num)
    group  <- c(rep("high", length(time_high)), rep("low", length(time_low)))
    
    fit <- try(survival::survdiff(survival::Surv(time, status) ~ group), silent = TRUE)
    if (inherits(fit, "try-error")) return(NA_real_)
    
    chisq <- as.numeric(fit$chisq)
    stats::pchisq(chisq, df = 1, lower.tail = FALSE)
  }
  
  # ---- coerce ----
  high_mat <- to_num_matrix(expr_high, "expr_high")
  low_mat  <- to_num_matrix(expr_low,  "expr_low")
  
  # ---- align features (columns) ----
  common_features <- intersect(colnames(high_mat), colnames(low_mat))
  if (length(common_features) == 0L) stop("No overlapping feature_coord columns between expr_high and expr_low.", call. = FALSE)
  # the code stoped in the line before this comment
  
  high_mat <- high_mat[, common_features, drop = FALSE]
  low_mat  <- low_mat[,  common_features, drop = FALSE]
  
  # ---- compute mean/diff ----
  mean_high <- colMeans(high_mat, na.rm = TRUE)
  mean_low  <- colMeans(low_mat,  na.rm = TRUE)
  diff_mean <- mean_high - mean_low
  
  # ---- stats per feature ----
  t_p   <- rep(NA_real_, length(common_features))
  ks_p  <- rep(NA_real_, length(common_features))
  ks_d  <- rep(NA_real_, length(common_features))
  w_p   <- rep(NA_real_, length(common_features))
  
  for (i in seq_along(common_features)) {
    a <- high_mat[, i]
    b <- low_mat[,  i]
    
    t_p[i] <- safe_t_p(a, b)
    
    ks_out <- safe_ks(a, b)
    ks_p[i] <- ks_out[1]
    ks_d[i] <- ks_out[2]
    
    w_p[i] <- safe_wilcox_p(a, b)
  }
  
  # ---- log-rank p-value (group-level, not per-feature) ----
  # If clinical_high/clinical_low are provided with time/status columns, compute one LogRank_P and repeat it.
  logrank_p <- NA_real_
  if (!is.null(clinical_high) && !is.null(clinical_low) &&
      !is.null(surv_time_col) && !is.null(surv_status_col)) {
    
    if (!is.data.frame(clinical_high) || !is.data.frame(clinical_low)) {
      stop("clinical_high and clinical_low must be data.frames when provided.", call. = FALSE)
    }
    if (!(surv_time_col %in% names(clinical_high)) || !(surv_status_col %in% names(clinical_high))) {
      stop("clinical_high is missing surv_time_col and/or surv_status_col.", call. = FALSE)
    }
    if (!(surv_time_col %in% names(clinical_low)) || !(surv_status_col %in% names(clinical_low))) {
      stop("clinical_low is missing surv_time_col and/or surv_status_col.", call. = FALSE)
    }
    
    logrank_p <- safe_logrank_p(
      time_high   = clinical_high[[surv_time_col]],
      status_high = clinical_high[[surv_status_col]],
      time_low    = clinical_low[[surv_time_col]],
      status_low  = clinical_low[[surv_status_col]]
    )
  }
  
  # ---- build results (required header) ----
  out <- data.frame(
    feature_coord = as.character(common_features),
    mean_high     = as.numeric(mean_high),
    mean_low      = as.numeric(mean_low),
    diff_mean     = as.numeric(diff_mean),
    Ttest_P       = as.numeric(t_p),
    KS_P          = as.numeric(ks_p),
    KS_D          = as.numeric(ks_d),
    Wilcoxon_P    = as.numeric(w_p),
    LogRank_P     = rep(as.numeric(logrank_p), length(common_features)),
    stringsAsFactors = FALSE
  )
  
  # optional sanity check: required columns for "intermediate/significants/*"
  req <- c("feature_coord","mean_high","mean_low","diff_mean","Ttest_P","KS_P","KS_D","Wilcoxon_P","LogRank_P")
  miss <- setdiff(req, names(out))
  if (length(miss) > 0L) stop(paste("Internal error: missing columns:", paste(miss, collapse = ", ")), call. = FALSE)
  
  out
}

# Filter significants according to analysis design:
# keep feature_coord if at least 2/3 of (Ttest, KS, Wilcoxon)
# are significant (p <= max_p_value).
filter_significants_any_test <- function(out, max_p_value) {
  
  if (!is.data.frame(out)) {
    stop("out must be a data.frame", call. = FALSE)
  }
  
  if (!is.numeric(max_p_value) || length(max_p_value) != 1 || !is.finite(max_p_value)) {
    stop("max_p_value must be a single finite numeric value", call. = FALSE)
  }
  
  required_cols <- c(
    "feature_coord",
    "Ttest_P",
    "KS_P",
    "Wilcoxon_P"
  )
  
  missing_cols <- setdiff(required_cols, colnames(out))
  if (length(missing_cols) > 0L) {
    stop(
      paste("out is missing required columns:",
            paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  
  pmat_3 <- as.matrix(out[, c("Ttest_P", "KS_P", "Wilcoxon_P")])
  storage.mode(pmat_3) <- "numeric"
  
  sig_3 <- is.finite(pmat_3) & !is.na(pmat_3) & (pmat_3 <= max_p_value)
  
  keep <- (rowSums(sig_3) >= 2L)
  
  out[keep, , drop = FALSE]
}

# ---- specification-driven functions ----

vital_status_to_event <- function(vital_status_vec) {
  v <- tolower(as.character(vital_status_vec))
  out <- ifelse(v == "dead", 1L, 0L)
  as.integer(out)
}

coerce_expression_numeric_vector <- function(expr_vec, context = "expression vector") {
  if (is.data.frame(expr_vec)) {
    if (ncol(expr_vec) < 1L) {
      stop(sprintf("%s is an empty data.frame.", context), call. = FALSE)
    }
    if (ncol(expr_vec) == 1L) {
      expr_vec <- expr_vec[[1]]
    } else {
      mat <- as.matrix(expr_vec)
      storage.mode(mat) <- "numeric"
      return(rowMeans(mat, na.rm = TRUE))
    }
  }
  if (is.list(expr_vec) && !is.data.frame(expr_vec)) {
    lens <- lengths(expr_vec)
    if (all(lens == 1L)) {
      expr_vec <- unlist(expr_vec, use.names = FALSE)
    } else {
      stop(
        sprintf("%s is a list with non-scalar elements; unnest/flatten before running.", context),
        call. = FALSE
      )
    }
  }
  suppressWarnings(as.numeric(expr_vec))
}

coerce_expression_df_numeric <- function(expr_df, patient_id_col = "patient_id") {
  if (!is.data.frame(expr_df)) stop("Expression input must be a data.frame.", call. = FALSE)
  if (!(patient_id_col %in% names(expr_df))) stop("patient_id_col missing.", call. = FALSE)
  feat_cols <- setdiff(names(expr_df), patient_id_col)
  if (length(feat_cols) < 1L) return(expr_df)
  for (feat in feat_cols) {
    expr_df[[feat]] <- coerce_expression_numeric_vector(expr_df[[feat]], context = paste0("feature '", feat, "'"))
  }
  expr_df
}

build_expression_split_inputs <- function(
    high_survival_expr,
    low_survival_expr,
    high_survival_clin,
    low_survival_clin,
    patient_id_col = "patient_id",
    time_col = "days_to_last_follow_up",
    status_col = "vital_status"
) {
  as_df <- function(x) {
    if (is.matrix(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
    if (!is.data.frame(x)) stop("Expression inputs must be data.frames or matrices.", call. = FALSE)
    x
  }
  
  high_survival_expr <- as_df(high_survival_expr)
  low_survival_expr  <- as_df(low_survival_expr)
  
  if (!is.data.frame(high_survival_clin) || !is.data.frame(low_survival_clin)) {
    stop("Clinical inputs must be data.frames.", call. = FALSE)
  }
  
  if (!(patient_id_col %in% names(high_survival_expr)) ||
      !(patient_id_col %in% names(low_survival_expr))) {
    stop("patient_id_col missing from expression inputs.", call. = FALSE)
  }
  if (!(patient_id_col %in% names(high_survival_clin)) ||
      !(patient_id_col %in% names(low_survival_clin))) {
    stop("patient_id_col missing from clinical inputs.", call. = FALSE)
  }
  if (!(time_col %in% names(high_survival_clin)) ||
      !(status_col %in% names(high_survival_clin)) ||
      !(time_col %in% names(low_survival_clin)) ||
      !(status_col %in% names(low_survival_clin))) {
    stop("Clinical inputs must include time and status columns.", call. = FALSE)
  }
  
  feature_cols_high <- setdiff(names(high_survival_expr), patient_id_col)
  feature_cols_low  <- setdiff(names(low_survival_expr),  patient_id_col)
  if (!identical(feature_cols_high, feature_cols_low)) {
    stop("High/low expression inputs must have identical feature columns.", call. = FALSE)
  }

  high_survival_expr <- coerce_expression_df_numeric(high_survival_expr, patient_id_col = patient_id_col)
  low_survival_expr  <- coerce_expression_df_numeric(low_survival_expr,  patient_id_col = patient_id_col)
  
  # align expr/clin within each group by patient_id
  common_high <- intersect(high_survival_expr[[patient_id_col]], high_survival_clin[[patient_id_col]])
  common_low  <- intersect(low_survival_expr[[patient_id_col]],  low_survival_clin[[patient_id_col]])
  
  high_survival_expr <- high_survival_expr[high_survival_expr[[patient_id_col]] %in% common_high, , drop = FALSE]
  high_survival_clin <- high_survival_clin[high_survival_clin[[patient_id_col]] %in% common_high, , drop = FALSE]
  low_survival_expr  <- low_survival_expr[low_survival_expr[[patient_id_col]] %in% common_low, , drop = FALSE]
  low_survival_clin  <- low_survival_clin[low_survival_clin[[patient_id_col]] %in% common_low, , drop = FALSE]
  
  # order by time, then vital_status (dead ranks less)
  order_by_time_status <- function(clin_df) {
    time_vals <- suppressWarnings(as.numeric(clin_df[[time_col]]))
    status_rank <- ifelse(tolower(as.character(clin_df[[status_col]])) == "dead", 0, 1)
    order(time_vals, status_rank)
  }
  
  ord_h <- order_by_time_status(high_survival_clin)
  high_survival_clin <- high_survival_clin[ord_h, , drop = FALSE]
  high_survival_expr <- high_survival_expr[match(high_survival_clin[[patient_id_col]], high_survival_expr[[patient_id_col]]), , drop = FALSE]
  
  ord_l <- order_by_time_status(low_survival_clin)
  low_survival_clin <- low_survival_clin[ord_l, , drop = FALSE]
  low_survival_expr <- low_survival_expr[match(low_survival_clin[[patient_id_col]], low_survival_expr[[patient_id_col]]), , drop = FALSE]
  
  expr_df <- rbind(high_survival_expr, low_survival_expr)
  clin_df <- rbind(high_survival_clin, low_survival_clin)
  group_vec <- c(rep("HighSurvival", nrow(high_survival_expr)),
                 rep("LowSurvival",  nrow(low_survival_expr)))
  
  # final ordering based on combined clinical data
  ord_all <- order_by_time_status(clin_df)
  expr_df <- expr_df[ord_all, , drop = FALSE]
  clin_df <- clin_df[ord_all, , drop = FALSE]
  group_vec <- group_vec[ord_all]
  
  list(expr_df = expr_df, clin_df = clin_df, group_vec = group_vec)
}

compute_means_per_feature <- function(
    high_survival_expr,
    low_survival_expr,
    patient_id_col = "patient_id"
) {
  to_num_matrix <- function(df) {
    if (is.matrix(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (!is.data.frame(df)) stop("Expression input must be data.frame or matrix.", call. = FALSE)
    if (!(patient_id_col %in% names(df))) stop("patient_id_col missing.", call. = FALSE)
    feats <- setdiff(names(df), patient_id_col)
    m <- as.matrix(df[, feats, drop = FALSE])
    storage.mode(m) <- "numeric"
    m
  }
  
  m_high <- to_num_matrix(high_survival_expr)
  m_low  <- to_num_matrix(low_survival_expr)
  if (!identical(colnames(m_high), colnames(m_low))) {
    stop("High/low expression inputs must have identical feature columns.", call. = FALSE)
  }
  
  mean_high <- colMeans(m_high, na.rm = TRUE)
  mean_low  <- colMeans(m_low,  na.rm = TRUE)
  
  list(mean_high = mean_high, mean_low = mean_low)
}

compute_diff_mean_per_feature <- function(mean_high, mean_low) {
  mean_high <- as.numeric(mean_high)
  mean_low  <- as.numeric(mean_low)
  names(mean_high) <- names(mean_low)
  diff_mean <- mean_high - mean_low
  diff_mean
}

compute_fold_change_per_feature <- function(mean_high, mean_low, pseudo_count = 1) {
  fold_change <- (mean_high + pseudo_count) / (mean_low + pseudo_count)
  fold_change
}

compute_log2_fold_change_per_feature <- function(fold_change) {
  log2(fold_change)
}

compute_ttest_per_feature <- function(
    high_survival_expr,
    low_survival_expr,
    patient_id_col = "patient_id"
) {
  to_num_matrix <- function(df) {
    if (is.matrix(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (!is.data.frame(df)) stop("Expression input must be data.frame or matrix.", call. = FALSE)
    if (!(patient_id_col %in% names(df))) stop("patient_id_col missing.", call. = FALSE)
    feats <- setdiff(names(df), patient_id_col)
    m <- as.matrix(df[, feats, drop = FALSE])
    storage.mode(m) <- "numeric"
    m
  }
  
  m_high <- to_num_matrix(high_survival_expr)
  m_low  <- to_num_matrix(low_survival_expr)
  if (!identical(colnames(m_high), colnames(m_low))) {
    stop("High/low expression inputs must have identical feature columns.", call. = FALSE)
  }
  
  t_p <- rep(NA_real_, ncol(m_high))
  names(t_p) <- colnames(m_high)
  warn_count <- 0L
  
  for (i in seq_len(ncol(m_high))) {
    a <- m_high[, i]
    b <- m_low[,  i]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 2L || length(b) < 2L || (sd(a) == 0 && sd(b) == 0)) {
      warn_count <- warn_count + 1L
      next
    }
    out <- try(stats::t.test(a, b), silent = TRUE)
    if (inherits(out, "try-error")) {
      warn_count <- warn_count + 1L
    } else {
      t_p[i] <- as.numeric(out$p.value)
    }
  }
  
  if (warn_count > 0L) {
    warning(sprintf("T-test undefined for %d feature(s).", warn_count), call. = FALSE)
  }
  
  t_p
}

compute_ks_per_feature <- function(
    high_survival_expr,
    low_survival_expr,
    patient_id_col = "patient_id"
) {
  to_num_matrix <- function(df) {
    if (is.matrix(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (!is.data.frame(df)) stop("Expression input must be data.frame or matrix.", call. = FALSE)
    if (!(patient_id_col %in% names(df))) stop("patient_id_col missing.", call. = FALSE)
    feats <- setdiff(names(df), patient_id_col)
    m <- as.matrix(df[, feats, drop = FALSE])
    storage.mode(m) <- "numeric"
    m
  }
  
  m_high <- to_num_matrix(high_survival_expr)
  m_low  <- to_num_matrix(low_survival_expr)
  if (!identical(colnames(m_high), colnames(m_low))) {
    stop("High/low expression inputs must have identical feature columns.", call. = FALSE)
  }
  
  ks_p <- rep(NA_real_, ncol(m_high))
  ks_d <- rep(NA_real_, ncol(m_high))
  names(ks_p) <- colnames(m_high)
  names(ks_d) <- colnames(m_high)
  warn_count <- 0L
  
  for (i in seq_len(ncol(m_high))) {
    a <- m_high[, i]
    b <- m_low[,  i]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 1L || length(b) < 1L) {
      warn_count <- warn_count + 1L
      next
    }
    out <- try(stats::ks.test(a, b), silent = TRUE)
    if (inherits(out, "try-error")) {
      warn_count <- warn_count + 1L
    } else {
      ks_p[i] <- as.numeric(out$p.value)
      ks_d[i] <- as.numeric(out$statistic)
    }
  }
  
  if (warn_count > 0L) {
    warning(sprintf("KS test undefined for %d feature(s).", warn_count), call. = FALSE)
  }
  
  list(KS_D = ks_d, KS_P = ks_p)
}

compute_wilcoxon_per_feature <- function(
    high_survival_expr,
    low_survival_expr,
    patient_id_col = "patient_id"
) {
  to_num_matrix <- function(df) {
    if (is.matrix(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (!is.data.frame(df)) stop("Expression input must be data.frame or matrix.", call. = FALSE)
    if (!(patient_id_col %in% names(df))) stop("patient_id_col missing.", call. = FALSE)
    feats <- setdiff(names(df), patient_id_col)
    m <- as.matrix(df[, feats, drop = FALSE])
    storage.mode(m) <- "numeric"
    m
  }
  
  m_high <- to_num_matrix(high_survival_expr)
  m_low  <- to_num_matrix(low_survival_expr)
  if (!identical(colnames(m_high), colnames(m_low))) {
    stop("High/low expression inputs must have identical feature columns.", call. = FALSE)
  }
  
  w_p <- rep(NA_real_, ncol(m_high))
  names(w_p) <- colnames(m_high)
  warn_count <- 0L
  
  for (i in seq_len(ncol(m_high))) {
    a <- m_high[, i]
    b <- m_low[,  i]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 1L || length(b) < 1L) {
      warn_count <- warn_count + 1L
      next
    }
    out <- try(stats::wilcox.test(a, b, exact = FALSE), silent = TRUE)
    if (inherits(out, "try-error")) {
      warn_count <- warn_count + 1L
    } else {
      w_p[i] <- as.numeric(out$p.value)
    }
  }
  
  if (warn_count > 0L) {
    warning(sprintf("Wilcoxon test undefined for %d feature(s).", warn_count), call. = FALSE)
  }
  
  w_p
}

cox_fit_one_feature_stratified <- function(
    expr_vec,
    clin_df,
    group_vec,
    time_col = "days_to_last_follow_up",
    status_col = "vital_status"
) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    warning("Package 'survival' is not available.", call. = FALSE)
    return(list(Cox_beta = NA_real_, Cox_HR = NA_real_, Cox_P = NA_real_, Cox_n = NA_real_))
  }
  
  if (!is.data.frame(clin_df)) {
    stop("clin_df must be a data.frame.", call. = FALSE)
  }
  if (!(time_col %in% names(clin_df)) || !(status_col %in% names(clin_df))) {
    stop("clin_df missing time/status columns.", call. = FALSE)
  }
  
  time_vals <- suppressWarnings(as.numeric(clin_df[[time_col]]))
  status_vals <- vital_status_to_event(clin_df[[status_col]])
  expr_vec <- coerce_expression_numeric_vector(expr_vec, context = "ROC expression vector")
  group_vec <- as.character(group_vec)
  
  keep <- is.finite(time_vals) & is.finite(status_vals) & is.finite(expr_vec) & !is.na(group_vec)
  time_vals <- time_vals[keep]
  status_vals <- status_vals[keep]
  expr_vec <- expr_vec[keep]
  group_vec <- group_vec[keep]
  
  if (length(time_vals) < 2L) {
    warning("Not enough observations for Cox model.", call. = FALSE)
    return(list(Cox_beta = NA_real_, Cox_HR = NA_real_, Cox_P = NA_real_, Cox_n = length(time_vals)))
  }
  
  df <- data.frame(
    time = time_vals,
    status = status_vals,
    expr = expr_vec,
    group = group_vec,
    stringsAsFactors = FALSE
  )
  
  fit <- try(survival::coxph(survival::Surv(time, status) ~ expr + strata(group), data = df), silent = TRUE)
  if (inherits(fit, "try-error")) {
    warning("Cox model failed for feature.", call. = FALSE)
    return(list(Cox_beta = NA_real_, Cox_HR = NA_real_, Cox_P = NA_real_, Cox_n = nrow(df)))
  }
  
  coef_val <- as.numeric(stats::coef(fit)[1])
  hr_val <- exp(coef_val)
  p_val <- as.numeric(summary(fit)$coef[1, "Pr(>|z|)"])
  
  list(Cox_beta = coef_val, Cox_HR = hr_val, Cox_P = p_val, Cox_n = nrow(df))
}

cox_fit_all_features_stratified <- function(
    expr_df,
    clin_df,
    group_vec,
    patient_id_col = "patient_id",
    time_col = "days_to_last_follow_up",
    status_col = "vital_status"
) {
  if (!is.data.frame(expr_df)) stop("expr_df must be a data.frame.", call. = FALSE)
  if (!(patient_id_col %in% names(expr_df))) stop("patient_id_col missing from expr_df.", call. = FALSE)
  
  feature_cols <- setdiff(names(expr_df), patient_id_col)
  out <- data.frame(
    feature_coord = feature_cols,
    Cox_beta = NA_real_,
    Cox_HR = NA_real_,
    Cox_P = NA_real_,
    Cox_n = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(feature_cols)) {
    feat <- feature_cols[i]
    res <- cox_fit_one_feature_stratified(
      expr_vec = expr_df[[feat]],
      clin_df = clin_df,
      group_vec = group_vec,
      time_col = time_col,
      status_col = status_col
    )
    out$Cox_beta[i] <- res$Cox_beta
    out$Cox_HR[i]   <- res$Cox_HR
    out$Cox_P[i]    <- res$Cox_P
    out$Cox_n[i]    <- res$Cox_n
  }
  
  out
}

cox_fit_one_feature_unstratified <- function(
    expr_vec,
    clin_df,
    time_col = "days_to_last_follow_up",
    status_col = "vital_status"
) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    warning("Package 'survival' is not available.", call. = FALSE)
    return(list(Cox_beta = NA_real_, Cox_HR = NA_real_, Cox_P = NA_real_, Cox_n = NA_real_))
  }
  
  if (!is.data.frame(clin_df)) {
    stop("clin_df must be a data.frame.", call. = FALSE)
  }
  if (!(time_col %in% names(clin_df)) || !(status_col %in% names(clin_df))) {
    stop("clin_df missing time/status columns.", call. = FALSE)
  }
  
  time_vals <- suppressWarnings(as.numeric(clin_df[[time_col]]))
  status_vals <- vital_status_to_event(clin_df[[status_col]])
  expr_vec <- coerce_expression_numeric_vector(expr_vec, context = "Cox expression vector")
  
  keep <- is.finite(time_vals) & is.finite(status_vals) & is.finite(expr_vec)
  time_vals <- time_vals[keep]
  status_vals <- status_vals[keep]
  expr_vec <- expr_vec[keep]
  
  if (length(time_vals) < 2L) {
    warning("Not enough observations for Cox model.", call. = FALSE)
    return(list(Cox_beta = NA_real_, Cox_HR = NA_real_, Cox_P = NA_real_, Cox_n = length(time_vals)))
  }
  
  df <- data.frame(
    time = time_vals,
    status = status_vals,
    expr = expr_vec,
    stringsAsFactors = FALSE
  )
  
  fit <- try(survival::coxph(survival::Surv(time, status) ~ expr, data = df), silent = TRUE)
  if (inherits(fit, "try-error")) {
    warning("Cox model failed for feature.", call. = FALSE)
    return(list(Cox_beta = NA_real_, Cox_HR = NA_real_, Cox_P = NA_real_, Cox_n = nrow(df)))
  }
  
  coef_val <- as.numeric(stats::coef(fit)[1])
  hr_val <- exp(coef_val)
  p_val <- as.numeric(summary(fit)$coef[1, "Pr(>|z|)"])
  
  list(Cox_beta = coef_val, Cox_HR = hr_val, Cox_P = p_val, Cox_n = nrow(df))
}

cox_fit_all_features_unstratified <- function(
    expr_df,
    clin_df,
    patient_id_col = "patient_id",
    time_col = "days_to_last_follow_up",
    status_col = "vital_status"
) {
  if (!is.data.frame(expr_df)) stop("expr_df must be a data.frame.", call. = FALSE)
  if (!(patient_id_col %in% names(expr_df))) stop("patient_id_col missing from expr_df.", call. = FALSE)
  
  feature_cols <- setdiff(names(expr_df), patient_id_col)
  out <- data.frame(
    feature_coord = feature_cols,
    Cox_beta = NA_real_,
    Cox_HR = NA_real_,
    Cox_P = NA_real_,
    Cox_n = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(feature_cols)) {
    feat <- feature_cols[i]
    res <- cox_fit_one_feature_unstratified(
      expr_vec = expr_df[[feat]],
      clin_df = clin_df,
      time_col = time_col,
      status_col = status_col
    )
    out$Cox_beta[i] <- res$Cox_beta
    out$Cox_HR[i]   <- res$Cox_HR
    out$Cox_P[i]    <- res$Cox_P
    out$Cox_n[i]    <- res$Cox_n
  }
  
  out
}

roc_auc_and_cutoff_youden_per_feature <- function(expr_vec, group_vec) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    warning("Package 'pROC' is not available.", call. = FALSE)
    return(list(ROC_AUC = NA_real_, cutoff_value = NA_real_, youden_J = NA_real_))
  }
  
  expr_vec <- coerce_expression_numeric_vector(expr_vec, context = "ROC expression vector")
  group_vec <- as.character(group_vec)
  response <- factor(group_vec, levels = c("LowSurvival", "HighSurvival"))
  
  keep <- is.finite(expr_vec) & !is.na(response)
  expr_vec <- expr_vec[keep]
  response <- response[keep]
  
  if (length(unique(response)) < 2L) {
    warning("ROC requires two groups.", call. = FALSE)
    return(list(ROC_AUC = NA_real_, cutoff_value = NA_real_, youden_J = NA_real_))
  }
  
  roc_obj <- try(pROC::roc(response = response, predictor = expr_vec, quiet = TRUE), silent = TRUE)
  if (inherits(roc_obj, "try-error")) {
    warning("ROC failed for feature.", call. = FALSE)
    return(list(ROC_AUC = NA_real_, cutoff_value = NA_real_, youden_J = NA_real_))
  }
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  coords <- pROC::coords(roc_obj, "best", best.method = "youden",
                         ret = c("threshold", "sensitivity", "specificity"))

  extract_coord <- function(coords, name) {
    if (is.null(coords)) return(NA_real_)
    if (is.list(coords) && !is.data.frame(coords)) {
      if (!is.null(coords[[name]])) return(as.numeric(coords[[name]])[1])
    }
    if (is.matrix(coords) || is.data.frame(coords)) {
      if (!is.null(colnames(coords)) && name %in% colnames(coords)) {
        return(as.numeric(coords[1, name]))
      }
      if (!is.null(rownames(coords)) && name %in% rownames(coords)) {
        return(as.numeric(coords[name, 1]))
      }
    }
    as.numeric(coords[name])
  }
  
  threshold_val <- extract_coord(coords, "threshold")
  sens_val <- extract_coord(coords, "sensitivity")
  spec_val <- extract_coord(coords, "specificity")
  youden_val <- if (is.finite(sens_val) && is.finite(spec_val)) {
    sens_val + spec_val - 1
  } else {
    coords2 <- pROC::coords(roc_obj, "best", best.method = "youden", ret = c("youden"))
    extract_coord(coords2, "youden")
  }
  
  list(
    ROC_AUC = auc_val,
    cutoff_value = threshold_val,
    youden_J = youden_val
  )
}

roc_auc_and_cutoff_youden_all_features <- function(
    expr_df,
    group_vec,
    patient_id_col = "patient_id"
) {
  if (!is.data.frame(expr_df)) stop("expr_df must be a data.frame.", call. = FALSE)
  if (!(patient_id_col %in% names(expr_df))) stop("patient_id_col missing from expr_df.", call. = FALSE)
  
  feature_cols <- setdiff(names(expr_df), patient_id_col)
  out <- data.frame(
    feature_coord = feature_cols,
    ROC_AUC = NA_real_,
    cutoff_value = NA_real_,
    youden_J = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(feature_cols)) {
    feat <- feature_cols[i]
    res <- roc_auc_and_cutoff_youden_per_feature(expr_df[[feat]], group_vec)
    out$ROC_AUC[i] <- res$ROC_AUC
    out$cutoff_value[i] <- res$cutoff_value
    out$youden_J[i] <- res$youden_J
  }
  
  out
}

roc_auc_and_cutoff_youden_per_feature_event <- function(expr_vec, status_vec) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    warning("Package 'pROC' is not available.", call. = FALSE)
    return(list(ROC_AUC = NA_real_, cutoff_value = NA_real_, youden_J = NA_real_, ROC_direction = NA_character_))
  }
  
  expr_vec <- coerce_expression_numeric_vector(expr_vec, context = "ROC expression vector")
  status_vec <- vital_status_to_event(status_vec)
  
  keep <- is.finite(expr_vec) & is.finite(status_vec)
  expr_vec <- expr_vec[keep]
  status_vec <- status_vec[keep]
  
  if (length(unique(status_vec)) < 2L) {
    warning("ROC requires two outcome classes.", call. = FALSE)
    return(list(ROC_AUC = NA_real_, cutoff_value = NA_real_, youden_J = NA_real_, ROC_direction = NA_character_))
  }
  
  response <- factor(status_vec, levels = c(0, 1))
  roc_obj <- try(pROC::roc(response = response, predictor = expr_vec, quiet = TRUE, direction = "auto"), silent = TRUE)
  if (inherits(roc_obj, "try-error")) {
    warning("ROC failed for feature.", call. = FALSE)
    return(list(ROC_AUC = NA_real_, cutoff_value = NA_real_, youden_J = NA_real_, ROC_direction = NA_character_))
  }
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  coords <- pROC::coords(roc_obj, "best", best.method = "youden",
                         ret = c("threshold", "sensitivity", "specificity"))
  
  extract_coord <- function(coords, name) {
    if (is.null(coords)) return(NA_real_)
    if (is.list(coords) && !is.data.frame(coords)) {
      if (!is.null(coords[[name]])) return(as.numeric(coords[[name]])[1])
    }
    if (is.matrix(coords) || is.data.frame(coords)) {
      if (!is.null(colnames(coords)) && name %in% colnames(coords)) {
        return(as.numeric(coords[1, name]))
      }
      if (!is.null(rownames(coords)) && name %in% rownames(coords)) {
        return(as.numeric(coords[name, 1]))
      }
    }
    as.numeric(coords[name])
  }
  
  threshold_val <- extract_coord(coords, "threshold")
  sens_val <- extract_coord(coords, "sensitivity")
  spec_val <- extract_coord(coords, "specificity")
  youden_val <- if (is.finite(sens_val) && is.finite(spec_val)) {
    sens_val + spec_val - 1
  } else {
    coords2 <- pROC::coords(roc_obj, "best", best.method = "youden", ret = c("youden"))
    extract_coord(coords2, "youden")
  }
  
  list(
    ROC_AUC = auc_val,
    cutoff_value = threshold_val,
    youden_J = youden_val,
    ROC_direction = roc_obj$direction
  )
}

roc_auc_and_cutoff_youden_all_features_event <- function(
    expr_df,
    status_vec,
    patient_id_col = "patient_id"
) {
  if (!is.data.frame(expr_df)) stop("expr_df must be a data.frame.", call. = FALSE)
  if (!(patient_id_col %in% names(expr_df))) stop("patient_id_col missing from expr_df.", call. = FALSE)
  
  feature_cols <- setdiff(names(expr_df), patient_id_col)
  out <- data.frame(
    feature_coord = feature_cols,
    ROC_AUC = NA_real_,
    cutoff_value = NA_real_,
    youden_J = NA_real_,
    ROC_direction = NA_character_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(feature_cols)) {
    feat <- feature_cols[i]
    res <- roc_auc_and_cutoff_youden_per_feature_event(expr_df[[feat]], status_vec)
    out$ROC_AUC[i] <- res$ROC_AUC
    out$cutoff_value[i] <- res$cutoff_value
    out$youden_J[i] <- res$youden_J
    out$ROC_direction[i] <- res$ROC_direction
  }
  
  out
}

compute_error_rate_from_cutoff <- function(expr_vec, status_vec, cutoff_value, direction) {
  expr_vec <- coerce_expression_numeric_vector(expr_vec, context = "expression vector for error rate")
  status_vec <- vital_status_to_event(status_vec)
  
  if (!is.finite(cutoff_value)) return(NA_real_)
  
  keep <- is.finite(expr_vec) & is.finite(status_vec)
  expr_vec <- expr_vec[keep]
  status_vec <- status_vec[keep]
  
  if (length(expr_vec) == 0L) return(NA_real_)
  
  if (is.na(direction) || direction == "") direction <- "<"
  
  # pROC direction: "<" means cases have higher predictor; ">" means cases have lower predictor
  if (direction == "<") {
    pred_positive <- expr_vec >= cutoff_value
  } else {
    pred_positive <- expr_vec <= cutoff_value
  }
  
  actual_positive <- status_vec == 1
  mean(pred_positive != actual_positive)
}

split_patient_ids_by_expression <- function(
    expr_vec,
    patient_id_vec,
    cutoff_value,
    high_rule = ">"
) {
  expr_vec <- coerce_expression_numeric_vector(expr_vec, context = "expression vector for split")
  if (any(is.na(expr_vec))) stop("expr_vec contains NA values.", call. = FALSE)
  if (high_rule == ">") {
    high_mask <- expr_vec > cutoff_value
  } else if (high_rule == ">=") {
    high_mask <- expr_vec >= cutoff_value
  } else {
    stop("high_rule must be '>' or '>='.", call. = FALSE)
  }
  
  high_ids <- patient_id_vec[high_mask]
  low_ids  <- patient_id_vec[!high_mask]
  
  list(high_ids = high_ids, low_ids = low_ids)
}

split_by_expression_cutoff_objects <- function(
    expr_df,
    clin_df,
    feature_name,
    cutoff_value,
    patient_id_col = "patient_id",
    high_rule = ">"
) {
  if (!is.data.frame(expr_df) || !is.data.frame(clin_df)) {
    stop("expr_df and clin_df must be data.frames.", call. = FALSE)
  }
  if (!(patient_id_col %in% names(expr_df)) || !(patient_id_col %in% names(clin_df))) {
    stop("patient_id_col missing from inputs.", call. = FALSE)
  }
  if (!(feature_name %in% names(expr_df))) {
    stop("feature_name not found in expr_df.", call. = FALSE)
  }
  
  expr_ids <- expr_df[[patient_id_col]]
  clin_ids <- clin_df[[patient_id_col]]
  if (!setequal(expr_ids, clin_ids)) {
    stop("expr_df and clin_df must have identical patient_id sets.", call. = FALSE)
  }
  
  split_ids <- split_patient_ids_by_expression(
    expr_vec = expr_df[[feature_name]],
    patient_id_vec = expr_ids,
    cutoff_value = cutoff_value,
    high_rule = high_rule
  )
  
  high_expr_expr <- expr_df[expr_df[[patient_id_col]] %in% split_ids$high_ids, , drop = FALSE]
  low_expr_expr  <- expr_df[expr_df[[patient_id_col]] %in% split_ids$low_ids,  , drop = FALSE]
  high_expr_clin <- clin_df[clin_df[[patient_id_col]] %in% split_ids$high_ids, , drop = FALSE]
  low_expr_clin  <- clin_df[clin_df[[patient_id_col]] %in% split_ids$low_ids,  , drop = FALSE]
  
  list(
    high_expr_expr = high_expr_expr,
    low_expr_expr  = low_expr_expr,
    high_expr_clin = high_expr_clin,
    low_expr_clin  = low_expr_clin
  )
}

logrank_by_expression_groups <- function(
    high_expr_clin,
    low_expr_clin,
    time_col = "days_to_last_follow_up",
    status_col = "vital_status"
) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    warning("Package 'survival' is not available.", call. = FALSE)
    return(NA_real_)
  }
  
  if (!is.data.frame(high_expr_clin) || !is.data.frame(low_expr_clin)) {
    stop("Clinical inputs must be data.frames.", call. = FALSE)
  }
  if (!(time_col %in% names(high_expr_clin)) || !(status_col %in% names(high_expr_clin)) ||
      !(time_col %in% names(low_expr_clin))  || !(status_col %in% names(low_expr_clin))) {
    stop("Clinical inputs missing time/status columns.", call. = FALSE)
  }
  
  time_high <- suppressWarnings(as.numeric(high_expr_clin[[time_col]]))
  time_low  <- suppressWarnings(as.numeric(low_expr_clin[[time_col]]))
  status_high <- vital_status_to_event(high_expr_clin[[status_col]])
  status_low  <- vital_status_to_event(low_expr_clin[[status_col]])
  
  keep_h <- is.finite(time_high) & is.finite(status_high)
  keep_l <- is.finite(time_low) & is.finite(status_low)
  
  time_high <- time_high[keep_h]; status_high <- status_high[keep_h]
  time_low  <- time_low[keep_l];  status_low  <- status_low[keep_l]
  
  if (length(time_high) < 2L || length(time_low) < 2L) {
    warning("Not enough observations for log-rank test.", call. = FALSE)
    return(NA_real_)
  }
  
  time <- c(time_high, time_low)
  status <- c(status_high, status_low)
  group <- c(rep("high", length(time_high)), rep("low", length(time_low)))
  
  fit <- try(survival::survdiff(survival::Surv(time, status) ~ group), silent = TRUE)
  if (inherits(fit, "try-error")) {
    warning("Log-rank test failed.", call. = FALSE)
    return(NA_real_)
  }
  
  chisq <- as.numeric(fit$chisq)
  stats::pchisq(chisq, df = 1, lower.tail = FALSE)
}

compute_split_confusion_and_kappa <- function(
    high_survival_ids,
    low_survival_ids,
    high_expr_ids,
    low_expr_ids
) {
  universe <- intersect(
    c(high_survival_ids, low_survival_ids),
    c(high_expr_ids, low_expr_ids)
  )
  
  actual_high <- universe %in% high_survival_ids
  pred_high   <- universe %in% high_expr_ids
  
  TP <- sum(actual_high & pred_high)
  FP <- sum(!actual_high & pred_high)
  FN <- sum(actual_high & !pred_high)
  TN <- sum(!actual_high & !pred_high)
  
  confusion_matrix <- matrix(
    c(TP, FN, FP, TN),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Actual = c("HighSurvival", "LowSurvival"),
      Predicted = c("HighExpr", "LowExpr")
    )
  )
  
  accuracy <- if ((TP + FP + FN + TN) > 0) (TP + TN) / (TP + FP + FN + TN) else NA_real_
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  F1 <- if (is.finite(precision) && is.finite(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }
  Jaccard <- if ((TP + FP + FN) > 0) TP / (TP + FP + FN) else NA_real_
  
  total <- TP + FP + FN + TN
  po <- if (total > 0) (TP + TN) / total else NA_real_
  pe <- if (total > 0) {
    ((TP + FP) / total) * ((TP + FN) / total) +
      ((FN + TN) / total) * ((FP + TN) / total)
  } else {
    NA_real_
  }
  kappa <- if (is.finite(po) && is.finite(pe) && (1 - pe) != 0) (po - pe) / (1 - pe) else NA_real_
  
  list(
    confusion_matrix = confusion_matrix,
    counts = list(TP = TP, FP = FP, FN = FN, TN = TN),
    metrics = list(
      accuracy = accuracy,
      precision = precision,
      recall = recall,
      F1 = F1,
      Jaccard = Jaccard
    ),
    kappa = kappa
  )
}

compute_confusion_and_kappa_from_expr_tables <- function(
    high_survival_expr,
    low_survival_expr,
    high_expr_expr,
    low_expr_expr,
    patient_id_col = "patient_id"
) {
  if (!(patient_id_col %in% names(high_survival_expr)) ||
      !(patient_id_col %in% names(low_survival_expr)) ||
      !(patient_id_col %in% names(high_expr_expr)) ||
      !(patient_id_col %in% names(low_expr_expr))) {
    stop("patient_id_col missing from one or more inputs.", call. = FALSE)
  }
  
  compute_split_confusion_and_kappa(
    high_survival_ids = high_survival_expr[[patient_id_col]],
    low_survival_ids  = low_survival_expr[[patient_id_col]],
    high_expr_ids     = high_expr_expr[[patient_id_col]],
    low_expr_ids      = low_expr_expr[[patient_id_col]]
  )
}

run_all_feature_statistics_table <- function(
    high_survival_expr,
    low_survival_expr,
    patient_id_col = "patient_id",
    pseudo_count = 1
) {
  means <- compute_means_per_feature(
    high_survival_expr = high_survival_expr,
    low_survival_expr = low_survival_expr,
    patient_id_col = patient_id_col
  )
  diff_mean <- compute_diff_mean_per_feature(means$mean_high, means$mean_low)
  fold_change <- compute_fold_change_per_feature(means$mean_high, means$mean_low, pseudo_count = pseudo_count)
  log2FC <- compute_log2_fold_change_per_feature(fold_change)
  
  Ttest_P <- compute_ttest_per_feature(high_survival_expr, low_survival_expr, patient_id_col = patient_id_col)
  ks_out <- compute_ks_per_feature(high_survival_expr, low_survival_expr, patient_id_col = patient_id_col)
  Wilcox_P <- compute_wilcoxon_per_feature(high_survival_expr, low_survival_expr, patient_id_col = patient_id_col)
  
  data.frame(
    feature_coord = names(means$mean_high),
    mean_high = as.numeric(means$mean_high),
    mean_low = as.numeric(means$mean_low),
    diff_mean = as.numeric(diff_mean),
    fold_change = as.numeric(fold_change),
    log2FC = as.numeric(log2FC),
    Ttest_P = as.numeric(Ttest_P),
    KS_D = as.numeric(ks_out$KS_D),
    KS_P = as.numeric(ks_out$KS_P),
    Wilcox_P = as.numeric(Wilcox_P),
    stringsAsFactors = FALSE
  )
}

run_expression_defined_survival_suite_one_feature <- function(
    high_survival_expr,
    low_survival_expr,
    high_survival_clin,
    low_survival_clin,
    feature_name,
    patient_id_col = "patient_id",
    time_col = "days_to_last_follow_up",
    status_col = "vital_status",
    high_rule = ">"
) {
  split_inputs <- build_expression_split_inputs(
    high_survival_expr = high_survival_expr,
    low_survival_expr = low_survival_expr,
    high_survival_clin = high_survival_clin,
    low_survival_clin = low_survival_clin,
    patient_id_col = patient_id_col,
    time_col = time_col,
    status_col = status_col
  )
  
  expr_df <- split_inputs$expr_df
  clin_df <- split_inputs$clin_df
  group_vec <- split_inputs$group_vec
  
  if (!(feature_name %in% names(expr_df))) {
    stop("feature_name not found in expr_df.", call. = FALSE)
  }
  
  cox_res <- cox_fit_one_feature_stratified(
    expr_vec = expr_df[[feature_name]],
    clin_df = clin_df,
    group_vec = group_vec,
    time_col = time_col,
    status_col = status_col
  )
  
  roc_res <- roc_auc_and_cutoff_youden_per_feature(expr_df[[feature_name]], group_vec)
  
  split_objs <- split_by_expression_cutoff_objects(
    expr_df = expr_df,
    clin_df = clin_df,
    feature_name = feature_name,
    cutoff_value = roc_res$cutoff_value,
    patient_id_col = patient_id_col,
    high_rule = high_rule
  )
  
  logrank_p <- logrank_by_expression_groups(
    high_expr_clin = split_objs$high_expr_clin,
    low_expr_clin = split_objs$low_expr_clin,
    time_col = time_col,
    status_col = status_col
  )
  
  overlap <- compute_split_confusion_and_kappa(
    high_survival_ids = high_survival_expr[[patient_id_col]],
    low_survival_ids  = low_survival_expr[[patient_id_col]],
    high_expr_ids     = split_objs$high_expr_expr[[patient_id_col]],
    low_expr_ids      = split_objs$low_expr_expr[[patient_id_col]]
  )
  
  data.frame(
    feature_coord = feature_name,
    Cox_beta = cox_res$Cox_beta,
    Cox_HR = cox_res$Cox_HR,
    Cox_P = cox_res$Cox_P,
    Cox_n = cox_res$Cox_n,
    ROC_AUC = roc_res$ROC_AUC,
    cutoff_value = roc_res$cutoff_value,
    youden_J = roc_res$youden_J,
    LogRankExpr_P = logrank_p,
    TP = overlap$counts$TP,
    FP = overlap$counts$FP,
    FN = overlap$counts$FN,
    TN = overlap$counts$TN,
    accuracy = overlap$metrics$accuracy,
    precision = overlap$metrics$precision,
    recall = overlap$metrics$recall,
    F1 = overlap$metrics$F1,
    Jaccard = overlap$metrics$Jaccard,
    kappa = overlap$kappa,
    stringsAsFactors = FALSE
  )
}

run_expression_defined_survival_suite_all_features <- function(
    high_survival_expr,
    low_survival_expr,
    high_survival_clin,
    low_survival_clin,
    patient_id_col = "patient_id",
    time_col = "days_to_last_follow_up",
    status_col = "vital_status",
    high_rule = ">"
) {
  if (!is.data.frame(high_survival_expr) || !is.data.frame(low_survival_expr)) {
    stop("Expression inputs must be data.frames.", call. = FALSE)
  }
  if (!(patient_id_col %in% names(high_survival_expr)) ||
      !(patient_id_col %in% names(low_survival_expr))) {
    stop("patient_id_col missing from expression inputs.", call. = FALSE)
  }
  
  feature_cols <- setdiff(names(high_survival_expr), patient_id_col)
  out_list <- vector("list", length(feature_cols))
  
  for (i in seq_along(feature_cols)) {
    feat <- feature_cols[i]
    out_list[[i]] <- run_expression_defined_survival_suite_one_feature(
      high_survival_expr = high_survival_expr,
      low_survival_expr = low_survival_expr,
      high_survival_clin = high_survival_clin,
      low_survival_clin = low_survival_clin,
      feature_name = feat,
      patient_id_col = patient_id_col,
      time_col = time_col,
      status_col = status_col,
      high_rule = high_rule
    )
  }
  
  do.call(rbind, out_list)
}

run_expression_statistics_runner <- function(
    high_survival_expr,
    low_survival_expr,
    high_survival_clin,
    low_survival_clin,
    full_survival_expr,
    full_survival_clin,
    patient_id_col = "patient_id",
    time_col = "days_to_last_follow_up",
    status_col = "vital_status",
    pseudo_count = 1,
    high_rule = ">",
    max_p_value = 0.05,
    max_error_rate = 0.3
) {
  if (!(patient_id_col %in% names(high_survival_expr)) ||
      !(patient_id_col %in% names(low_survival_expr)) ||
      !(patient_id_col %in% names(full_survival_expr))) {
    stop("patient_id_col missing from expression inputs.", call. = FALSE)
  }
  if (!is.numeric(max_p_value) || length(max_p_value) != 1L || !is.finite(max_p_value)) {
    stop("max_p_value must be a single finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(max_error_rate) || length(max_error_rate) != 1L ||
      !is.finite(max_error_rate) || max_error_rate < 0 || max_error_rate > 1) {
    stop("max_error_rate must be a single numeric value between 0 and 1.", call. = FALSE)
  }
  
  drop_sample_id <- function(df) {
    if ("sample_id" %in% names(df)) df[["sample_id"]] <- NULL
    df
  }
  
  collapse_duplicate_features <- function(df, patient_id_col) {
    if (is.matrix(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (!is.data.frame(df)) stop("Expression input must be data.frame or matrix.", call. = FALSE)
    if (!(patient_id_col %in% names(df))) stop("patient_id_col missing.", call. = FALSE)
    feat_names <- setdiff(names(df), patient_id_col)
    dup_mask <- duplicated(feat_names) | duplicated(feat_names, fromLast = TRUE)
    if (!any(dup_mask)) return(df)
    dup_feats <- unique(feat_names[dup_mask])
    warning(sprintf("Duplicate feature names detected (%d). Collapsing by row-mean.", length(dup_feats)), call. = FALSE)
    warning(sprintf("Example duplicates: %s", paste(head(dup_feats, 10), collapse = ", ")), call. = FALSE)
    unique_feats <- unique(feat_names)
    out <- data.frame(df[[patient_id_col]], stringsAsFactors = FALSE)
    names(out)[1] <- patient_id_col
    for (feat in unique_feats) {
      cols <- which(names(df) == feat)
      if (length(cols) == 1L) {
        out[[feat]] <- df[[feat]]
      } else {
        mat <- as.matrix(df[, cols, drop = FALSE])
        storage.mode(mat) <- "numeric"
        out[[feat]] <- rowMeans(mat, na.rm = TRUE)
      }
    }
    out
  }
  
  # clean expression inputs
  high_survival_expr <- drop_sample_id(high_survival_expr)
  low_survival_expr  <- drop_sample_id(low_survival_expr)
  full_survival_expr <- drop_sample_id(full_survival_expr)
  
  high_survival_expr <- collapse_duplicate_features(high_survival_expr, patient_id_col)
  low_survival_expr  <- collapse_duplicate_features(low_survival_expr, patient_id_col)
  full_survival_expr <- collapse_duplicate_features(full_survival_expr, patient_id_col)
  
  high_survival_expr <- coerce_expression_df_numeric(high_survival_expr, patient_id_col = patient_id_col)
  low_survival_expr  <- coerce_expression_df_numeric(low_survival_expr,  patient_id_col = patient_id_col)
  full_survival_expr <- coerce_expression_df_numeric(full_survival_expr, patient_id_col = patient_id_col)
  
  # align full cohort expr/clin
  common_full <- intersect(full_survival_expr[[patient_id_col]], full_survival_clin[[patient_id_col]])
  full_survival_expr <- full_survival_expr[full_survival_expr[[patient_id_col]] %in% common_full, , drop = FALSE]
  full_survival_clin <- full_survival_clin[full_survival_clin[[patient_id_col]] %in% common_full, , drop = FALSE]
  
  # ensure no NA in time/status for full cohort
  keep_full <- !is.na(full_survival_clin[[time_col]]) & !is.na(full_survival_clin[[status_col]])
  full_survival_clin <- full_survival_clin[keep_full, , drop = FALSE]
  full_survival_expr <- full_survival_expr[full_survival_expr[[patient_id_col]] %in% full_survival_clin[[patient_id_col]], , drop = FALSE]
  
  feature_cols <- setdiff(names(high_survival_expr), patient_id_col)
  if (length(feature_cols) < 1L) stop("No feature columns found in expression inputs.", call. = FALSE)
  
  # Step 2: t-test, KS, Wilcoxon (high vs low survival)
  means <- compute_means_per_feature(high_survival_expr, low_survival_expr, patient_id_col = patient_id_col)
  diff_mean <- compute_diff_mean_per_feature(means$mean_high, means$mean_low)
  Ttest_P <- compute_ttest_per_feature(high_survival_expr, low_survival_expr, patient_id_col = patient_id_col)
  ks_out <- compute_ks_per_feature(high_survival_expr, low_survival_expr, patient_id_col = patient_id_col)
  Wilcoxon_P <- compute_wilcoxon_per_feature(high_survival_expr, low_survival_expr, patient_id_col = patient_id_col)
  
  stats_df <- data.frame(
    feature_coord = names(means$mean_high),
    mean_high = as.numeric(means$mean_high),
    mean_low = as.numeric(means$mean_low),
    diff_mean = as.numeric(diff_mean),
    Ttest_P = as.numeric(Ttest_P),
    KS_D = as.numeric(ks_out$KS_D),
    KS_P = as.numeric(ks_out$KS_P),
    Wilcoxon_P = as.numeric(Wilcoxon_P),
    stringsAsFactors = FALSE
  )
  
  # Step 3: keep at least 2/3 of (Ttest, Wilcoxon, KS)
  step3_df <- filter_significants_any_test(stats_df, max_p_value)
  
  # Step 4: Cox on full cohort (before survival split)
  cox_df <- cox_fit_all_features_unstratified(
    expr_df = full_survival_expr,
    clin_df = full_survival_clin,
    patient_id_col = patient_id_col,
    time_col = time_col,
    status_col = status_col
  )
  
  # Step 6: ROC/Youden on full cohort (event = vital_status)
  status_vec <- full_survival_clin[[status_col]]
  roc_df <- roc_auc_and_cutoff_youden_all_features_event(
    expr_df = full_survival_expr,
    status_vec = status_vec,
    patient_id_col = patient_id_col
  )
  
  # Step 7: error rate per feature using ROC cutoff
  error_rate <- rep(NA_real_, length(feature_cols))
  names(error_rate) <- feature_cols
  for (i in seq_along(feature_cols)) {
    feat <- feature_cols[i]
    cutoff_val <- roc_df$cutoff_value[roc_df$feature_coord == feat]
    direction_val <- roc_df$ROC_direction[roc_df$feature_coord == feat]
    error_rate[i] <- compute_error_rate_from_cutoff(
      expr_vec = full_survival_expr[[feat]],
      status_vec = status_vec,
      cutoff_value = cutoff_val,
      direction = direction_val
    )
  }
  error_df <- data.frame(
    feature_coord = feature_cols,
    error_rate = as.numeric(error_rate),
    stringsAsFactors = FALSE
  )

  # Candidates for log-rank:
  # must pass Step 3 and Step 8 (error rate threshold).
  step8_candidates <- feature_cols[
    feature_cols %in% step3_df$feature_coord &
      is.finite(error_rate[feature_cols]) &
      (error_rate[feature_cols] <= max_error_rate)
  ]
  
  # Step 9: Log-rank per Step-8 feature using ROC cutoff (full cohort)
  logrank_p <- rep(NA_real_, length(feature_cols))
  names(logrank_p) <- feature_cols
  for (feat in step8_candidates) {
    cutoff_val <- roc_df$cutoff_value[roc_df$feature_coord == feat]
    if (!is.finite(cutoff_val)) next
    split_objs <- split_by_expression_cutoff_objects(
      expr_df = full_survival_expr,
      clin_df = full_survival_clin,
      feature_name = feat,
      cutoff_value = cutoff_val,
      patient_id_col = patient_id_col,
      high_rule = high_rule
    )
    logrank_p[i] <- logrank_by_expression_groups(
      high_expr_clin = split_objs$high_expr_clin,
      low_expr_clin  = split_objs$low_expr_clin,
      time_col = time_col,
      status_col = status_col
    )
  }
  logrank_df <- data.frame(
    feature_coord = feature_cols,
    LogRank_P = as.numeric(logrank_p),
    stringsAsFactors = FALSE
  )
  
  # Merge all stats
  all_stats_df <- merge(stats_df, cox_df, by = "feature_coord", all = TRUE)
  all_stats_df <- merge(all_stats_df, roc_df, by = "feature_coord", all = TRUE)
  all_stats_df <- merge(all_stats_df, error_df, by = "feature_coord", all = TRUE)
  all_stats_df <- merge(all_stats_df, logrank_df, by = "feature_coord", all = TRUE)
  
  # Step 5: filter Cox_P
  step5_df <- all_stats_df[
    all_stats_df$feature_coord %in% step3_df$feature_coord &
      is.finite(all_stats_df$Cox_P) & (all_stats_df$Cox_P <= max_p_value),
    , drop = FALSE
  ]
  
  # Step 8: filter error rate on all Step-3 features (independent of Cox pass)
  step8_df <- all_stats_df[
    all_stats_df$feature_coord %in% step3_df$feature_coord &
      is.finite(all_stats_df$error_rate) & (all_stats_df$error_rate <= max_error_rate),
    , drop = FALSE
  ]
  
  # Step 10: filter log-rank p-value on Step-8 features
  step10_df <- all_stats_df[
    all_stats_df$feature_coord %in% step8_df$feature_coord &
      is.finite(all_stats_df$LogRank_P) & (all_stats_df$LogRank_P <= max_p_value),
    , drop = FALSE
  ]

  # Final keep rule:
  # - always keep Cox-significant features from Step 5
  # - also keep Step-10 survivors from the error/log-rank branch
  keep_features <- union(step5_df$feature_coord, step10_df$feature_coord)
  filtered_df <- all_stats_df[all_stats_df$feature_coord %in% keep_features, , drop = FALSE]
  
  list(
    all_stats_df = all_stats_df,
    filtered_df = filtered_df
  )
}
