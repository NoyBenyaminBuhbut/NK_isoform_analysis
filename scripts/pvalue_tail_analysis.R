#!/usr/bin/env Rscript

parse_args <- function(args) {
  defaults <- list(
    input = file.path("results", "presto", "pancan", "pancan_T1_vs_T3_presto.csv"),
    pval_col = "pval",
    out_dir = file.path("results", "tail_analysis", "pancan"),
    stability_fraction = 0.999
  )

  if (length(args) == 0L) {
    return(defaults)
  }

  out <- defaults
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--input", "--pval-col", "--out-dir", "--stability-fraction")) {
      if (i == length(args)) {
        stop("Missing value for argument: ", key, call. = FALSE)
      }
      value <- args[[i + 1L]]
      if (key == "--input") out$input <- value
      if (key == "--pval-col") out$pval_col <- value
      if (key == "--out-dir") out$out_dir <- value
      if (key == "--stability-fraction") out$stability_fraction <- as.numeric(value)
      i <- i + 2L
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
  }

  if (!is.finite(out$stability_fraction) || out$stability_fraction <= 0 || out$stability_fraction > 1) {
    stop("--stability-fraction must be in (0, 1].", call. = FALSE)
  }

  out
}

read_input_table <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }

  if (requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::fread(path, data.table = FALSE, check.names = FALSE, showProgress = FALSE))
  }

  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

prepare_input <- function(df, pval_col) {
  if (!(pval_col %in% names(df))) {
    stop("Missing p-value column: ", pval_col, call. = FALSE)
  }

  pvals <- suppressWarnings(as.numeric(df[[pval_col]]))
  keep <- is.finite(pvals) & (pvals >= 0) & (pvals <= 1)
  filtered_df <- df[keep, , drop = FALSE]
  filtered_pvals <- pvals[keep]

  if (length(filtered_pvals) < 10L) {
    stop("Need at least 10 finite p-values in [0, 1].", call. = FALSE)
  }

  list(df = filtered_df, pvals = filtered_pvals)
}

build_density_curve <- function(pvals) {
  dens <- stats::density(
    pvals,
    na.rm = TRUE,
    from = 0,
    to = max(pvals),
    n = 4096
  )

  list(
    x = dens$x,
    y = dens$y
  )
}

mirror_average_curve <- function(x_vals, y_vals, peak_x) {
  mirrored_y <- stats::approx(
    x = x_vals,
    y = y_vals,
    xout = (2 * peak_x) - x_vals,
    rule = 1
  )$y

  corrected_y <- rowMeans(cbind(y_vals, mirrored_y), na.rm = TRUE)
  corrected_y[!is.finite(corrected_y)] <- y_vals[!is.finite(corrected_y)]
  corrected_y
}

cumtrapz_from_zero <- function(x_vals, y_vals) {
  out <- numeric(length(x_vals))
  if (length(x_vals) < 2L) {
    return(out)
  }

  for (i in 2:length(x_vals)) {
    out[[i]] <- out[[i - 1L]] + ((y_vals[[i - 1L]] + y_vals[[i]]) / 2) * (x_vals[[i]] - x_vals[[i - 1L]])
  }
  out
}

compute_stability_fraction <- function(positive_flags) {
  n <- length(positive_flags)
  out <- numeric(n)
  running_pos <- 0L
  running_n <- 0L

  for (i in seq.int(from = n, to = 1L, by = -1L)) {
    running_n <- running_n + 1L
    if (isTRUE(positive_flags[[i]])) {
      running_pos <- running_pos + 1L
    }
    out[[i]] <- running_pos / running_n
  }

  out
}

find_stable_tail_threshold <- function(curve_x, diff_y, peak_x, stability_fraction) {
  tail_keep <- curve_x <= peak_x
  tail_x <- curve_x[tail_keep]
  tail_diff <- diff_y[tail_keep]

  if (length(tail_x) < 2L) {
    return(list(
      threshold_p = NA_real_,
      threshold_index = NA_integer_,
      cumulative_integral = numeric(0),
      stability = numeric(0)
    ))
  }

  cumulative_integral <- cumtrapz_from_zero(tail_x, tail_diff)
  positive_flags <- cumulative_integral > 0
  stability <- compute_stability_fraction(positive_flags)

  candidate_idx <- which(positive_flags & (stability >= stability_fraction))
  threshold_idx <- if (length(candidate_idx) < 1L) NA_integer_ else max(candidate_idx)
  threshold_p <- if (is.na(threshold_idx)) NA_real_ else tail_x[[threshold_idx]]

  list(
    threshold_p = threshold_p,
    threshold_index = threshold_idx,
    tail_x = tail_x,
    tail_diff = tail_diff,
    cumulative_integral = cumulative_integral,
    stability = stability
  )
}

raw_p_threshold_positions <- c(0.05, 0.01, 0.001)

interpolate_difference_at_points <- function(curve_x, diff_y, point_x) {
  stats::approx(curve_x, diff_y, xout = point_x, rule = 2)$y
}

build_tail_member_table <- function(df, pvals, pval_col, threshold_p, diff_curve_x, diff_curve_y, stability_fraction) {
  if (!is.finite(threshold_p)) {
    return(NULL)
  }

  point_diff <- interpolate_difference_at_points(diff_curve_x, diff_curve_y, pvals)
  member_mask <- is.finite(pvals) & (pvals <= threshold_p)
  if (!any(member_mask)) {
    return(NULL)
  }

  out_df <- df[member_mask, , drop = FALSE]
  out_df$stable_tail_threshold_p <- threshold_p
  out_df$stable_tail_diff_value <- point_diff[member_mask]
  out_df$stable_tail_rule <- paste0("cumulative_integral_positive_and_future_positive_fraction_ge_", stability_fraction)
  out_df <- out_df[order(suppressWarnings(as.numeric(out_df[[pval_col]])), na.last = TRUE), , drop = FALSE]
  rownames(out_df) <- NULL
  out_df
}

plot_full_curve <- function(curve_x, original_y, corrected_y, peak_x, threshold_p, out_file) {
  y_max <- max(c(original_y, corrected_y), na.rm = TRUE)
  grDevices::png(filename = out_file, width = 1400, height = 900, res = 150)

  graphics::plot(
    curve_x, original_y,
    type = "l",
    lwd = 3,
    col = "#0B5FFF",
    xlab = "p-value",
    ylab = "Density",
    main = "Stable significance-tail analysis (full)",
    ylim = c(0, y_max * 1.05)
  )
  graphics::lines(curve_x, corrected_y, lwd = 3, col = "#D1495B")
  graphics::abline(v = peak_x, lty = 2, lwd = 2, col = "#222222")
  if (is.finite(threshold_p)) {
    graphics::abline(v = threshold_p, lty = 3, lwd = 2, col = "#2A9D8F")
  }
  graphics::abline(v = raw_p_threshold_positions, lty = 3, col = "#999999")
  graphics::legend(
    "topright",
    legend = c("Original density", "Mirrored mean", "Peak", "Stable tail threshold", "p thresholds"),
    col = c("#0B5FFF", "#D1495B", "#222222", "#2A9D8F", "#999999"),
    lty = c(1, 1, 2, 3, 3),
    lwd = c(3, 3, 2, 2, 1),
    bty = "n"
  )
  grDevices::dev.off()
}

plot_zoom_curve <- function(curve_x, original_y, corrected_y, peak_x, threshold_p, out_file) {
  zoom_upper <- min(max(curve_x), max(0.05, if (is.finite(threshold_p)) threshold_p * 1.5 else 0.05))
  zoom_keep <- curve_x <= zoom_upper
  y_max <- max(c(original_y[zoom_keep], corrected_y[zoom_keep]), na.rm = TRUE)

  grDevices::png(filename = out_file, width = 1400, height = 900, res = 150)
  graphics::plot(
    curve_x, original_y,
    type = "l",
    lwd = 3,
    col = "#0B5FFF",
    xlab = "p-value",
    ylab = "Density",
    main = "Stable significance-tail analysis (tail zoom)",
    xlim = c(0, zoom_upper),
    ylim = c(0, y_max * 1.05)
  )
  graphics::lines(curve_x, corrected_y, lwd = 3, col = "#D1495B")
  graphics::abline(v = peak_x, lty = 2, lwd = 2, col = "#222222")
  if (is.finite(threshold_p)) {
    graphics::abline(v = threshold_p, lty = 3, lwd = 2, col = "#2A9D8F")
    shade_keep <- curve_x <= threshold_p & (original_y > corrected_y)
    idx <- which(shade_keep)
    if (length(idx) >= 2L) {
      split_groups <- split(idx, cumsum(c(1L, diff(idx) != 1L)))
      for (group_idx in split_groups) {
        if (length(group_idx) < 2L) next
        xs <- curve_x[group_idx]
        graphics::polygon(
          x = c(xs, rev(xs)),
          y = c(original_y[group_idx], rev(corrected_y[group_idx])),
          border = NA,
          col = "#F4A26180"
        )
      }
    }
  }
  graphics::abline(v = raw_p_threshold_positions, lty = 3, col = "#999999")
  grDevices::dev.off()
}

plot_residuals <- function(tail_x, tail_diff, cumulative_integral, stability, threshold_p, stability_fraction, out_file) {
  grDevices::png(filename = out_file, width = 1400, height = 1000, res = 150)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

  graphics::plot(
    tail_x, tail_diff,
    type = "l",
    lwd = 3,
    col = "#264653",
    xlab = "p-value",
    ylab = "Original - mirrored mean",
    main = "Residual density difference"
  )
  graphics::abline(h = 0, col = "#444444")
  if (is.finite(threshold_p)) {
    graphics::abline(v = threshold_p, lty = 3, lwd = 2, col = "#2A9D8F")
  }
  graphics::abline(v = raw_p_threshold_positions, lty = 3, col = "#999999")

  graphics::plot(
    tail_x, cumulative_integral,
    type = "l",
    lwd = 3,
    col = "#7A3E9D",
    xlab = "p-value",
    ylab = "Integral from 0 to p",
    main = paste0("Cumulative tail integral and stability (target ", stability_fraction, ")")
  )
  graphics::abline(h = 0, col = "#444444")
  if (is.finite(threshold_p)) {
    graphics::abline(v = threshold_p, lty = 3, lwd = 2, col = "#2A9D8F")
  }
  par(new = TRUE)
  graphics::plot(
    tail_x, stability,
    type = "l",
    lwd = 2,
    col = "#E76F51",
    axes = FALSE,
    xlab = "",
    ylab = "",
    ylim = c(0, 1)
  )
  graphics::axis(side = 4, col.axis = "#E76F51")
  graphics::mtext("Future positive fraction", side = 4, line = 2, col = "#E76F51")
  graphics::abline(h = stability_fraction, lty = 2, col = "#E76F51")
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

  input_res <- prepare_input(read_input_table(args$input), args$pval_col)
  df <- input_res$df
  pvals <- input_res$pvals

  curve <- build_density_curve(pvals)
  peak_idx <- which.max(curve$y)
  peak_x <- curve$x[[peak_idx]]
  corrected_y <- mirror_average_curve(curve$x, curve$y, peak_x)
  diff_y <- curve$y - corrected_y

  stable_res <- find_stable_tail_threshold(
    curve_x = curve$x,
    diff_vals = diff_y,
    peak_x = peak_x,
    stability_fraction = args$stability_fraction
  )

  threshold_p <- stable_res$threshold_p
  threshold_integral <- if (is.na(stable_res$threshold_index)) NA_real_ else stable_res$cumulative_integral[[stable_res$threshold_index]]
  threshold_stability <- if (is.na(stable_res$threshold_index)) NA_real_ else stable_res$stability[[stable_res$threshold_index]]
  tail_exists <- is.finite(threshold_p)

  analysis_id <- "stable_tail_density_raw_p"
  full_plot_file <- file.path(args$out_dir, paste0(analysis_id, "_full.png"))
  zoom_plot_file <- file.path(args$out_dir, paste0(analysis_id, "_zoom.png"))
  residual_plot_file <- file.path(args$out_dir, paste0(analysis_id, "_residual.png"))
  tail_member_file <- file.path(args$out_dir, paste0(analysis_id, "_tail_members.csv"))

  plot_full_curve(curve$x, curve$y, corrected_y, peak_x, threshold_p, full_plot_file)
  plot_zoom_curve(curve$x, curve$y, corrected_y, peak_x, threshold_p, zoom_plot_file)
  plot_residuals(stable_res$tail_x, stable_res$tail_diff, stable_res$cumulative_integral, stable_res$stability, threshold_p, args$stability_fraction, residual_plot_file)

  if (tail_exists) {
    tail_members <- build_tail_member_table(
      df = df,
      pvals = pvals,
      pval_col = args$pval_col,
      threshold_p = threshold_p,
      diff_curve_x = curve$x,
      diff_curve_y = diff_y,
      stability_fraction = args$stability_fraction
    )
    if (!is.null(tail_members)) {
      utils::write.csv(tail_members, tail_member_file, row.names = FALSE)
    }
  } else {
    tail_member_file <- NA_character_
  }

  summary_df <- data.frame(
    analysis_id = analysis_id,
    curve_method = "density",
    x_scale = "raw_p",
    n_features = length(pvals),
    peak_x = peak_x,
    stable_tail_threshold_p = threshold_p,
    cumulative_integral_at_threshold = threshold_integral,
    future_positive_fraction_at_threshold = threshold_stability,
    required_future_positive_fraction = args$stability_fraction,
    tail_exists = tail_exists,
    full_plot_file = full_plot_file,
    zoom_plot_file = zoom_plot_file,
    residual_plot_file = residual_plot_file,
    tail_member_file = tail_member_file,
    stringsAsFactors = FALSE
  )

  summary_file <- file.path(args$out_dir, "pvalue_tail_summary.csv")
  utils::write.csv(summary_df, summary_file, row.names = FALSE)

  message("Wrote summary: ", summary_file)
  message("Plot directory: ", args$out_dir)
  print(summary_df)
}

main()
