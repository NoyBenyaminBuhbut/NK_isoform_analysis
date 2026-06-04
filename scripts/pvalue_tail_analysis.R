#!/usr/bin/env Rscript

parse_args <- function(args) {
  defaults <- list(
    input = file.path("results", "presto", "pancan", "pancan_T1_vs_T3_presto.csv"),
    pval_col = "pval",
    out_dir = file.path("results", "tail_analysis", "pancan")
  )

  if (length(args) == 0L) {
    return(defaults)
  }

  out <- defaults
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--input", "--pval-col", "--out-dir")) {
      if (i == length(args)) {
        stop("Missing value for argument: ", key, call. = FALSE)
      }
      value <- args[[i + 1L]]
      if (key == "--input") out$input <- value
      if (key == "--pval-col") out$pval_col <- value
      if (key == "--out-dir") out$out_dir <- value
      i <- i + 2L
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
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

prepare_pvalues <- function(df, pval_col) {
  if (!(pval_col %in% names(df))) {
    stop("Missing p-value column: ", pval_col, call. = FALSE)
  }

  pvals <- suppressWarnings(as.numeric(df[[pval_col]]))
  keep <- is.finite(pvals) & (pvals >= 0) & (pvals <= 1)
  pvals <- pvals[keep]

  if (length(pvals) < 10L) {
    stop("Need at least 10 finite p-values in [0, 1].", call. = FALSE)
  }

  pvals
}

transform_pvalues <- function(pvals, x_scale) {
  if (x_scale == "raw_p") {
    return(pvals)
  }
  if (x_scale == "neglog10_p") {
    return(-log10(pmax(pvals, .Machine$double.xmin)))
  }
  stop("Unsupported x_scale: ", x_scale, call. = FALSE)
}

build_histogram_curve <- function(x_vals) {
  hist_obj <- graphics::hist(x_vals, breaks = "FD", plot = FALSE, include.lowest = TRUE, right = TRUE)
  widths <- diff(hist_obj$breaks)
  y_vals <- hist_obj$counts / sum(hist_obj$counts)

  list(
    x = hist_obj$mids,
    y = y_vals,
    widths = widths,
    method = "histogram",
    x_min = min(hist_obj$breaks),
    x_max = max(hist_obj$breaks)
  )
}

build_density_curve <- function(x_vals) {
  dens <- stats::density(x_vals, na.rm = TRUE, from = min(x_vals), to = max(x_vals), n = 2048)
  list(
    x = dens$x,
    y = dens$y,
    widths = c(diff(dens$x), tail(diff(dens$x), 1L)),
    method = "density",
    x_min = min(dens$x),
    x_max = max(dens$x)
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

find_tail_crossing <- function(x_vals, diff_vals, peak_x, x_scale) {
  tail_keep <- if (x_scale == "raw_p") x_vals <= peak_x else x_vals >= peak_x
  tail_x <- x_vals[tail_keep]
  tail_diff <- diff_vals[tail_keep]

  if (length(tail_x) < 2L) {
    return(NA_real_)
  }

  crossings <- numeric(0)
  for (i in seq_len(length(tail_x) - 1L)) {
    x1 <- tail_x[[i]]
    x2 <- tail_x[[i + 1L]]
    d1 <- tail_diff[[i]]
    d2 <- tail_diff[[i + 1L]]

    if (!is.finite(d1) || !is.finite(d2)) {
      next
    }

    if (d1 == 0) {
      crossings <- c(crossings, x1)
      next
    }

    if (d1 * d2 < 0 || d2 == 0) {
      crossing <- x1 + (0 - d1) * (x2 - x1) / (d2 - d1)
      crossings <- c(crossings, crossing)
    }
  }

  if (length(crossings) < 1L) {
    return(NA_real_)
  }

  if (x_scale == "raw_p") {
    return(max(crossings))
  }

  min(crossings)
}

tail_integral <- function(x_vals, diff_vals, crossing_x, x_scale) {
  if (!is.finite(crossing_x)) {
    return(NA_real_)
  }

  tail_keep <- if (x_scale == "raw_p") x_vals <= crossing_x else x_vals >= crossing_x
  tail_x <- x_vals[tail_keep]
  tail_diff <- diff_vals[tail_keep]

  ord <- order(tail_x)
  tail_x <- tail_x[ord]
  tail_diff <- tail_diff[ord]

  if (length(tail_x) < 2L) {
    return(0)
  }

  sum(diff(tail_x) * (tail_diff[-length(tail_diff)] + tail_diff[-1L]) / 2)
}

crossing_to_raw_p <- function(crossing_x, x_scale) {
  if (!is.finite(crossing_x)) {
    return(NA_real_)
  }
  if (x_scale == "raw_p") {
    return(crossing_x)
  }
  10^(-crossing_x)
}

plot_curve_analysis <- function(curve, corrected_y, peak_x, crossing_x, analysis_label, x_scale, out_file) {
  y_max <- max(c(curve$y, corrected_y), na.rm = TRUE)
  graphics::png(filename = out_file, width = 1400, height = 900, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::plot(
    curve$x, curve$y,
    type = "l",
    lwd = 3,
    col = "#0B5FFF",
    xlab = if (x_scale == "raw_p") "p-value" else expression(-log[10](p)),
    ylab = "Frequency / density",
    main = analysis_label,
    ylim = c(0, y_max * 1.05)
  )
  graphics::lines(curve$x, corrected_y, lwd = 3, col = "#D1495B")
  graphics::abline(v = peak_x, lty = 2, lwd = 2, col = "#222222")

  if (is.finite(crossing_x)) {
    graphics::abline(v = crossing_x, lty = 3, lwd = 2, col = "#2A9D8F")
  }

  graphics::legend(
    "topright",
    legend = c("Original", "Mirrored mean", "Peak", "Tail crossing"),
    col = c("#0B5FFF", "#D1495B", "#222222", "#2A9D8F"),
    lty = c(1, 1, 2, 3),
    lwd = c(3, 3, 2, 2),
    bty = "n"
  )
}

analyze_curve <- function(pvals, x_scale, curve_method) {
  x_vals <- transform_pvalues(pvals, x_scale)
  curve <- if (curve_method == "histogram") build_histogram_curve(x_vals) else build_density_curve(x_vals)

  peak_idx <- which.max(curve$y)
  peak_x <- curve$x[[peak_idx]]
  corrected_y <- mirror_average_curve(curve$x, curve$y, peak_x)
  diff_y <- curve$y - corrected_y

  crossing_x <- find_tail_crossing(curve$x, diff_y, peak_x, x_scale)
  integral_value <- tail_integral(curve$x, diff_y, crossing_x, x_scale)
  tail_exists <- is.finite(integral_value) && (integral_value > 0)

  list(
    curve = curve,
    corrected_y = corrected_y,
    peak_x = peak_x,
    crossing_x = crossing_x,
    crossing_raw_p = crossing_to_raw_p(crossing_x, x_scale),
    tail_integral = integral_value,
    tail_exists = tail_exists,
    x_scale = x_scale,
    curve_method = curve_method
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

  df <- read_input_table(args$input)
  pvals <- prepare_pvalues(df, args$pval_col)

  analyses <- list(
    analyze_curve(pvals, x_scale = "raw_p", curve_method = "histogram"),
    analyze_curve(pvals, x_scale = "neglog10_p", curve_method = "histogram"),
    analyze_curve(pvals, x_scale = "raw_p", curve_method = "density"),
    analyze_curve(pvals, x_scale = "neglog10_p", curve_method = "density")
  )

  summary_df <- do.call(rbind, lapply(analyses, function(res) {
    analysis_id <- paste(res$curve_method, res$x_scale, sep = "_")

    plot_file <- file.path(args$out_dir, paste0(analysis_id, ".png"))
    plot_curve_analysis(
      curve = res$curve,
      corrected_y = res$corrected_y,
      peak_x = res$peak_x,
      crossing_x = res$crossing_x,
      analysis_label = paste("Tail analysis:", analysis_id),
      x_scale = res$x_scale,
      out_file = plot_file
    )

    data.frame(
      analysis_id = analysis_id,
      curve_method = res$curve_method,
      x_scale = res$x_scale,
      n_features = length(pvals),
      peak_x = res$peak_x,
      crossing_x = res$crossing_x,
      crossing_raw_p = res$crossing_raw_p,
      tail_integral = res$tail_integral,
      tail_exists = res$tail_exists,
      plot_file = plot_file,
      stringsAsFactors = FALSE
    )
  }))

  summary_file <- file.path(args$out_dir, "pvalue_tail_summary.csv")
  utils::write.csv(summary_df, summary_file, row.names = FALSE)

  message("Wrote summary: ", summary_file)
  message("Plot directory: ", args$out_dir)
  print(summary_df)
}

main()
