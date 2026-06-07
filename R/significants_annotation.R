# ============================================================
# Global parser (NO regex)
# Accepts: "chr:start-end" or "chr:start-end:strand"
# Strand default: "*"
# Output: data.frame(chr,start,end,strand)
# ============================================================
parse_coord_global <- function(x) {
  x <- as.character(x)
  parts <- strsplit(x, ":", fixed = TRUE)
  
  chr <- vapply(parts, function(p) if (length(p) >= 1) p[[1]] else NA_character_, character(1))
  rest <- vapply(parts, function(p) if (length(p) >= 2) p[[2]] else NA_character_, character(1))
  strand <- vapply(parts, function(p) if (length(p) >= 3) p[[3]] else "*", character(1))
  
  se <- strsplit(rest, "-", fixed = TRUE)
  start <- suppressWarnings(as.integer(vapply(se, function(p) if (length(p) >= 1) p[[1]] else NA_character_, character(1))))
  end   <- suppressWarnings(as.integer(vapply(se, function(p) if (length(p) >= 2) p[[2]] else NA_character_, character(1))))
  
  data.frame(chr = chr, start = start, end = end, strand = strand, stringsAsFactors = FALSE)
}

# ============================================================
# 1) annotate_features_by_exon_overlap (NO regex)
# Adds: gene_id, exon_rank, exon_coord
# Replicates rows when a feature overlaps multiple exons
# ============================================================
annotate_features_by_exon_overlap <- function(
    out_filt,
    biomart_exon_df,
    feature_coord_col = "feature_coord",
    exon_coord_col    = "exon_coord",
    gene_id_col       = "gene_id",
    exon_rank_col     = "exon_rank_in_transcript",
    ignore_strand     = TRUE
) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Packages GenomicRanges, IRanges, S4Vectors are required.", call. = FALSE)
  }
  
  # guards
  stopifnot(is.data.frame(out_filt), is.data.frame(biomart_exon_df))
  if (!(feature_coord_col %in% names(out_filt))) stop("Missing column in out_filt: ", feature_coord_col, call. = FALSE)
  req_ex <- c(exon_coord_col, gene_id_col)
  miss <- setdiff(req_ex, names(biomart_exon_df))
  if (length(miss) > 0) stop("Missing columns in biomart_exon_df: ", paste(miss, collapse = ", "), call. = FALSE)
  has_exon_rank <- exon_rank_col %in% names(biomart_exon_df)

  # No significant features to annotate: return an empty table with expected columns.
  if (nrow(out_filt) == 0L) {
    out0 <- out_filt[0, , drop = FALSE]
    out0$gene_id <- character(0)
    out0$exon_rank <- integer(0)
    out0$exon_coord <- character(0)
    return(out0)
  }
  
  # parse coords (NO regex)
  feat_tbl <- parse_coord_global(out_filt[[feature_coord_col]])
  exon_tbl <- parse_coord_global(biomart_exon_df[[exon_coord_col]])
  
  # validate parsed coords (basic)
  if (any(is.na(feat_tbl$chr) | is.na(feat_tbl$start) | is.na(feat_tbl$end))) {
    stop("Failed to parse feature_coord for some rows (expected chr:start-end[:strand]).", call. = FALSE)
  }
  if (any(is.na(exon_tbl$chr) | is.na(exon_tbl$start) | is.na(exon_tbl$end))) {
    stop("Failed to parse exon_coord for some rows (expected chr:start-end[:strand]).", call. = FALSE)
  }
  
  # build GRanges
  feat_gr <- GenomicRanges::GRanges(
    seqnames = feat_tbl$chr,
    ranges   = IRanges::IRanges(start = feat_tbl$start, end = feat_tbl$end),
    strand   = if (ignore_strand) rep("*", nrow(feat_tbl)) else feat_tbl$strand
  )
  
  exon_gr <- GenomicRanges::GRanges(
    seqnames = exon_tbl$chr,
    ranges   = IRanges::IRanges(start = exon_tbl$start, end = exon_tbl$end),
    strand   = rep("*", nrow(exon_tbl))
  )
  
  hits <- GenomicRanges::findOverlaps(feat_gr, exon_gr, ignore.strand = TRUE)
  
  if (length(hits) == 0L) {
    out0 <- out_filt[0, , drop = FALSE]
    out0$gene_id <- character(0)
    out0$exon_rank <- integer(0)
    out0$exon_coord <- character(0)
    return(out0)
  }
  
  q <- S4Vectors::queryHits(hits)
  s <- S4Vectors::subjectHits(hits)
  
  out_rep  <- out_filt[q, , drop = FALSE]
  exon_sub <- biomart_exon_df[s, , drop = FALSE]
  
  # ensure exon_coord always has strand; if missing -> add ":*"
  exon_coord <- as.character(exon_sub[[exon_coord_col]])
  exon_tbl2 <- parse_coord_global(exon_coord)
  exon_coord <- paste0(exon_tbl2$chr, ":", exon_tbl2$start, "-", exon_tbl2$end, ":", exon_tbl2$strand)
  
  out_rep$gene_id   <- as.character(exon_sub[[gene_id_col]])
  out_rep$exon_rank <- if (has_exon_rank) suppressWarnings(as.integer(exon_sub[[exon_rank_col]])) else NA_integer_
  out_rep$exon_coord <- exon_coord
  
  rownames(out_rep) <- NULL
  out_rep
}

# ============================================================
# 2) add_exon_overlap_metrics (NO regex)
# Adds:
#   exon_bp_in_feature_frac
#   exon_feature_illustration
# ============================================================
add_exon_overlap_metrics <- function(
    df,
    feature_coord_col = "feature_coord",
    exon_coord_col    = "exon_coord"
) {
  stopifnot(is.data.frame(df))
  if (!(feature_coord_col %in% names(df))) stop("Missing feature coord column: ", feature_coord_col, call. = FALSE)
  if (!(exon_coord_col %in% names(df))) stop("Missing exon coord column: ", exon_coord_col, call. = FALSE)

  if (nrow(df) == 0L) {
    df$exon_bp_in_feature_frac <- numeric(0)
    df$exon_feature_illustration <- character(0)
    return(df)
  }
  
  f <- parse_coord_global(df[[feature_coord_col]])
  e <- parse_coord_global(df[[exon_coord_col]])
  
  ov_start <- pmax(f$start, e$start, na.rm = TRUE)
  ov_end   <- pmin(f$end,   e$end,   na.rm = TRUE)
  ov_len   <- pmax(0L, ov_end - ov_start + 1L)
  
  exon_len <- pmax(0L, e$end - e$start + 1L)
  df$exon_bp_in_feature_frac <- ov_len / exon_len
  
  left_mark  <- ifelse(f$start < e$start, "-", ifelse(f$start > e$start, "!", "*"))
  right_mark <- ifelse(f$end   > e$end,   "-", ifelse(f$end   < e$end,   "!", "*"))
  df$exon_feature_illustration <- paste0(left_mark, "█ █ █", right_mark)
  
  df
}
# ============================================================
# separate experiment fron control
# negative control can have an empty data frame
# ============================================================

split_results_by_gene_type <- function(
    out_annot,
    genes_config = "config/genes.csv",
    gene_id_col = "gene_id",
    type_col    = "type"
) {
  # ---- sanity: out_annot ----
  if (!is.data.frame(out_annot)) {
    stop("out_annot must be a data.frame", call. = FALSE)
  }
  if (!gene_id_col %in% names(out_annot)) {
    stop(sprintf("Missing '%s' column in out_annot", gene_id_col), call. = FALSE)
  }
  
  genes_df <- normalize_gene_type_map(genes_config)
  
  # ---- sanity: genes_df ----
  if (!gene_id_col %in% names(genes_df)) {
    stop(sprintf("Missing '%s' column in genes_config", gene_id_col), call. = FALSE)
  }
  if (!type_col %in% names(genes_df)) {
    stop(sprintf("Missing '%s' column in genes_config", type_col), call. = FALSE)
  }
  
  # ---- normalize ----
  out_annot[[gene_id_col]] <- as.character(out_annot[[gene_id_col]])
  genes_df[[gene_id_col]]  <- as.character(genes_df[[gene_id_col]])
  genes_df[[type_col]]     <- as.character(genes_df[[type_col]])
  
  # ---- map gene_id -> type ----
  type_map <- setNames(genes_df[[type_col]], genes_df[[gene_id_col]])
  
  missing_genes <- setdiff(unique(out_annot[[gene_id_col]]), names(type_map))
  if (length(missing_genes) > 0L) {
    stop(
      sprintf(
        "The following gene_id values are missing from genes_config: %s",
        paste(missing_genes, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # ---- attach type ----
  out_annot$type <- type_map[out_annot[[gene_id_col]]]
  
  # ---- split (empty is OK) ----
  out_experiment <- out_annot[out_annot$type == "experiment", , drop = FALSE]
  out_positive   <- out_annot[out_annot$type == "positive",   , drop = FALSE]
  out_negative   <- out_annot[out_annot$type == "negative",   , drop = FALSE]
  out_normalization <- out_annot[out_annot$type == "normalization", , drop = FALSE]
  
  out_experiment$type <- NULL
  out_positive$type   <- NULL
  out_negative$type   <- NULL
  out_normalization$type <- NULL
  
  list(
    experiment = out_experiment,
    positive   = out_positive,
    negative   = out_negative,
    normalization = out_normalization
  )
}
save_significants_by_type <- function(
    split_results,
    cohort_id,
    base_dir = "results/significants/cohorts"
) {
  # ---- sanity ----
  if (!is.list(split_results)) {
    stop("split_results must be a list (output of split_results_by_gene_type)", call. = FALSE)
  }
  
  if (missing(cohort_id) || length(cohort_id) != 1L || is.na(cohort_id) || !nzchar(cohort_id)) {
    stop("cohort_id must be a single non-empty string", call. = FALSE)
  }
  
  # expected types (keys of the list)
  allowed_types <- c("experiment", "positive", "negative", "normalization")
  bad_types <- setdiff(names(split_results), allowed_types)
  if (length(bad_types) > 0L) {
    stop(
      sprintf("Unexpected result types: %s", paste(bad_types, collapse = ", ")),
      call. = FALSE
    )
  }
  
  # ---- create cohort directory ----
  cohort_dir <- file.path(base_dir, cohort_id)
  if (!dir.exists(cohort_dir)) {
    dir.create(cohort_dir, recursive = TRUE)
  }
  
  # ---- save each type ----
  for (type in names(split_results)) {
    df <- split_results[[type]]
    
    if (!is.data.frame(df)) {
      stop(sprintf("Result for type '%s' is not a data.frame", type), call. = FALSE)
    }
    
    # empty data.frame is allowed → still write header-only CSV
    out_file <- file.path(
      cohort_dir,
      sprintf("%s_%s_sign.csv", cohort_id, type)
    )
    
    write.csv(
      df,
      file      = out_file,
      row.names = FALSE,
      quote     = TRUE
    )
  }
  
  invisible(TRUE)
}
