#!/usr/bin/env Rscript

script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_path_arg) > 0L) sub("^--file=", "", script_path_arg[[1L]]) else file.path("scripts", "prepare_isoform_normalization.R")
script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE)
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

source(file.path(repo_root, "R", "utility.R"))

ensure_biomart_annotation_available <- function(
    genes_config = "config/genes.csv",
    out_dir = file.path("intermediate", "biomart")
) {
  biomart_override <- Sys.getenv("BIOMART_ISOFORM_INFO_RDA", unset = "")
  if (nzchar(biomart_override) && file.exists(path.expand(biomart_override))) {
    return(path.expand(biomart_override))
  }

  expected_files <- c(
    file.path(out_dir, "NK_genes_info.RDA"),
    file.path(out_dir, "NK_exons_info.RDA"),
    file.path(out_dir, "NK_isoform_info.RDA"),
    file.path(out_dir, "mapped_exons_to_isoforms.RDA")
  )

  if (all(file.exists(expected_files))) {
    return(expected_files)
  }

  if (!is.character(genes_config) || length(genes_config) != 1L) {
    stop(
      "genes_config must be a path string when biomart annotation needs to be generated locally.",
      call. = FALSE
    )
  }
  if (!file.exists(genes_config)) {
    stop("genes_config file not found: ", genes_config, call. = FALSE)
  }

  biomart_script <- file.path(repo_root, "R", "biomart_hg19_info.R")
  if (!file.exists(biomart_script)) {
    stop("Biomart export script not found: ", biomart_script, call. = FALSE)
  }

  source(biomart_script)
  if (!exists("export_biomart_tables_hg19", mode = "function")) {
    stop("export_biomart_tables_hg19() was not loaded from: ", biomart_script, call. = FALSE)
  }

  export_biomart_tables_hg19(
    genes_csv_path = genes_config,
    out_dir = out_dir
  )

  if (!all(file.exists(expected_files))) {
    stop(
      "Biomart export completed but required annotation files are still missing from: ",
      out_dir,
      call. = FALSE
    )
  }

  write_manifest_file(
    folder_path = out_dir,
    manifest_entries = list(
      genes_config = normalizePath(genes_config, winslash = "/", mustWork = TRUE),
      output_files = normalizePath(expected_files, winslash = "/", mustWork = TRUE),
      transformation = "Exported BioMart GRCh37 gene, exon, isoform, and exon-to-isoform annotation tables for isoform-to-gene mapping.",
      created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    ),
    file_name = "manifest.txt"
  )

  expected_files
}

normalize_feature_core <- function(feature_ids) {
  sub("\\..*$", "", as.character(feature_ids))
}

cohort_label_from_id <- function(cohort_id) {
  if (tolower(as.character(cohort_id)) == "pancan") "pancan" else toupper(as.character(cohort_id))
}

normalization_label_from_id <- function(normalization_gene) {
  gene_label <- as.character(normalization_gene)
  if (!nzchar(gene_label)) {
    stop("normalization_gene must be a non-empty string.", call. = FALSE)
  }
  if (tolower(gene_label) == "asitself") {
    return("asitself")
  }
  toupper(gene_label)
}

required_gene_ids_for_normalization <- function(filtered_map_df, normalization_gene_label) {
  if (normalization_gene_label == "asitself") {
    unique(filtered_map_df$gene_id)
  } else {
    normalization_gene_label
  }
}

candidate_gene_roots <- function() {
  c(
    file.path("data", "cohorts"),
    file.path("GDCdata", "cohorts"),
    path.expand("~/data/cohorts")
  )
}

discover_pancan_gene_expression_files <- function() {
  roots <- unique(candidate_gene_roots())
  roots <- roots[dir.exists(roots)]
  if (length(roots) < 1L) {
    stop(
      "Could not find cohort gene-expression roots for pancan. Tried:\n",
      paste(candidate_gene_roots(), collapse = "\n"),
      call. = FALSE
    )
  }

  gene_files <- character(0)
  for (root_dir in roots) {
    cohort_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
    if (length(cohort_dirs) < 1L) {
      next
    }
    for (cohort_dir in cohort_dirs) {
      genes_dir <- file.path(cohort_dir, "genes")
      if (dir.exists(genes_dir)) {
        gene_files <- c(gene_files, get_single_file_in_dir(genes_dir))
      }
    }
  }

  gene_files <- unique(gene_files[file.exists(gene_files)])
  if (length(gene_files) < 1L) {
    stop("No pancan cohort gene-expression files were found.", call. = FALSE)
  }

  gene_files
}

validate_no_duplicate_ids <- function(ids, label) {
  ids <- as.character(ids)
  dup_ids <- unique(ids[duplicated(ids)])
  if (length(dup_ids) > 0L) {
    stop(
      "Duplicated ", label, ": ",
      paste(head(dup_ids, 10L), collapse = ", "),
      call. = FALSE
    )
  }
}

candidate_gene_matrix_from_feature_rows <- function(gene_df, needed_gene_ids) {
  gene_ids <- trimws(as.character(gene_df[[1L]]))
  needed_gene_ids_upper <- toupper(as.character(needed_gene_ids))
  sample_ids_raw <- names(gene_df)[-1L]
  sample_ids_canonical <- canonicalize_matrix_sample_ids(sample_ids_raw)

  validate_no_duplicate_ids(sample_ids_canonical, "gene-expression sample IDs after canonicalization")

  mat <- as.matrix(gene_df[, -1L, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- gene_ids

  list(
    orientation = "feature_rows",
    gene_ids = gene_ids,
    sample_ids_raw = sample_ids_raw,
    sample_ids_canonical = sample_ids_canonical,
    matrix = mat,
    needed_gene_hits = sum(unique(needed_gene_ids_upper) %in% toupper(gene_ids), na.rm = TRUE)
  )
}

candidate_gene_matrix_from_sample_rows <- function(gene_df, needed_gene_ids) {
  sample_ids_raw <- trimws(as.character(gene_df[[1L]]))
  sample_ids_canonical <- canonicalize_matrix_sample_ids(sample_ids_raw)
  validate_no_duplicate_ids(sample_ids_canonical, "gene-expression sample IDs after canonicalization")

  gene_ids <- names(gene_df)[-1L]
  needed_gene_ids_upper <- toupper(as.character(needed_gene_ids))
  mat <- t(as.matrix(gene_df[, -1L, drop = FALSE]))
  storage.mode(mat) <- "numeric"
  rownames(mat) <- gene_ids
  colnames(mat) <- sample_ids_raw

  list(
    orientation = "sample_rows",
    gene_ids = gene_ids,
    sample_ids_raw = sample_ids_raw,
    sample_ids_canonical = sample_ids_canonical,
    matrix = mat,
    needed_gene_hits = sum(unique(needed_gene_ids_upper) %in% toupper(gene_ids), na.rm = TRUE)
  )
}

choose_gene_matrix_candidate <- function(gene_df, target_sample_ids_canonical, needed_gene_ids) {
  candidates <- list(
    candidate_gene_matrix_from_feature_rows(gene_df, needed_gene_ids),
    candidate_gene_matrix_from_sample_rows(gene_df, needed_gene_ids)
  )

  scores <- vapply(candidates, function(candidate) {
    sample_overlap <- sum(candidate$sample_ids_canonical %in% target_sample_ids_canonical, na.rm = TRUE)
    sample_overlap * 100000L + candidate$needed_gene_hits
  }, numeric(1))

  best_idx <- which.max(scores)
  best_candidate <- candidates[[best_idx]]
  best_overlap <- sum(best_candidate$sample_ids_canonical %in% target_sample_ids_canonical, na.rm = TRUE)
  if (!is.finite(best_overlap) || best_overlap < 1L) {
    stop(
      "Could not orient gene-expression table to a sample-compatible matrix. ",
      "No overlapping sample IDs were detected."
      , call. = FALSE
    )
  }

  best_candidate
}

standardize_gene_expression_matrix <- function(gene_file, target_sample_ids_canonical, needed_gene_ids) {
  gene_df <- read_cluster_table(gene_file)
  if (!is.data.frame(gene_df) || ncol(gene_df) < 2L) {
    stop("Gene-expression table must have at least two columns: ", gene_file, call. = FALSE)
  }

  candidate <- choose_gene_matrix_candidate(gene_df, target_sample_ids_canonical, needed_gene_ids)

  gene_ids <- trimws(as.character(candidate$gene_ids))
  gene_ids_upper <- toupper(gene_ids)
  needed_gene_ids_chr <- as.character(needed_gene_ids)
  needed_gene_ids_upper <- toupper(needed_gene_ids_chr)
  needed_gene_map <- setNames(needed_gene_ids_chr, needed_gene_ids_upper)
  keep_gene_rows <- !is.na(gene_ids) & nzchar(gene_ids) & gene_ids_upper %in% needed_gene_ids_upper
  if (!any(keep_gene_rows)) {
    return(list(
      file = gene_file,
      orientation = candidate$orientation,
      gene_matrix = candidate$matrix[0, , drop = FALSE],
      sample_ids_raw = candidate$sample_ids_raw,
      sample_ids_canonical = candidate$sample_ids_canonical
    ))
  }

  gene_matrix <- candidate$matrix[keep_gene_rows, , drop = FALSE]
  gene_ids <- gene_ids[keep_gene_rows]
  gene_ids_upper <- gene_ids_upper[keep_gene_rows]
  rownames(gene_matrix) <- unname(needed_gene_map[gene_ids_upper])

  validate_no_duplicate_ids(rownames(gene_matrix), paste0("gene IDs in gene-expression file ", gene_file))

  list(
    file = gene_file,
    orientation = candidate$orientation,
    gene_matrix = gene_matrix,
    sample_ids_raw = candidate$sample_ids_raw,
    sample_ids_canonical = candidate$sample_ids_canonical
  )
}

combine_gene_expression_matrices <- function(standardized_gene_matrices, needed_gene_ids) {
  standardized_gene_matrices <- standardized_gene_matrices[vapply(standardized_gene_matrices, function(x) ncol(x$gene_matrix) > 0L, logical(1))]
  if (length(standardized_gene_matrices) < 1L) {
    stop("No standardized gene-expression matrices were available to combine.", call. = FALSE)
  }

  sample_ids_canonical_all <- unlist(lapply(standardized_gene_matrices, `[[`, "sample_ids_canonical"), use.names = FALSE)
  validate_no_duplicate_ids(sample_ids_canonical_all, "combined gene-expression sample IDs after canonicalization")

  combined <- matrix(
    NA_real_,
    nrow = length(needed_gene_ids),
    ncol = length(sample_ids_canonical_all),
    dimnames = list(needed_gene_ids, sample_ids_canonical_all)
  )

  col_offset <- 0L
  sample_ids_raw_all <- character(length(sample_ids_canonical_all))
  sample_ids_canonical_out <- character(length(sample_ids_canonical_all))
  source_files <- character(length(sample_ids_canonical_all))

  for (standardized in standardized_gene_matrices) {
    gene_matrix <- standardized$gene_matrix
    n_cols <- ncol(gene_matrix)
    idx <- seq.int(col_offset + 1L, col_offset + n_cols)

    if (nrow(gene_matrix) > 0L) {
      combined[rownames(gene_matrix), idx] <- gene_matrix
    }

    sample_ids_raw_all[idx] <- standardized$sample_ids_raw
    sample_ids_canonical_out[idx] <- standardized$sample_ids_canonical
    source_files[idx] <- standardized$file
    col_offset <- col_offset + n_cols
  }

  list(
    gene_matrix = combined,
    sample_ids_raw = sample_ids_raw_all,
    sample_ids_canonical = sample_ids_canonical_out,
    source_files = source_files
  )
}

load_required_gene_expression_matrix <- function(cohort_label, isoform_sample_ids_canonical, needed_gene_ids) {
  if (tolower(cohort_label) == "pancan") {
    gene_files <- discover_pancan_gene_expression_files()
  } else {
    gene_files <- resolve_gene_matrix_file(cohort_label)
  }

  standardized <- lapply(gene_files, function(gene_file) {
    standardize_gene_expression_matrix(
      gene_file = gene_file,
      target_sample_ids_canonical = isoform_sample_ids_canonical,
      needed_gene_ids = needed_gene_ids
    )
  })

  combined <- combine_gene_expression_matrices(standardized, needed_gene_ids)
  overlap_count <- sum(combined$sample_ids_canonical %in% isoform_sample_ids_canonical, na.rm = TRUE)
  if (overlap_count < 1L) {
    stop("No overlapping sample IDs were found between isoform and gene-expression matrices.", call. = FALSE)
  }

  list(
    gene_matrix = combined$gene_matrix,
    sample_ids_raw = combined$sample_ids_raw,
    sample_ids_canonical = combined$sample_ids_canonical,
    source_files = unique(gene_files)
  )
}

align_gene_expression_matrix_to_isoforms <- function(gene_expr, isoform_sample_ids_raw) {
  isoform_sample_ids_canonical <- canonicalize_matrix_sample_ids(isoform_sample_ids_raw)
  gene_col_idx <- match(isoform_sample_ids_canonical, gene_expr$sample_ids_canonical)

  aligned_gene_matrix <- matrix(
    NA_real_,
    nrow = nrow(gene_expr$gene_matrix),
    ncol = length(isoform_sample_ids_raw),
    dimnames = list(rownames(gene_expr$gene_matrix), isoform_sample_ids_raw)
  )

  matched_samples <- !is.na(gene_col_idx)
  if (any(matched_samples)) {
    aligned_gene_matrix[, matched_samples] <- gene_expr$gene_matrix[, gene_col_idx[matched_samples], drop = FALSE]
  }

  list(
    aligned_gene_matrix = aligned_gene_matrix,
    isoform_sample_ids_raw = isoform_sample_ids_raw,
    isoform_sample_ids_canonical = isoform_sample_ids_canonical,
    matched_sample_count = sum(matched_samples),
    unmatched_sample_ids = isoform_sample_ids_raw[!matched_samples]
  )
}

load_and_filter_isoform_expression <- function(cohort_label, genes_config = "config/genes.csv") {
  ensure_biomart_annotation_available(genes_config = genes_config)
  isoform_matrix_file <- resolve_isoform_matrix_file(cohort_label)
  isoform_expr_df <- read_cluster_table(isoform_matrix_file)
  if (!is.data.frame(isoform_expr_df) || ncol(isoform_expr_df) < 2L) {
    stop("Isoform-expression table must have at least two columns.", call. = FALSE)
  }

  feature_ids <- as.character(isoform_expr_df[[1L]])
  if (anyNA(feature_ids) || any(!nzchar(feature_ids))) {
    stop("Isoform-expression table contains missing feature identifiers.", call. = FALSE)
  }

  isoform_info <- load_biomart_isoform_info()
  genes_df <- read_genes_config(genes_config)
  allowed_gene_ids <- unique(genes_df$gene_id)

  feature_core <- normalize_feature_core(feature_ids)
  map_idx <- match(feature_core, isoform_info$isoform_core)
  mapped_gene_id <- isoform_info$gene_id[map_idx]
  mapped_isoform_id <- isoform_info$isoform_id[map_idx]
  keep_rows <- !is.na(mapped_gene_id) & mapped_gene_id %in% allowed_gene_ids

  if (!any(keep_rows)) {
    stop("No isoform rows matched genes_config after transcript-to-gene mapping.", call. = FALSE)
  }

  filtered_expr_df <- isoform_expr_df[keep_rows, , drop = FALSE]
  filtered_map_df <- data.frame(
    feature_id = feature_ids[keep_rows],
    feature_core = feature_core[keep_rows],
    gene_id = mapped_gene_id[keep_rows],
    isoform_id = mapped_isoform_id[keep_rows],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (anyNA(filtered_map_df$gene_id) || any(!nzchar(filtered_map_df$gene_id))) {
    stop("Filtered isoform rows contain missing parent-gene mappings.", call. = FALSE)
  }

  list(
    isoform_matrix_file = isoform_matrix_file,
    isoform_expr_df = filtered_expr_df,
    filtered_map_df = filtered_map_df,
    biomart_isoform_info_file = resolve_biomart_isoform_info_file(),
    genes_config_file = if (is.character(genes_config) && length(genes_config) == 1L) genes_config else NA_character_
  )
}

convert_log2_tpm_df_to_tpm <- function(expr_df) {
  feature_ids <- as.character(expr_df[[1L]])
  value_matrix <- as.matrix(expr_df[, -1L, drop = FALSE])
  storage.mode(value_matrix) <- "numeric"
  converted <- (2^value_matrix) - 0.001
  dimnames(converted) <- list(feature_ids, names(expr_df)[-1L])

  converted_df <- data.frame(
    feature_id = feature_ids,
    converted,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  list(
    expr_df = converted_df,
    value_matrix = converted
  )
}

normalize_by_single_gene <- function(isoform_tpm_matrix, aligned_gene_matrix, normalization_gene_label) {
  gene_ids_upper <- toupper(rownames(aligned_gene_matrix))
  gene_match_idx <- match(toupper(normalization_gene_label), gene_ids_upper)
  if (is.na(gene_match_idx)) {
    stop(
      "Normalization gene was not found in the gene-expression matrix: ",
      normalization_gene_label,
      call. = FALSE
    )
  }

  denominator <- aligned_gene_matrix[gene_match_idx, ]
  denominator[!is.finite(denominator) | denominator <= 0] <- NA_real_

  normalized_matrix <- sweep(isoform_tpm_matrix, 2L, denominator, "/")
  normalized_matrix[!is.finite(normalized_matrix)] <- NA_real_

  list(
    normalized_matrix = normalized_matrix,
    denominator_gene = rownames(aligned_gene_matrix)[gene_match_idx],
    missing_denominator_samples = colnames(isoform_tpm_matrix)[is.na(denominator)]
  )
}

normalize_as_itself <- function(isoform_tpm_matrix, filtered_map_df, aligned_gene_matrix) {
  gene_ids_upper <- toupper(rownames(aligned_gene_matrix))
  feature_gene_upper <- toupper(filtered_map_df$gene_id)
  gene_row_idx <- match(feature_gene_upper, gene_ids_upper)

  if (anyNA(filtered_map_df$gene_id) || any(!nzchar(filtered_map_df$gene_id))) {
    stop("Every filtered isoform must have a parent gene for 'asitself' normalization.", call. = FALSE)
  }

  normalized_matrix <- matrix(
    NA_real_,
    nrow = nrow(isoform_tpm_matrix),
    ncol = ncol(isoform_tpm_matrix),
    dimnames = dimnames(isoform_tpm_matrix)
  )

  missing_parent_gene_ids <- unique(filtered_map_df$gene_id[is.na(gene_row_idx)])
  present_gene_ids <- unique(filtered_map_df$gene_id[!is.na(gene_row_idx)])

  for (gene_id in present_gene_ids) {
    gene_rows <- which(filtered_map_df$gene_id == gene_id)
    gene_idx <- gene_row_idx[gene_rows[[1L]]]
    denominator <- aligned_gene_matrix[gene_idx, ]
    valid_denominator <- is.finite(denominator) & denominator > 0

    if (any(valid_denominator)) {
      normalized_matrix[gene_rows, valid_denominator] <-
        isoform_tpm_matrix[gene_rows, valid_denominator, drop = FALSE] /
        rep(denominator[valid_denominator], each = length(gene_rows))
    }
  }

  normalized_matrix[!is.finite(normalized_matrix)] <- NA_real_

  list(
    normalized_matrix = normalized_matrix,
    missing_parent_gene_ids = missing_parent_gene_ids
  )
}

save_prepared_matrix_artifact <- function(folder_path, file_name, expr_df, filtered_map_df, metadata) {
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  }

  output_file <- file.path(folder_path, file_name)
  saveRDS(
    list(
      expr_df = expr_df,
      feature_map = filtered_map_df,
      metadata = metadata
    ),
    output_file
  )
  output_file
}

prepare_isoform_normalization <- function(
    cohort_id,
    normalization_gene,
    genes_config = "config/genes.csv",
    intermediate_root = "intermediate"
) {
  cohort_label <- cohort_label_from_id(cohort_id)
  normalization_label <- normalization_label_from_id(normalization_gene)

  isoform_loaded <- load_and_filter_isoform_expression(cohort_label, genes_config = genes_config)
  isoform_sample_ids_raw <- names(isoform_loaded$isoform_expr_df)[-1L]
  isoform_sample_ids_canonical <- canonicalize_matrix_sample_ids(isoform_sample_ids_raw)
  validate_no_duplicate_ids(isoform_sample_ids_canonical, "isoform-expression sample IDs after canonicalization")

  converted <- convert_log2_tpm_df_to_tpm(isoform_loaded$isoform_expr_df)

  isoform_tpm_dir <- file.path(intermediate_root, "isoform_tpm", cohort_label)
  isoform_tpm_metadata <- list(
    cohort_id = cohort_label,
    normalization_gene = normalization_label,
    input_files = c(
      normalizePath(isoform_loaded$isoform_matrix_file, winslash = "/", mustWork = TRUE),
      normalizePath(isoform_loaded$biomart_isoform_info_file, winslash = "/", mustWork = TRUE),
      if (!is.na(isoform_loaded$genes_config_file)) normalizePath(isoform_loaded$genes_config_file, winslash = "/", mustWork = TRUE) else "config/genes.csv"
    ),
    transformation = "Filtered isoform matrix to genes in config/genes.csv, then converted values from log2(TPM + 0.001) to TPM using TPM = 2^value - 0.001.",
    output_file = file.path(isoform_tpm_dir, paste0(cohort_label, "_isoform_tpm_filtered.rds")),
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  isoform_tpm_output_file <- save_prepared_matrix_artifact(
    folder_path = isoform_tpm_dir,
    file_name = paste0(cohort_label, "_isoform_tpm_filtered.rds"),
    expr_df = converted$expr_df,
    filtered_map_df = isoform_loaded$filtered_map_df,
    metadata = isoform_tpm_metadata
  )
  isoform_tpm_manifest <- write_manifest_file(isoform_tpm_dir, isoform_tpm_metadata)

  needed_gene_ids <- required_gene_ids_for_normalization(isoform_loaded$filtered_map_df, normalization_label)
  gene_expr <- load_required_gene_expression_matrix(
    cohort_label = cohort_label,
    isoform_sample_ids_canonical = isoform_sample_ids_canonical,
    needed_gene_ids = needed_gene_ids
  )
  aligned_gene <- align_gene_expression_matrix_to_isoforms(gene_expr, isoform_sample_ids_raw)

  if (aligned_gene$matched_sample_count < 1L) {
    stop("No sample IDs overlapped between isoform and gene-expression matrices after alignment.", call. = FALSE)
  }

  if (normalization_label == "asitself") {
    normalized <- normalize_as_itself(
      isoform_tpm_matrix = converted$value_matrix,
      filtered_map_df = isoform_loaded$filtered_map_df,
      aligned_gene_matrix = aligned_gene$aligned_gene_matrix
    )
    transformation_description <- paste(
      "Filtered isoform TPM values were divided sample-wise by the TPM of each isoform's mapped parent gene.",
      "Values with missing, zero, or negative parent-gene TPM denominators were set to NA."
    )
    missing_gene_ids <- normalized$missing_parent_gene_ids
  } else {
    normalized <- normalize_by_single_gene(
      isoform_tpm_matrix = converted$value_matrix,
      aligned_gene_matrix = aligned_gene$aligned_gene_matrix,
      normalization_gene_label = normalization_label
    )
    transformation_description <- paste(
      "Filtered isoform TPM values were divided sample-wise by the TPM of normalization gene",
      normalized$denominator_gene, ".",
      "Values with missing, zero, or negative normalization-gene TPM denominators were set to NA."
    )
    missing_gene_ids <- character(0)
  }

  normalized_expr_df <- data.frame(
    feature_id = isoform_loaded$filtered_map_df$feature_id,
    normalized$normalized_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  normalization_dir <- file.path(intermediate_root, "normalization", normalization_label, cohort_label)
  normalization_metadata <- list(
    cohort_id = cohort_label,
    normalization_gene = normalization_label,
    input_files = c(
      normalizePath(isoform_loaded$isoform_matrix_file, winslash = "/", mustWork = TRUE),
      normalizePath(isoform_tpm_output_file, winslash = "/", mustWork = TRUE),
      normalizePath(isoform_loaded$biomart_isoform_info_file, winslash = "/", mustWork = TRUE),
      normalizePath(gene_expr$source_files, winslash = "/", mustWork = TRUE)
    ),
    transformation = transformation_description,
    output_file = file.path(normalization_dir, paste0(cohort_label, "_normalized_isoform_matrix.rds")),
    unmatched_isoform_samples = aligned_gene$unmatched_sample_ids,
    missing_gene_denominators = missing_gene_ids,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  normalization_output_file <- save_prepared_matrix_artifact(
    folder_path = normalization_dir,
    file_name = paste0(cohort_label, "_normalized_isoform_matrix.rds"),
    expr_df = normalized_expr_df,
    filtered_map_df = isoform_loaded$filtered_map_df,
    metadata = normalization_metadata
  )
  normalization_manifest <- write_manifest_file(normalization_dir, normalization_metadata)

  list(
    cohort_id = cohort_label,
    normalization_gene = normalization_label,
    isoform_matrix_file = isoform_loaded$isoform_matrix_file,
    gene_expression_files = gene_expr$source_files,
    biomart_isoform_info_file = isoform_loaded$biomart_isoform_info_file,
    isoform_tpm_output_file = isoform_tpm_output_file,
    isoform_tpm_manifest = isoform_tpm_manifest,
    normalization_output_file = normalization_output_file,
    normalization_manifest = normalization_manifest,
    filtered_map_df = isoform_loaded$filtered_map_df,
    isoform_tpm_expr_df = converted$expr_df,
    normalized_expr_df = normalized_expr_df
  )
}

run_prepare_isoform_normalization_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) {
    stop(
      "Usage: Rscript scripts/prepare_isoform_normalization.R <cohort_id> <normalization_gene>",
      call. = FALSE
    )
  }

  cohort_id <- args[[1L]]
  normalization_gene <- args[[2L]]
  result <- prepare_isoform_normalization(cohort_id = cohort_id, normalization_gene = normalization_gene)

  message("Prepared normalized isoform matrix: ", result$normalization_output_file)
  invisible(result)
}

if (!isTRUE(getOption("nk_isoform_analysis.skip_prepare_isoform_normalization_cli", FALSE))) {
  run_prepare_isoform_normalization_cli()
}
