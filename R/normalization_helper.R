source("R/utility.R")

nh_default_cohort_data_roots <- function() {
  c(
    file.path("data", "cohorts"),
    file.path("GDCdata", "cohorts"),
    path.expand("~/data/cohorts")
  )
}

nh_canonicalize_sample_ids <- function(sample_ids) {
  sample_ids <- trimws(as.character(sample_ids))
  matched <- grepl("^(([^-]+-){3}[0-9]{2}).*$", sample_ids, perl = TRUE)
  canonical <- sample_ids
  canonical[matched] <- sub("^((?:[^-]+-){3}[0-9]{2}).*$", "\\1", sample_ids[matched], perl = TRUE)
  fallback <- !matched & nchar(sample_ids) >= 15L
  canonical[fallback] <- substr(sample_ids[fallback], 1L, 15L)
  canonical
}

nh_assert_sample_level_ids <- function(sample_ids) {
  sample_ids <- trimws(as.character(sample_ids))
  patient_level_tcga <- grepl("^TCGA-[^-]+-[^-]+$", sample_ids)
  if (any(patient_level_tcga)) {
    bad_ids <- unique(sample_ids[patient_level_tcga])
    stop(
      "Isoform matrix contains TCGA patient IDs instead of sample IDs: ",
      paste(head(bad_ids, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

nh_detect_table_sep <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)

  first_line <- readLines(con, n = 1L, warn = FALSE)
  if (length(first_line) < 1L || !nzchar(first_line[[1L]])) {
    stop("Input table is empty: ", path, call. = FALSE)
  }

  counts <- c(
    tab = lengths(regmatches(first_line, gregexpr("\t", first_line, fixed = TRUE))),
    comma = lengths(regmatches(first_line, gregexpr(",", first_line, fixed = TRUE))),
    semicolon = lengths(regmatches(first_line, gregexpr(";", first_line, fixed = TRUE)))
  )

  c(tab = "\t", comma = ",", semicolon = ";")[[names(which.max(counts))]]
}

nh_read_table <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }

  if (requireNamespace("data.table", quietly = TRUE) &&
      (!grepl("\\.gz$", path, ignore.case = TRUE) || requireNamespace("R.utils", quietly = TRUE))) {
    return(data.table::fread(path, data.table = FALSE, check.names = FALSE, showProgress = FALSE))
  }

  sep <- nh_detect_table_sep(path)
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)

  utils::read.table(
    con,
    sep = sep,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = "\""
  )
}

nh_read_genes_config <- function(genes_config = "config/genes.csv") {
  if (!file.exists(genes_config)) {
    stop("Missing genes config file: ", genes_config, call. = FALSE)
  }

  genes_df <- utils::read.csv(genes_config, stringsAsFactors = FALSE, check.names = FALSE)
  assert_required_columns(genes_df, c("gene_id", "type"), df_name = "genes_config")
  genes_df$gene_id <- trimws(as.character(genes_df$gene_id))
  genes_df$type <- trimws(as.character(genes_df$type))
  genes_df <- genes_df[nzchar(genes_df$gene_id), , drop = FALSE]
  genes_df
}

nh_read_gene_symbol_aliases <- function(alias_csv = "config/gene_symbol_aliases.csv") {
  if (!file.exists(alias_csv)) {
    return(data.frame(
      query_symbol = character(0),
      resolved_symbol = character(0),
      ensembl_gene_id = character(0),
      assembly = character(0),
      source = character(0),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  alias_df <- utils::read.csv(alias_csv, stringsAsFactors = FALSE, check.names = FALSE)
  assert_required_columns(
    alias_df,
    c("query_symbol", "resolved_symbol"),
    df_name = "gene_symbol_aliases"
  )

  alias_df$query_symbol <- trimws(as.character(alias_df$query_symbol))
  alias_df$resolved_symbol <- trimws(as.character(alias_df$resolved_symbol))
  if ("assembly" %in% names(alias_df)) {
    alias_df$assembly <- trimws(as.character(alias_df$assembly))
  } else {
    alias_df$assembly <- NA_character_
  }

  alias_df <- alias_df[
    nzchar(alias_df$query_symbol) & nzchar(alias_df$resolved_symbol),
    ,
    drop = FALSE
  ]
  alias_df
}

nh_resolve_gene_symbols_against_matrix <- function(
    requested_symbols,
    available_symbols,
    assembly = "GRCh37",
    alias_csv = "config/gene_symbol_aliases.csv"
) {
  requested_symbols <- trimws(as.character(requested_symbols))
  available_symbols <- trimws(as.character(available_symbols))

  resolved_symbols <- requested_symbols
  used_alias <- rep(FALSE, length(requested_symbols))

  alias_df <- nh_read_gene_symbol_aliases(alias_csv = alias_csv)
  if (nrow(alias_df) < 1L) {
    return(list(
      resolved_symbols = resolved_symbols,
      used_alias = used_alias
    ))
  }

  if (!is.null(assembly) && nzchar(as.character(assembly)) && "assembly" %in% names(alias_df)) {
    keep_assembly <- is.na(alias_df$assembly) | !nzchar(alias_df$assembly) | toupper(alias_df$assembly) == toupper(as.character(assembly))
    alias_df <- alias_df[keep_assembly, , drop = FALSE]
  }

  if (nrow(alias_df) < 1L) {
    return(list(
      resolved_symbols = resolved_symbols,
      used_alias = used_alias
    ))
  }

  available_upper <- toupper(available_symbols)

  for (idx in seq_along(requested_symbols)) {
    req <- requested_symbols[[idx]]
    if (!nzchar(req)) {
      next
    }
    if (toupper(req) %in% available_upper) {
      next
    }

    alias_rows <- alias_df[toupper(alias_df$query_symbol) == toupper(req), , drop = FALSE]
    if (nrow(alias_rows) < 1L) {
      next
    }

    candidates <- unique(alias_rows$resolved_symbol)
    candidate_match <- candidates[toupper(candidates) %in% available_upper]
    if (length(candidate_match) == 1L) {
      resolved_symbols[[idx]] <- candidate_match[[1L]]
      used_alias[[idx]] <- TRUE
    }
  }

  list(
    resolved_symbols = resolved_symbols,
    used_alias = used_alias
  )
}

nh_load_biomart_isoform_info <- function(
    biomart_isoform_rda = file.path("intermediate", "biomart", "NK_isoform_info.RDA")
) {
  biomart_file <- path.expand(as.character(biomart_isoform_rda[[1L]]))
  if (!file.exists(biomart_file)) {
    stop("Missing biomart isoform annotation file: ", biomart_file, call. = FALSE)
  }

  biomart_env <- new.env(parent = emptyenv())
  load(biomart_file, envir = biomart_env)
  if (!exists("NK_isoform_info", envir = biomart_env, inherits = FALSE)) {
    stop("NK_isoform_info object is missing from: ", biomart_file, call. = FALSE)
  }

  biomart_isoform_df <- get("NK_isoform_info", envir = biomart_env, inherits = FALSE)
  assert_required_columns(biomart_isoform_df, c("gene_id", "isoform_id"), df_name = "NK_isoform_info")
  biomart_isoform_df
}

nh_extract_isoform_core <- function(feature_ids) {
  feature_ids <- trimws(as.character(feature_ids))
  has_enst <- grepl("ENST[0-9]+", feature_ids)
  out_ids <- feature_ids
  out_ids[has_enst] <- sub(".*(ENST[0-9]+).*", "\\1", feature_ids[has_enst])
  out_ids[!has_enst] <- sub("\\..*$", "", feature_ids[!has_enst])
  out_ids
}

nh_resolve_cohort_file <- function(
    cohort_id,
    subdirs,
    cohort_data_roots = nh_default_cohort_data_roots(),
    file_pattern = NULL,
    preferred_patterns = NULL
) {
  cohort_id <- trimws(as.character(cohort_id))
  if (!nzchar(cohort_id)) {
    stop("cohort_id must be a non-empty string.", call. = FALSE)
  }

  cohort_tokens <- unique(c(cohort_id, toupper(cohort_id), tolower(cohort_id)))
  candidate_roots <- unique(path.expand(as.character(cohort_data_roots)))
  candidate_dirs <- unique(unlist(lapply(cohort_tokens, function(tok) {
    unlist(lapply(candidate_roots, function(root_dir) {
      file.path(root_dir, tok, subdirs)
    }), use.names = FALSE)
  })))

  for (dir_path in candidate_dirs) {
    if (dir.exists(dir_path)) {
      entries <- list.files(
        path = dir_path,
        pattern = file_pattern,
        all.files = FALSE,
        full.names = TRUE,
        recursive = FALSE,
        include.dirs = FALSE,
        no.. = TRUE
      )

      if (length(entries) < 1L) {
        next
      }

      if (length(entries) == 1L) {
        return(entries[[1L]])
      }

      if (!is.null(preferred_patterns)) {
        for (pat in preferred_patterns) {
          preferred <- grep(pat, basename(entries), value = FALSE)
          if (length(preferred) == 1L) {
            return(entries[[preferred]])
          }
        }
      }

      stop(
        "Expected exactly one file in directory but found ", length(entries), " files: ", dir_path, "\n",
        paste0(" - ", basename(entries), collapse = "\n"),
        call. = FALSE
      )
    }
  }

  stop(
    "Could not find a cohort input directory for cohort_id='", cohort_id,
    "' in subdir(s): ", paste(subdirs, collapse = ", "),
    "\nTried:\n", paste(candidate_dirs, collapse = "\n"),
    call. = FALSE
  )
}

nh_get_calculated_tpm_dir <- function(
    cohort_id,
    calculated_tpm_root = file.path("intermediate", "calculated_TPM", "cohort")
) {
  file.path(calculated_tpm_root, toupper(trimws(as.character(cohort_id))))
}

nh_get_calculated_tpm_file <- function(
    cohort_id,
    calculated_tpm_root = file.path("intermediate", "calculated_TPM", "cohort")
) {
  file.path(
    nh_get_calculated_tpm_dir(cohort_id, calculated_tpm_root = calculated_tpm_root),
    paste0(toupper(trimws(as.character(cohort_id))), "_isoform_TPM.csv")
  )
}

nh_read_prepared_tpm_matrix <- function(
    cohort_id,
    calculated_tpm_root = file.path("intermediate", "calculated_TPM", "cohort")
) {
  tpm_file <- nh_get_calculated_tpm_file(cohort_id, calculated_tpm_root = calculated_tpm_root)
  if (!file.exists(tpm_file)) {
    stop("Prepared TPM matrix does not exist: ", tpm_file, call. = FALSE)
  }
  nh_read_table(tpm_file)
}

nh_filter_isoform_matrix_by_config <- function(
    expr_df,
    genes_config = "config/genes.csv",
    biomart_isoform_rda = file.path("intermediate", "biomart", "NK_isoform_info.RDA")
) {
  if (!is.data.frame(expr_df) || ncol(expr_df) < 2L) {
    stop("expr_df must contain one feature column and at least one sample column.", call. = FALSE)
  }

  genes_df <- nh_read_genes_config(genes_config = genes_config)
  biomart_isoform_df <- nh_load_biomart_isoform_info(biomart_isoform_rda = biomart_isoform_rda)

  biomart_isoform_df$gene_id <- trimws(as.character(biomart_isoform_df$gene_id))
  biomart_isoform_df$isoform_id <- trimws(as.character(biomart_isoform_df$isoform_id))
  biomart_isoform_df$isoform_core <- nh_extract_isoform_core(biomart_isoform_df$isoform_id)
  biomart_isoform_df <- biomart_isoform_df[!duplicated(biomart_isoform_df$isoform_core), , drop = FALSE]

  allowed_gene_ids <- unique(genes_df$gene_id)
  feature_ids <- as.character(expr_df[[1L]])
  feature_cores <- nh_extract_isoform_core(feature_ids)
  mapped_gene_ids <- biomart_isoform_df$gene_id[match(feature_cores, biomart_isoform_df$isoform_core)]

  keep <- !is.na(mapped_gene_ids) & nzchar(mapped_gene_ids) & mapped_gene_ids %in% allowed_gene_ids
  filtered_df <- expr_df[keep, , drop = FALSE]
  if (nrow(filtered_df) < 1L) {
    stop("No isoforms remained after filtering with genes config and biomart annotation.", call. = FALSE)
  }

  attr(filtered_df, "mapped_gene_ids") <- mapped_gene_ids[keep]
  filtered_df
}

nh_read_gene_matrix <- function(
    cohort_id,
    cohort_data_roots = nh_default_cohort_data_roots()
) {
  gene_file <- nh_resolve_cohort_file(
    cohort_id = cohort_id,
    subdirs = "genes",
    cohort_data_roots = cohort_data_roots,
    file_pattern = "\\.(csv|tsv|txt|gz)$",
    preferred_patterns = c("_gene_expr\\.tsv\\.gz$", "_gene_expr\\.csv\\.gz$", "\\.tsv\\.gz$", "\\.csv\\.gz$")
  )
  gene_df <- nh_read_table(gene_file)
  list(file = gene_file, data = gene_df)
}

nh_write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, path, row.names = FALSE)
  invisible(path)
}

convert_log2_tpm_to_tpm <- function(
    cohort_id,
    cohort_data_roots = nh_default_cohort_data_roots(),
    calculated_tpm_root = file.path("intermediate", "calculated_TPM", "cohort"),
    pseudocount = 0.001,
    feature_col = 1L
) {
  cohort_id <- trimws(as.character(cohort_id))
  if (!nzchar(cohort_id)) {
    stop("cohort_id must be a non-empty string.", call. = FALSE)
  }
  if (tolower(cohort_id) == "pancan") {
    stop("convert_log2_tpm_to_tpm() does not support pancan.", call. = FALSE)
  }

  isoform_file <- nh_resolve_cohort_file(
    cohort_id = cohort_id,
    subdirs = c("isoform", "isoforms"),
    cohort_data_roots = cohort_data_roots
  )
  expr_df <- nh_read_table(isoform_file)
  nh_assert_sample_level_ids(names(expr_df)[-feature_col])

  if (!is.data.frame(expr_df)) {
    stop("expr_df must be a data.frame.", call. = FALSE)
  }
  if (ncol(expr_df) < 2L) {
    stop("expr_df must contain one feature column and at least one sample column.", call. = FALSE)
  }
  if (!is.numeric(feature_col) || length(feature_col) != 1L || feature_col < 1L || feature_col > ncol(expr_df)) {
    stop("feature_col must point to a valid column in expr_df.", call. = FALSE)
  }

  out_df <- expr_df
  sample_cols <- setdiff(seq_len(ncol(out_df)), feature_col)
  for (col_idx in sample_cols) {
    values <- suppressWarnings(as.numeric(out_df[[col_idx]]))
    if (all(is.na(values))) {
      stop(
        "All values became NA after numeric coercion for sample column: ",
        names(out_df)[[col_idx]],
        call. = FALSE
      )
    }
    out_df[[col_idx]] <- (2^values) - pseudocount
  }

  output_file <- nh_get_calculated_tpm_file(cohort_id, calculated_tpm_root = calculated_tpm_root)
  nh_write_csv(out_df, output_file)

  attr(out_df, "isoform_file") <- isoform_file
  attr(out_df, "output_file") <- output_file
  out_df
}

normalize_isoform_by_gene <- function(
    cohort_id,
    normalization_gene,
    cohort_data_roots = nh_default_cohort_data_roots(),
    calculated_tpm_root = file.path("intermediate", "calculated_TPM", "cohort"),
    genes_config = "config/genes.csv",
    biomart_isoform_rda = file.path("intermediate", "biomart", "NK_isoform_info.RDA"),
    normalization_root = file.path("intermediate", "cohorts"),
    gene_symbol_alias_csv = "config/gene_symbol_aliases.csv"
) {
  cohort_id <- trimws(as.character(cohort_id))
  normalization_gene <- trimws(as.character(normalization_gene))

  if (!nzchar(cohort_id)) {
    stop("cohort_id must be a non-empty string.", call. = FALSE)
  }
  if (!nzchar(normalization_gene)) {
    stop("normalization_gene must be a non-empty string.", call. = FALSE)
  }
  if (tolower(cohort_id) == "pancan") {
    stop("normalize_isoform_by_gene() does not support pancan.", call. = FALSE)
  }
  if (tolower(normalization_gene) == "asitself") {
    stop("normalize_isoform_by_gene() cannot be used with normalization_gene='asitself'.", call. = FALSE)
  }

  isoform_tpm_df <- nh_read_prepared_tpm_matrix(
    cohort_id = cohort_id,
    calculated_tpm_root = calculated_tpm_root
  )
  isoform_tpm_df <- nh_filter_isoform_matrix_by_config(
    expr_df = isoform_tpm_df,
    genes_config = genes_config,
    biomart_isoform_rda = biomart_isoform_rda
  )

  gene_info <- nh_read_gene_matrix(
    cohort_id = cohort_id,
    cohort_data_roots = cohort_data_roots
  )
  gene_df <- gene_info$data
  gene_ids <- trimws(as.character(gene_df[[1L]]))
  resolved_gene <- nh_resolve_gene_symbols_against_matrix(
    requested_symbols = normalization_gene,
    available_symbols = gene_ids,
    assembly = "GRCh37",
    alias_csv = gene_symbol_alias_csv
  )
  normalization_gene_resolved <- resolved_gene$resolved_symbols[[1L]]
  if (isTRUE(resolved_gene$used_alias[[1L]])) {
    message(
      "Resolved normalization_gene alias: ",
      normalization_gene, " -> ", normalization_gene_resolved
    )
  }

  gene_match_idx <- which(toupper(gene_ids) == toupper(normalization_gene_resolved))
  if (length(gene_match_idx) < 1L) {
    stop("normalization_gene was not found in the gene expression table: ", normalization_gene, call. = FALSE)
  }
  if (length(gene_match_idx) > 1L) {
    stop("normalization_gene appears more than once in the gene expression table: ", normalization_gene, call. = FALSE)
  }

  isoform_sample_names <- names(isoform_tpm_df)[-1L]
  gene_sample_names <- names(gene_df)[-1L]
  isoform_sample_ids <- nh_canonicalize_sample_ids(isoform_sample_names)
  gene_sample_ids <- nh_canonicalize_sample_ids(gene_sample_names)

  if (anyDuplicated(isoform_sample_ids)) {
    dup_ids <- unique(isoform_sample_ids[duplicated(isoform_sample_ids)])
    stop("Duplicated isoform sample IDs after canonicalization: ", paste(head(dup_ids, 10L), collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(gene_sample_ids)) {
    dup_ids <- unique(gene_sample_ids[duplicated(gene_sample_ids)])
    stop("Duplicated gene sample IDs after canonicalization: ", paste(head(dup_ids, 10L), collapse = ", "), call. = FALSE)
  }

  missing_gene_samples <- setdiff(isoform_sample_ids, gene_sample_ids)
  if (length(missing_gene_samples) > 0L) {
    warning(
      "The gene expression table is missing sample(s) needed for normalization. ",
      "These samples will be set to NA: ",
      paste(head(missing_gene_samples, 10L), collapse = ", "),
      call. = FALSE,
      immediate. = TRUE
    )
  }

  normalized_df <- isoform_tpm_df
  for (sample_idx in seq_along(isoform_sample_ids)) {
    sample_id <- isoform_sample_ids[[sample_idx]]
    isoform_col_idx <- sample_idx + 1L
    isoform_values <- suppressWarnings(as.numeric(isoform_tpm_df[[isoform_col_idx]]))
    gene_match_col <- match(sample_id, gene_sample_ids)

    if (is.na(gene_match_col)) {
      normalized_df[[isoform_col_idx]] <- rep(NA_real_, length(isoform_values))
      next
    }

    gene_col_idx <- gene_match_col + 1L
    denominator <- suppressWarnings(as.numeric(gene_df[[gene_col_idx]][[gene_match_idx]]))

    if (!is.finite(denominator) || denominator <= 0) {
      normalized_df[[isoform_col_idx]] <- rep(NA_real_, length(isoform_values))
    } else {
      normalized_df[[isoform_col_idx]] <- isoform_values / denominator
    }
  }

  output_file <- file.path(
    normalization_root,
    toupper(cohort_id),
    "normalization",
    normalization_gene,
    paste0(toupper(cohort_id), "_normalized_isoform_matrix.csv")
  )
  nh_write_csv(normalized_df, output_file)

  attr(normalized_df, "prepared_tpm_file") <- nh_get_calculated_tpm_file(cohort_id, calculated_tpm_root = calculated_tpm_root)
  attr(normalized_df, "gene_file") <- gene_info$file
  attr(normalized_df, "output_file") <- output_file
  normalized_df
}

normalize_isoform_asitself <- function(
    cohort_id,
    cohort_data_roots = nh_default_cohort_data_roots(),
    calculated_tpm_root = file.path("intermediate", "calculated_TPM", "cohort"),
    genes_config = "config/genes.csv",
    biomart_isoform_rda = file.path("intermediate", "biomart", "NK_isoform_info.RDA"),
    normalization_root = file.path("intermediate", "cohorts"),
    gene_symbol_alias_csv = "config/gene_symbol_aliases.csv"
) {
  cohort_id <- trimws(as.character(cohort_id))

  if (!nzchar(cohort_id)) {
    stop("cohort_id must be a non-empty string.", call. = FALSE)
  }
  if (tolower(cohort_id) == "pancan") {
    stop("normalize_isoform_asitself() does not support pancan.", call. = FALSE)
  }

  isoform_tpm_df <- nh_read_prepared_tpm_matrix(
    cohort_id = cohort_id,
    calculated_tpm_root = calculated_tpm_root
  )
  isoform_tpm_df <- nh_filter_isoform_matrix_by_config(
    expr_df = isoform_tpm_df,
    genes_config = genes_config,
    biomart_isoform_rda = biomart_isoform_rda
  )
  mapped_gene_ids <- attr(isoform_tpm_df, "mapped_gene_ids")

  gene_info <- nh_read_gene_matrix(
    cohort_id = cohort_id,
    cohort_data_roots = cohort_data_roots
  )
  gene_df <- gene_info$data
  gene_ids <- trimws(as.character(gene_df[[1L]]))
  resolved_gene_map <- nh_resolve_gene_symbols_against_matrix(
    requested_symbols = mapped_gene_ids,
    available_symbols = gene_ids,
    assembly = "GRCh37",
    alias_csv = gene_symbol_alias_csv
  )
  resolved_mapped_gene_ids <- resolved_gene_map$resolved_symbols
  if (any(resolved_gene_map$used_alias)) {
    alias_pairs <- unique(paste(mapped_gene_ids[resolved_gene_map$used_alias], "->", resolved_mapped_gene_ids[resolved_gene_map$used_alias]))
    message(
      "Resolved parent gene alias(es) for 'asitself': ",
      paste(alias_pairs, collapse = ", ")
    )
  }

  gene_row_idx <- match(resolved_mapped_gene_ids, gene_ids)
  missing_genes <- unique(mapped_gene_ids[is.na(gene_row_idx)])
  if (length(missing_genes) > 0L) {
    stop(
      "The gene expression table is missing parent gene(s) needed for 'asitself': ",
      paste(head(missing_genes, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  isoform_sample_names <- names(isoform_tpm_df)[-1L]
  gene_sample_names <- names(gene_df)[-1L]
  isoform_sample_ids <- nh_canonicalize_sample_ids(isoform_sample_names)
  gene_sample_ids <- nh_canonicalize_sample_ids(gene_sample_names)

  if (anyDuplicated(isoform_sample_ids)) {
    dup_ids <- unique(isoform_sample_ids[duplicated(isoform_sample_ids)])
    stop("Duplicated isoform sample IDs after canonicalization: ", paste(head(dup_ids, 10L), collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(gene_sample_ids)) {
    dup_ids <- unique(gene_sample_ids[duplicated(gene_sample_ids)])
    stop("Duplicated gene sample IDs after canonicalization: ", paste(head(dup_ids, 10L), collapse = ", "), call. = FALSE)
  }

  missing_gene_samples <- setdiff(isoform_sample_ids, gene_sample_ids)
  if (length(missing_gene_samples) > 0L) {
    warning(
      "The gene expression table is missing sample(s) needed for 'asitself'. ",
      "These samples will be set to NA: ",
      paste(head(missing_gene_samples, 10L), collapse = ", "),
      call. = FALSE,
      immediate. = TRUE
    )
  }

  normalized_df <- isoform_tpm_df
  for (sample_idx in seq_along(isoform_sample_ids)) {
    sample_id <- isoform_sample_ids[[sample_idx]]
    isoform_col_idx <- sample_idx + 1L
    isoform_values <- suppressWarnings(as.numeric(isoform_tpm_df[[isoform_col_idx]]))
    gene_match_col <- match(sample_id, gene_sample_ids)

    if (is.na(gene_match_col)) {
      normalized_df[[isoform_col_idx]] <- rep(NA_real_, length(isoform_values))
      next
    }

    gene_col_idx <- gene_match_col + 1L
    parent_gene_values <- suppressWarnings(as.numeric(gene_df[[gene_col_idx]][gene_row_idx]))
    bad_denominator <- !is.finite(parent_gene_values) | parent_gene_values <= 0
    normalized_values <- isoform_values / parent_gene_values
    normalized_values[bad_denominator] <- NA_real_
    normalized_df[[isoform_col_idx]] <- normalized_values
  }

  output_file <- file.path(
    normalization_root,
    toupper(cohort_id),
    "normalization",
    "asitself",
    paste0(toupper(cohort_id), "_normalized_isoform_matrix.csv")
  )
  nh_write_csv(normalized_df, output_file)

  attr(normalized_df, "prepared_tpm_file") <- nh_get_calculated_tpm_file(cohort_id, calculated_tpm_root = calculated_tpm_root)
  attr(normalized_df, "gene_file") <- gene_info$file
  attr(normalized_df, "output_file") <- output_file
  normalized_df
}
