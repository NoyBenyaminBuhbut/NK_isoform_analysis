# scripts/run_significants_one_cohort.R

source("R/utility.R")

ensure_presto_installed <- function(
    github_repo = "immunogenomics/presto",
    cran_repo = Sys.getenv("R_CRAN_MIRROR", unset = "https://cloud.r-project.org")
) {
  required_cran_packages <- c(
    "Rcpp",
    "data.table",
    "dplyr",
    "tidyr",
    "purrr",
    "tibble",
    "rlang",
    "RcppArmadillo"
  )

  if (requireNamespace("presto", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  message("Package 'presto' not found. Attempting cluster-side installation from GitHub: ", github_repo)

  if (!requireNamespace("remotes", quietly = TRUE)) {
    message("Package 'remotes' not found. Installing from CRAN: ", cran_repo)
    utils::install.packages("remotes", repos = cran_repo, quiet = FALSE)
  }

  missing_required <- required_cran_packages[
    !vapply(required_cran_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_required) > 0L) {
    message(
      "Installing required CRAN dependencies for presto: ",
      paste(missing_required, collapse = ", ")
    )
    utils::install.packages(missing_required, repos = cran_repo, quiet = FALSE)
  }

  remotes::install_github(github_repo, upgrade = "never", dependencies = FALSE, quiet = FALSE)

  if (!requireNamespace("presto", quietly = TRUE)) {
    stop(
      "Automatic installation of 'presto' failed. Tried remotes::install_github('",
      github_repo,
      "').",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

detect_table_sep <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)

  first_line <- readLines(con, n = 1L, warn = FALSE)
  if (length(first_line) < 1L || !nzchar(first_line[[1]])) {
    stop("Input table is empty: ", path, call. = FALSE)
  }

  counts <- c(
    tab = lengths(regmatches(first_line, gregexpr("\t", first_line, fixed = TRUE))),
    comma = lengths(regmatches(first_line, gregexpr(",", first_line, fixed = TRUE))),
    semicolon = lengths(regmatches(first_line, gregexpr(";", first_line, fixed = TRUE)))
  )

  c(tab = "\t", comma = ",", semicolon = ";")[[names(which.max(counts))]]
}

read_cluster_table <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }

  message("Reading table: ", path)

  can_use_fread <- requireNamespace("data.table", quietly = TRUE) &&
    (!grepl("\\.gz$", path, ignore.case = TRUE) || requireNamespace("R.utils", quietly = TRUE))

  if (can_use_fread) {
    df <- data.table::fread(path, data.table = FALSE, check.names = FALSE, showProgress = FALSE)
  } else {
    if (grepl("\\.gz$", path, ignore.case = TRUE) && requireNamespace("data.table", quietly = TRUE)) {
      message("Falling back to base R reader for gzipped table because R.utils is not installed.")
    }
    sep <- detect_table_sep(path)
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
    on.exit(close(con), add = TRUE)
    df <- utils::read.table(
      con,
      sep = sep,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = "",
      quote = "\""
    )
  }

  if (!is.data.frame(df)) {
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  }

  df
}

discover_pancan_clinical_files <- function() {
  candidate_roots <- c(
    file.path("data", "cohorts"),
    file.path("GDCdata", "cohorts"),
    path.expand("~/data/cohorts")
  )
  candidate_roots <- unique(candidate_roots[dir.exists(candidate_roots)])

  if (length(candidate_roots) < 1L) {
    stop(
      "Could not find cohort clinical roots for pancan. Tried:\n",
      paste(c(file.path("data", "cohorts"), file.path("GDCdata", "cohorts"), path.expand("~/data/cohorts")), collapse = "\n"),
      call. = FALSE
    )
  }

  files <- character(0)
  for (root_dir in candidate_roots) {
    clinical_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
    if (length(clinical_dirs) < 1L) {
      next
    }
    for (cohort_dir in clinical_dirs) {
      clinical_dir <- file.path(cohort_dir, "clinical")
      if (dir.exists(clinical_dir)) {
        files <- c(files, get_single_file_in_dir(clinical_dir, pattern = "\\.(csv|tsv|txt|gz)$"))
      }
    }
  }

  files <- unique(files[file.exists(files)])
  if (length(files) < 1L) {
    stop("No cohort clinical files were found for pancan assembly.", call. = FALSE)
  }

  files
}

infer_cohort_id_from_clinical_df <- function(clin_df) {
  candidate_cols <- c(
    "cohort_id",
    "cohort",
    "project_id",
    "project",
    "study",
    "acronym",
    "cancer_type_abbreviation"
  )
  present_cols <- intersect(candidate_cols, names(clin_df))
  if (length(present_cols) < 1L) {
    return(rep(NA_character_, nrow(clin_df)))
  }

  for (col_name in present_cols) {
    values <- trimws(as.character(clin_df[[col_name]]))
    values[!nzchar(values)] <- NA_character_
    if (all(!is.na(values))) {
      return(values)
    }
  }

  values <- trimws(as.character(clin_df[[present_cols[[1L]]]]))
  values[!nzchar(values)] <- NA_character_
  values
}

load_pancan_clinical_df <- function(
    pancan_clinical_file = Sys.getenv("PANCAN_CLINICAL_FILE", unset = ""),
    merged_output_file = file.path("intermediate", "splitted_cohorts", "clin_based", "pancan", "pancan_merged_clinical.csv")
) {
  if (nzchar(pancan_clinical_file)) {
    clinical_file <- path.expand(pancan_clinical_file)
    if (!file.exists(clinical_file)) {
      stop("PANCAN_CLINICAL_FILE does not exist: ", clinical_file, call. = FALSE)
    }

    clin_df <- read_cluster_table(clinical_file)
    clin_df <- normalize_clinical_df_for_split(clin_df)
    cohort_values <- infer_cohort_id_from_clinical_df(clin_df)
    if (anyNA(cohort_values) || !all(nzchar(cohort_values))) {
      stop(
        "Pancan clinical table must contain a cohort-identifying column such as cohort_id, cohort, project_id, or acronym.",
        call. = FALSE
      )
    }
    clin_df$cohort_id <- as.character(cohort_values)
    return(list(
      clinical_df = clin_df,
      clinical_file = clinical_file,
      auto_built = FALSE
    ))
  }

  clinical_files <- discover_pancan_clinical_files()
  message("Auto-building pancan clinical table from ", length(clinical_files), " cohort clinical files.")
  clin_list <- vector("list", length(clinical_files))

  for (i in seq_along(clinical_files)) {
    file_path <- clinical_files[[i]]
    clin_df <- read_cluster_table(file_path)
    clin_df <- normalize_clinical_df_for_split(clin_df)
    clin_df$cohort_id <- basename(dirname(dirname(file_path)))
    clin_list[[i]] <- clin_df
  }

  merged_df <- do.call(rbind, clin_list)
  merged_df$cohort_id <- as.character(merged_df$cohort_id)
  merged_df <- merged_df[!duplicated(merged_df$bcr_patient_barcode), , drop = FALSE]
  rownames(merged_df) <- NULL

  dir.create(dirname(merged_output_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(merged_df, merged_output_file, row.names = FALSE)
  message(
    "Merged pancan clinical table has ", nrow(merged_df), " rows across ",
    length(unique(merged_df$cohort_id)), " cohorts. Saved to: ", merged_output_file
  )

  list(
    clinical_df = merged_df,
    clinical_file = merged_output_file,
    auto_built = TRUE
  )
}

resolve_matrix_file <- function(
    cohort_id,
    pancan_matrix_file = Sys.getenv(
      "PANCAN_TRANSCRIPT_TABLE",
      unset = "~/tcga_transcripts/TcgaTargetGtex_rsem_isoform_tpm.gz"
    )
) {
  cohort_id <- as.character(cohort_id)
  if (!nzchar(cohort_id)) {
    stop("cohort_id must be a non-empty string.", call. = FALSE)
  }

  if (tolower(cohort_id) == "pancan") {
    matrix_file <- path.expand(pancan_matrix_file)
    if (!file.exists(matrix_file)) {
      stop(
        "Pancan transcript table not found: ", matrix_file,
        "\nSet PANCAN_TRANSCRIPT_TABLE to the raw transcript table path on the cluster.",
        call. = FALSE
      )
    }
    return(matrix_file)
  }

  cohort_tokens <- unique(c(cohort_id, toupper(cohort_id), tolower(cohort_id)))
  candidate_dirs <- unique(unlist(lapply(cohort_tokens, function(tok) {
    c(
      file.path("data", "cohorts", tok, "isoform"),
      file.path("data", "cohorts", tok, "isoforms"),
      file.path("GDCdata", "cohorts", tok, "isoform"),
      file.path("GDCdata", "cohorts", tok, "isoforms")
    )
  })))

  for (dir_path in candidate_dirs) {
    if (dir.exists(dir_path)) {
      return(get_single_file_in_dir(dir_path))
    }
  }

  stop(
    "Could not find an isoform matrix directory for cohort_id='", cohort_id, "'. Tried:\n",
    paste(candidate_dirs, collapse = "\n"),
    call. = FALSE
  )
}

resolve_clinical_file <- function(
    cohort_id,
    pancan_clinical_file = Sys.getenv("PANCAN_CLINICAL_FILE", unset = "")
) {
  cohort_id <- as.character(cohort_id)
  if (!nzchar(cohort_id)) {
    stop("cohort_id must be a non-empty string.", call. = FALSE)
  }

  if (tolower(cohort_id) == "pancan" && nzchar(pancan_clinical_file)) {
    clinical_file <- path.expand(pancan_clinical_file)
    if (!file.exists(clinical_file)) {
      stop("PANCAN_CLINICAL_FILE does not exist: ", clinical_file, call. = FALSE)
    }
    return(clinical_file)
  }

  cohort_tokens <- unique(c(cohort_id, toupper(cohort_id), tolower(cohort_id)))
  candidate_dirs <- unique(unlist(lapply(cohort_tokens, function(tok) {
    c(
      file.path("data", "cohorts", tok, "clinical"),
      file.path("GDCdata", "cohorts", tok, "clinical")
    )
  })))

  for (dir_path in candidate_dirs) {
    if (dir.exists(dir_path)) {
      return(get_single_file_in_dir(dir_path))
    }
  }

  stop(
    "Could not find a clinical directory for cohort_id='", cohort_id, "'. Tried:\n",
    paste(candidate_dirs, collapse = "\n"),
    if (tolower(cohort_id) == "pancan") {
      "\nSet PANCAN_CLINICAL_FILE if pancan clinical data lives outside the standard cohort tree."
    } else {
      ""
    },
    call. = FALSE
  )
}

normalize_clinical_df_for_split <- function(clin_df) {
  if (!is.data.frame(clin_df)) {
    stop("clin_df must be a data.frame.", call. = FALSE)
  }

  if (!("bcr_patient_barcode" %in% names(clin_df))) {
    if ("submitter_id" %in% names(clin_df)) {
      clin_df$bcr_patient_barcode <- clin_df$submitter_id
    } else if ("case_submitter_id" %in% names(clin_df)) {
      clin_df$bcr_patient_barcode <- clin_df$case_submitter_id
    } else if ("patient_id" %in% names(clin_df)) {
      clin_df$bcr_patient_barcode <- clin_df$patient_id
    }
  }

  if (!("days_to_last_follow_up" %in% names(clin_df))) {
    if ("days_to_last_known_disease_status" %in% names(clin_df)) {
      clin_df$days_to_last_follow_up <- clin_df$days_to_last_known_disease_status
    } else if ("days_to_death" %in% names(clin_df)) {
      clin_df$days_to_last_follow_up <- clin_df$days_to_death
    }
  }

  if (!("bcr_patient_barcode" %in% names(clin_df))) {
    stop("Clinical table is missing a patient ID column.", call. = FALSE)
  }
  if (!("vital_status" %in% names(clin_df))) {
    stop("Clinical table is missing vital_status.", call. = FALSE)
  }
  if (!("days_to_last_follow_up" %in% names(clin_df))) {
    stop("Clinical table is missing days_to_last_follow_up.", call. = FALSE)
  }

  clin_df$bcr_patient_barcode <- as.character(clin_df$bcr_patient_barcode)
  clin_df$vital_status <- as.character(clin_df$vital_status)
  clin_df$days_to_last_follow_up <- suppressWarnings(as.numeric(clin_df$days_to_last_follow_up))

  keep <- !is.na(clin_df$bcr_patient_barcode) &
    nzchar(clin_df$bcr_patient_barcode) &
    !is.na(clin_df$vital_status) &
    nzchar(clin_df$vital_status) &
    !is.na(clin_df$days_to_last_follow_up)

  clin_df <- clin_df[keep, , drop = FALSE]
  clin_df <- clin_df[!duplicated(clin_df$bcr_patient_barcode), , drop = FALSE]
  rownames(clin_df) <- NULL

  if (nrow(clin_df) < 2L) {
    stop("Not enough clinical rows after normalization.", call. = FALSE)
  }

  clin_df
}

split_clinical_into_tertiles <- function(clin_df) {
  clin_df <- clin_df[order(clin_df$days_to_last_follow_up, na.last = TRUE), , drop = FALSE]

  n <- nrow(clin_df)
  subgroup_n <- as.integer(floor(n / 3))
  if (subgroup_n < 1L) subgroup_n <- 1L
  if (subgroup_n > n) subgroup_n <- n

  list(
    low_survival_clin = clin_df[seq_len(subgroup_n), , drop = FALSE],
    high_survival_clin = clin_df[seq.int(from = max(1L, n - subgroup_n + 1L), to = n), , drop = FALSE]
  )
}

split_clinical_into_tertiles_by_cohort <- function(clin_df, cohort_col = "cohort_id") {
  if (!(cohort_col %in% names(clin_df))) {
    stop("Clinical table is missing cohort column '", cohort_col, "'.", call. = FALSE)
  }

  cohort_values <- trimws(as.character(clin_df[[cohort_col]]))
  if (any(is.na(cohort_values)) || any(!nzchar(cohort_values))) {
    stop("Clinical table contains missing cohort labels required for pancan stratified splitting.", call. = FALSE)
  }

  split_dfs <- split(clin_df, cohort_values, drop = TRUE)
  low_list <- vector("list", length(split_dfs))
  high_list <- vector("list", length(split_dfs))
  names(low_list) <- names(split_dfs)
  names(high_list) <- names(split_dfs)
  skipped_cohorts <- character(0)

  for (cohort_name in names(split_dfs)) {
    cohort_df <- split_dfs[[cohort_name]]
    if (nrow(cohort_df) < 3L) {
      warning(
        "Skipping cohort '", cohort_name,
        "' because it has fewer than 3 clinical rows after normalization."
      )
      skipped_cohorts <- c(skipped_cohorts, cohort_name)
      next
    }

    tertiles <- split_clinical_into_tertiles(cohort_df)
    low_list[[cohort_name]] <- tertiles$low_survival_clin
    high_list[[cohort_name]] <- tertiles$high_survival_clin
  }

  low_list <- low_list[!vapply(low_list, is.null, logical(1))]
  high_list <- high_list[!vapply(high_list, is.null, logical(1))]

  if (length(low_list) < 1L || length(high_list) < 1L) {
    stop("No pancan cohorts remained after stratified tertile filtering.", call. = FALSE)
  }

  low_survival_clin <- do.call(rbind, low_list)
  high_survival_clin <- do.call(rbind, high_list)
  rownames(low_survival_clin) <- NULL
  rownames(high_survival_clin) <- NULL

  message(
    "Built pancan T1/T3 splits from ", length(low_list), " cohorts",
    if (length(skipped_cohorts) > 0L) {
      paste0(" (skipped ", length(skipped_cohorts), " small cohorts).")
    } else {
      "."
    }
  )

  list(
    low_survival_clin = low_survival_clin,
    high_survival_clin = high_survival_clin
  )
}

patient_ids_to_tier_samples <- function(patient_ids) {
  patient_ids <- as.character(patient_ids)
  patient_ids <- patient_ids[!is.na(patient_ids) & nzchar(patient_ids)]
  unique(paste0(patient_ids, "-01"))
}

canonicalize_matrix_sample_ids <- function(sample_ids) {
  sample_ids <- trimws(as.character(sample_ids))
  matched <- grepl("^(([^-]+-){3}[0-9]{2}).*$", sample_ids, perl = TRUE)
  canonical <- sample_ids
  canonical[matched] <- sub("^((?:[^-]+-){3}[0-9]{2}).*$", "\\1", sample_ids[matched], perl = TRUE)

  fallback <- !matched & nchar(sample_ids) >= 15L
  canonical[fallback] <- substr(sample_ids[fallback], 1L, 15L)
  canonical
}

prepare_presto_inputs <- function(expr_df, T1_ids, T3_ids) {
  if (!is.data.frame(expr_df) || ncol(expr_df) < 3L) {
    stop("Expression table must contain one feature column and at least two sample columns.", call. = FALSE)
  }

  feature_ids <- as.character(expr_df[[1]])
  feature_ids[is.na(feature_ids) | !nzchar(feature_ids)] <- paste0("feature_", which(is.na(feature_ids) | !nzchar(feature_ids)))
  feature_ids <- make.unique(feature_ids)

  sample_ids_raw <- names(expr_df)[-1]
  sample_ids_canonical <- canonicalize_matrix_sample_ids(sample_ids_raw)

  if (anyDuplicated(sample_ids_canonical)) {
    dup_ids <- unique(sample_ids_canonical[duplicated(sample_ids_canonical)])
    stop(
      "Duplicated sample IDs after canonicalization: ",
      paste(head(dup_ids, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  y <- rep(NA_character_, length(sample_ids_canonical))
  y[sample_ids_canonical %in% T1_ids] <- "T1"
  y[sample_ids_canonical %in% T3_ids] <- "T3"
  keep <- !is.na(y)

  if (sum(y == "T1", na.rm = TRUE) < 1L || sum(y == "T3", na.rm = TRUE) < 1L) {
    stop("Could not match both T1 and T3 sample IDs to matrix columns.", call. = FALSE)
  }

  X <- as.matrix(expr_df[, 1L + which(keep), drop = FALSE])
  storage.mode(X) <- "numeric"
  rownames(X) <- feature_ids

  keep_features <- rowSums(is.finite(X)) > 0L
  X <- X[keep_features, , drop = FALSE]

  if (nrow(X) < 1L) {
    stop("No analyzable features remained after numeric coercion.", call. = FALSE)
  }

  list(
    X = X,
    y = y[keep],
    sample_ids_raw = sample_ids_raw[keep],
    sample_ids_canonical = sample_ids_canonical[keep]
  )
}

save_tier_lists <- function(base_dir, cohort_label, T1_ids, T3_ids) {
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

  tier_lists <- list(T1 = T1_ids, T3 = T3_ids)
  saveRDS(tier_lists, file.path(base_dir, paste0(cohort_label, "_T1_T3_lists.rds")))
  utils::write.table(
    data.frame(sample_id = T1_ids, stringsAsFactors = FALSE),
    file = file.path(base_dir, paste0(cohort_label, "_T1_samples.txt")),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  utils::write.table(
    data.frame(sample_id = T3_ids, stringsAsFactors = FALSE),
    file = file.path(base_dir, paste0(cohort_label, "_T3_samples.txt")),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}

run_significants_one_cohort <- function(cohort_id, max_p_value, max_error_rate = 0.3,
                                        cohorts_csv = "config/cohorts.csv",
                                        biomart_exon_df = NULL, biomart_gene_df = NULL,
                                        biomart_isoform_df = NULL, biomart_map_df = NULL) {
  if (missing(cohort_id) || !nzchar(as.character(cohort_id))) {
    stop("cohort_id must be provided.", call. = FALSE)
  }
  if (missing(max_p_value) || !is.numeric(max_p_value) || length(max_p_value) != 1L || !is.finite(max_p_value)) {
    stop("max_p_value must be a single finite numeric value.", call. = FALSE)
  }
  ensure_presto_installed()

  cohort_label <- if (tolower(cohort_id) == "pancan") "pancan" else toupper(as.character(cohort_id))
  message("Running isoform analysis for cohort_id=", cohort_label)

  matrix_file <- resolve_matrix_file(cohort_label)
  message("Using matrix file: ", matrix_file)
  expr_df <- read_cluster_table(matrix_file)
  message("Expression table dimensions: ", nrow(expr_df), " features x ", max(0L, ncol(expr_df) - 1L), " samples")

  if (cohort_label == "pancan") {
    pancan_clin <- load_pancan_clinical_df()
    clin_df <- pancan_clin$clinical_df
    clinical_file <- pancan_clin$clinical_file
    split_res <- split_clinical_into_tertiles_by_cohort(clin_df, cohort_col = "cohort_id")
  } else {
    clinical_file <- resolve_clinical_file(cohort_label)
    clin_df <- read_cluster_table(clinical_file)
    clin_df <- normalize_clinical_df_for_split(clin_df)
    split_res <- split_clinical_into_tertiles(clin_df)
  }
  message("Using clinical file: ", clinical_file)

  low_survival_clin <- split_res$low_survival_clin
  high_survival_clin <- split_res$high_survival_clin
  message("T1 patient count: ", nrow(low_survival_clin), "; T3 patient count: ", nrow(high_survival_clin))

  split_out_dir <- file.path("intermediate", "splitted_cohorts", "clin_based", cohort_label)
  dir.create(split_out_dir, recursive = TRUE, showWarnings = FALSE)
  save(low_survival_clin, file = file.path(split_out_dir, paste0(cohort_label, "_low_survival_clin.RDA")))
  save(high_survival_clin, file = file.path(split_out_dir, paste0(cohort_label, "_high_survival_clin.RDA")))

  T1_ids <- patient_ids_to_tier_samples(low_survival_clin$bcr_patient_barcode)
  T3_ids <- patient_ids_to_tier_samples(high_survival_clin$bcr_patient_barcode)
  message("T1 sample count: ", length(T1_ids), "; T3 sample count: ", length(T3_ids))
  save_tier_lists(split_out_dir, cohort_label, T1_ids, T3_ids)

  presto_inputs <- prepare_presto_inputs(expr_df, T1_ids, T3_ids)
  message("Matched matrix samples: ", length(presto_inputs$y), " columns for presto")

  presto_res <- presto::wilcoxauc(
    X = presto_inputs$X,
    y = presto_inputs$y,
    groups_use = c("T1", "T3")
  )
  message("Presto result rows: ", nrow(presto_res))

  presto_sig <- presto_res[is.finite(presto_res$padj) & (presto_res$padj <= max_p_value), , drop = FALSE]

  presto_out_dir <- file.path("results", "presto", cohort_label)
  dir.create(presto_out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(presto_res, file.path(presto_out_dir, paste0(cohort_label, "_T1_vs_T3_presto.csv")), row.names = FALSE)
  utils::write.csv(presto_sig, file.path(presto_out_dir, paste0(cohort_label, "_T1_vs_T3_presto_padj_le_", max_p_value, ".csv")), row.names = FALSE)
  message("Wrote presto outputs to: ", presto_out_dir)

  invisible(list(
    cohort_id = cohort_label,
    matrix_file = matrix_file,
    clinical_file = clinical_file,
    T1 = T1_ids,
    T3 = T3_ids,
    presto_results = presto_res,
    presto_significant = presto_sig
  ))
}
