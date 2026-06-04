# scripts/run_significants_one_cohort.R

source("R/utility.R")

ensure_presto_installed <- function(
    github_repo = "immunogenomics/presto",
    cran_repo = Sys.getenv("R_CRAN_MIRROR", unset = "https://cloud.r-project.org")
) {
  if (requireNamespace("presto", quietly = TRUE)) {
    return(invisible(TRUE))
  }

  message("Package 'presto' not found. Attempting cluster-side installation from GitHub: ", github_repo)

  if (!requireNamespace("devtools", quietly = TRUE)) {
    message("Package 'devtools' not found. Installing from CRAN: ", cran_repo)
    utils::install.packages("devtools", repos = cran_repo, quiet = FALSE)
  }

  devtools::install_github(github_repo, upgrade = "never", dependencies = TRUE, quiet = FALSE)

  if (!requireNamespace("presto", quietly = TRUE)) {
    stop(
      "Automatic installation of 'presto' failed. Tried devtools::install_github('",
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

  if (requireNamespace("data.table", quietly = TRUE)) {
    df <- data.table::fread(path, data.table = FALSE, check.names = FALSE, showProgress = FALSE)
  } else {
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

  matrix_file <- resolve_matrix_file(cohort_label)
  clinical_file <- resolve_clinical_file(cohort_label)

  expr_df <- read_cluster_table(matrix_file)
  clin_df <- read_cluster_table(clinical_file)
  clin_df <- normalize_clinical_df_for_split(clin_df)

  split_res <- split_clinical_into_tertiles(clin_df)
  low_survival_clin <- split_res$low_survival_clin
  high_survival_clin <- split_res$high_survival_clin

  split_out_dir <- file.path("intermediate", "splitted_cohorts", "clin_based", cohort_label)
  dir.create(split_out_dir, recursive = TRUE, showWarnings = FALSE)
  save(low_survival_clin, file = file.path(split_out_dir, paste0(cohort_label, "_low_survival_clin.RDA")))
  save(high_survival_clin, file = file.path(split_out_dir, paste0(cohort_label, "_high_survival_clin.RDA")))

  T1_ids <- patient_ids_to_tier_samples(low_survival_clin$bcr_patient_barcode)
  T3_ids <- patient_ids_to_tier_samples(high_survival_clin$bcr_patient_barcode)
  save_tier_lists(split_out_dir, cohort_label, T1_ids, T3_ids)

  presto_inputs <- prepare_presto_inputs(expr_df, T1_ids, T3_ids)

  presto_res <- presto::wilcoxauc(
    X = presto_inputs$X,
    y = presto_inputs$y,
    groups_use = c("T1", "T3")
  )

  presto_sig <- presto_res[is.finite(presto_res$padj) & (presto_res$padj <= max_p_value), , drop = FALSE]

  presto_out_dir <- file.path("results", "presto", cohort_label)
  dir.create(presto_out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(presto_res, file.path(presto_out_dir, paste0(cohort_label, "_T1_vs_T3_presto.csv")), row.names = FALSE)
  utils::write.csv(presto_sig, file.path(presto_out_dir, paste0(cohort_label, "_T1_vs_T3_presto_padj_le_", max_p_value, ".csv")), row.names = FALSE)

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
