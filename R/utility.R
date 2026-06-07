#############################
# R/utility.R
#
# Utility functions
# Purpose:
# Enforce pipeline contracts, invariants, and pre-statistical hygiene.
# No statistics, no logging, no checkpoints, no plotting.
#############################

# ==============================================================================
# Filesystem & guards
# ==============================================================================

# Create all static (and future canonical) directories defined in stat_path.md
init_static_paths <- function(stat_path_file = "config/stat_path.md") {
  # scripts/init_static_paths.R
  # Create static directory structure from config/static_path.md
  # Idempotent: creates directories only if they do not exist
  
  spec_file <- stat_path_file
  
  if (!file.exists(spec_file)) {
    stop("Static path specification file not found: ", spec_file)
  }
  
  lines <- readLines(spec_file, warn = FALSE)
  
  # clean lines
  paths <- trimws(lines)
  paths <- paths[paths != ""]
  paths <- paths[!startsWith(paths, "#")]
  
  if (length(paths) == 0) {
    message("No static paths found in ", spec_file)
    quit(status = 0)
  }
  
  created <- character(0)
  existing <- character(0)
  
  for (p in paths) {
    if (!dir.exists(p)) {
      dir.create(p, recursive = TRUE, showWarnings = FALSE)
      created <- c(created, p)
    } else {
      existing <- c(existing, p)
    }
  }
  
  # summary
  if (length(created) > 0) {
    message("Created directories:")
    for (p in created) message("  + ", p)
  }
  
  if (length(existing) > 0) {
    message("Already existed:")
    for (p in existing) message("  = ", p)
  }
  
  message("Static path initialization completed.")
  
}

# Assert that a path exists and is a file or directory
assert_exists <- function(path, what = c("file", "dir")) {
  what <- match.arg(what)
  
  if (!file.exists(path)) {
    stop(sprintf(
      "Path does not exist (%s): %s",
      what, path
    ), call. = FALSE)
  }
  
  if (what == "file" && dir.exists(path)) {
    stop(sprintf(
      "Expected a file but found a directory: %s",
      path
    ), call. = FALSE)
  }
  
  if (what == "dir" && !dir.exists(path)) {
    stop(sprintf(
      "Expected a directory but found a file: %s",
      path
    ), call. = FALSE)
  }
  
  invisible(TRUE)
}

# Assert that a directory contains exactly one file and return its path
get_single_file_in_dir <- function(dir_path, pattern = NULL) {
  
  # Spec-aligned guard: directory must exist
  assert_exists(dir_path, what = "dir")
  
  # List immediate children (no recursion)
  entries <- list.files(
    path = dir_path,
    pattern = pattern,
    all.files = FALSE,
    full.names = TRUE,
    recursive = FALSE,
    include.dirs = TRUE,
    no.. = TRUE
  )
  
  # Separate files vs directories (subdirectories are not allowed here)
  is_dir <- dir.exists(entries)
  dirs  <- entries[is_dir]
  files <- entries[!is_dir]
  
  if (length(dirs) > 0L) {
    stop(
      paste0(
        "Expected exactly one file in directory but found subdirectory(ies): ",
        dir_path, "\n",
        paste0(" - ", basename(dirs), collapse = "\n")
      ),
      call. = FALSE
    )
  }
  
  if (length(files) == 0L) {
    msg <- if (is.null(pattern)) {
      sprintf("Expected exactly one file in directory but found 0 files: %s", dir_path)
    } else {
      sprintf("Expected exactly one file matching pattern in directory but found 0 files: %s (pattern=%s)",
              dir_path, pattern)
    }
    stop(msg, call. = FALSE)
  }
  
  if (length(files) > 1L) {
    msg <- if (is.null(pattern)) {
      paste0(
        "Expected exactly one file in directory but found ", length(files), " files: ", dir_path, "\n",
        paste0(" - ", basename(files), collapse = "\n")
      )
    } else {
      paste0(
        "Expected exactly one file matching pattern in directory but found ", length(files), " files: ",
        dir_path, " (pattern=", pattern, ")\n",
        paste0(" - ", basename(files), collapse = "\n")
      )
    }
    stop(msg, call. = FALSE)
  }
  
  # Final guard: ensure the returned path is a file (not a directory)
  assert_exists(files[[1]], what = "file")
  
  files[[1]]
}

detect_table_sep <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }

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

canonicalize_matrix_sample_ids <- function(sample_ids) {
  sample_ids <- trimws(as.character(sample_ids))
  matched <- grepl("^(([^-]+-){3}[0-9]{2}).*$", sample_ids, perl = TRUE)
  canonical <- sample_ids
  canonical[matched] <- sub("^((?:[^-]+-){3}[0-9]{2}).*$", "\\1", sample_ids[matched], perl = TRUE)

  fallback <- !matched & nchar(sample_ids) >= 15L
  canonical[fallback] <- substr(sample_ids[fallback], 1L, 15L)
  canonical
}

resolve_isoform_matrix_file <- function(
    cohort_id,
    pancan_matrix_file = Sys.getenv(
      "PANCAN_TRANSCRIPT_TABLE",
      unset = "~/tcga_transcripts/TcgaTargetGtex_rsem_isoform_tpm.gz"
    ),
    override_matrix_file = Sys.getenv("ISOFORM_MATRIX_FILE", unset = "")
) {
  cohort_id <- as.character(cohort_id)
  if (!nzchar(cohort_id)) {
    stop("cohort_id must be a non-empty string.", call. = FALSE)
  }

  if (nzchar(override_matrix_file)) {
    matrix_file <- path.expand(override_matrix_file)
    if (!file.exists(matrix_file)) {
      stop("ISOFORM_MATRIX_FILE does not exist: ", matrix_file, call. = FALSE)
    }
    return(matrix_file)
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
      file.path("GDCdata", "cohorts", tok, "isoforms"),
      file.path(path.expand("~/data/cohorts"), tok, "isoform"),
      file.path(path.expand("~/data/cohorts"), tok, "isoforms")
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

resolve_gene_matrix_file <- function(
    cohort_id,
    override_gene_matrix_file = Sys.getenv("GENE_MATRIX_FILE", unset = "")
) {
  cohort_id <- as.character(cohort_id)
  if (!nzchar(cohort_id)) {
    stop("cohort_id must be a non-empty string.", call. = FALSE)
  }

  if (nzchar(override_gene_matrix_file)) {
    gene_matrix_file <- path.expand(override_gene_matrix_file)
    if (!file.exists(gene_matrix_file)) {
      stop("GENE_MATRIX_FILE does not exist: ", gene_matrix_file, call. = FALSE)
    }
    return(gene_matrix_file)
  }

  cohort_tokens <- unique(c(cohort_id, toupper(cohort_id), tolower(cohort_id)))
  candidate_dirs <- unique(unlist(lapply(cohort_tokens, function(tok) {
    c(
      file.path("data", "cohorts", tok, "genes"),
      file.path("GDCdata", "cohorts", tok, "genes"),
      file.path(path.expand("~/data/cohorts"), tok, "genes")
    )
  })))

  for (dir_path in candidate_dirs) {
    if (dir.exists(dir_path)) {
      return(get_single_file_in_dir(dir_path))
    }
  }

  stop(
    "Could not find a gene-expression directory for cohort_id='", cohort_id, "'. Tried:\n",
    paste(candidate_dirs, collapse = "\n"),
    call. = FALSE
  )
}

resolve_biomart_isoform_info_file <- function(
    biomart_isoform_info_file = Sys.getenv(
      "BIOMART_ISOFORM_INFO_RDA",
      unset = file.path("..", "NK_exon_analysis", "intermediate", "biomart", "NK_isoform_info.RDA")
    )
) {
  resolved_path <- path.expand(biomart_isoform_info_file)
  if (!file.exists(resolved_path)) {
    stop(
      "Biomart isoform info file not found: ", resolved_path,
      "\nSet BIOMART_ISOFORM_INFO_RDA to the transcript-to-gene RDA path.",
      call. = FALSE
    )
  }
  resolved_path
}

load_biomart_isoform_info <- function(
    biomart_isoform_info_file = Sys.getenv(
      "BIOMART_ISOFORM_INFO_RDA",
      unset = file.path("..", "NK_exon_analysis", "intermediate", "biomart", "NK_isoform_info.RDA")
    )
) {
  resolved_path <- resolve_biomart_isoform_info_file(biomart_isoform_info_file)
  env <- new.env(parent = emptyenv())
  load(resolved_path, envir = env)

  if (!exists("NK_isoform_info", envir = env, inherits = FALSE)) {
    stop("Expected object 'NK_isoform_info' in: ", resolved_path, call. = FALSE)
  }

  isoform_info <- get("NK_isoform_info", envir = env, inherits = FALSE)
  if (!is.data.frame(isoform_info)) {
    stop("NK_isoform_info must be a data.frame in: ", resolved_path, call. = FALSE)
  }

  required_cols <- c("gene_id", "isoform_id")
  missing_cols <- setdiff(required_cols, names(isoform_info))
  if (length(missing_cols) > 0L) {
    stop(
      "NK_isoform_info is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  isoform_info$gene_id <- as.character(isoform_info$gene_id)
  isoform_info$isoform_id <- as.character(isoform_info$isoform_id)
  isoform_info$isoform_core <- sub("\\..*$", "", isoform_info$isoform_id)

  isoform_info <- isoform_info[
    !is.na(isoform_info$gene_id) & nzchar(isoform_info$gene_id) &
      !is.na(isoform_info$isoform_core) & nzchar(isoform_info$isoform_core),
    ,
    drop = FALSE
  ]

  dup_ids <- unique(isoform_info$isoform_core[duplicated(isoform_info$isoform_core)])
  if (length(dup_ids) > 0L) {
    stop(
      "Biomart isoform info has duplicated isoform identifiers: ",
      paste(head(dup_ids, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  isoform_info
}

read_genes_config <- function(genes_config = "config/genes.csv") {
  if (is.character(genes_config) && length(genes_config) == 1L) {
    if (!file.exists(genes_config)) {
      stop("genes_config file not found: ", genes_config, call. = FALSE)
    }
    genes_df <- read.csv(genes_config, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.data.frame(genes_config)) {
    genes_df <- genes_config
  } else {
    stop("genes_config must be either a data.frame or a single file path string.", call. = FALSE)
  }

  required_cols <- c("gene_id", "type")
  missing_cols <- setdiff(required_cols, names(genes_df))
  if (length(missing_cols) > 0L) {
    stop(
      "genes_config is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  genes_df$gene_id <- trimws(as.character(genes_df$gene_id))
  genes_df$type <- trimws(as.character(genes_df$type))
  genes_df <- genes_df[!is.na(genes_df$gene_id) & nzchar(genes_df$gene_id), , drop = FALSE]

  allowed_types <- c("experiment", "positive", "negative", "normalization")
  bad_types <- setdiff(unique(genes_df$type), allowed_types)
  if (length(bad_types) > 0L) {
    stop(
      "Invalid gene types in genes_config: ",
      paste(bad_types, collapse = ", "),
      call. = FALSE
    )
  }

  genes_df
}

normalize_gene_type_map <- function(genes_df) {
  genes_df <- read_genes_config(genes_df)
  priority_map <- c(experiment = 1L, positive = 2L, negative = 3L, normalization = 4L)
  genes_df$type_priority <- unname(priority_map[genes_df$type])

  split_rows <- split(genes_df, genes_df$gene_id, drop = TRUE)
  out_gene_id <- character(length(split_rows))
  out_type <- character(length(split_rows))
  idx <- 1L

  for (gene_name in names(split_rows)) {
    gene_rows <- split_rows[[gene_name]]
    non_norm_types <- unique(gene_rows$type[gene_rows$type != "normalization"])
    if (length(non_norm_types) > 1L) {
      stop(
        "Gene appears in multiple analysis categories in genes_config: ",
        gene_name, " -> ", paste(non_norm_types, collapse = ", "),
        call. = FALSE
      )
    }

    gene_rows <- gene_rows[order(gene_rows$type_priority), , drop = FALSE]
    out_gene_id[[idx]] <- gene_name
    out_type[[idx]] <- gene_rows$type[[1L]]
    idx <- idx + 1L
  }

  data.frame(
    gene_id = out_gene_id,
    type = out_type,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_manifest_file <- function(folder_path, manifest_entries, file_name = "manifest.txt") {
  if (missing(folder_path) || !nzchar(folder_path)) {
    stop("folder_path must be a non-empty string.", call. = FALSE)
  }
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  }
  if (!is.list(manifest_entries) || is.null(names(manifest_entries))) {
    stop("manifest_entries must be a named list.", call. = FALSE)
  }

  lines <- character(0)
  for (entry_name in names(manifest_entries)) {
    value <- manifest_entries[[entry_name]]
    if (length(value) < 1L || all(is.na(value))) {
      lines <- c(lines, paste0(entry_name, ": NA"))
      next
    }

    value_chr <- as.character(value)
    value_chr[is.na(value_chr)] <- "NA"
    if (length(value_chr) == 1L) {
      lines <- c(lines, paste0(entry_name, ": ", value_chr))
    } else {
      lines <- c(lines, paste0(entry_name, ":"))
      lines <- c(lines, paste0("  - ", value_chr))
    }
  }

  manifest_path <- file.path(folder_path, file_name)
  writeLines(lines, manifest_path, useBytes = TRUE)
  manifest_path
}


# Assert that required columns exist in a data.frame

assert_required_columns <- function(df, required_cols, df_name = "data.frame") {
  if (!is.data.frame(df)) {
    stop(sprintf("%s must be a data.frame", df_name), call. = FALSE)
  }
  
  if (!is.character(required_cols)) {
    stop("required_cols must be a character vector", call. = FALSE)
  }
  
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0L) {
    stop(
      sprintf(
        "%s is missing required columns: %s",
        df_name,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

# ==============================================================================
# Config resolution (single choke point)
# ==============================================================================

# Return a configured path/artifact value for a given cohort_id and config header
# NOTE: future canonical token expansion (<COHORT>) should be implemented here only
cohort_config_path <- function(
    cohort_id,
    dir_header,
    cohorts_csv = "config/cohorts.csv"
) {
  
  # ---- guards ----
  if (missing(cohort_id) || !nzchar(as.character(cohort_id))) {
    stop("cohort_id must be a non-empty string.", call. = FALSE)
  }
  if (missing(dir_header) || !nzchar(as.character(dir_header))) {
    stop("dir_header must be a non-empty string (a column name from config/cohorts.csv).", call. = FALSE)
  }
  if (!file.exists(cohorts_csv)) {
    stop(sprintf("Missing cohorts config file: %s", cohorts_csv), call. = FALSE)
  }
  
  # ---- read config ----
  cohorts_df <- read.csv(cohorts_csv, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!("cohort_id" %in% names(cohorts_df))) {
    stop("config/cohorts.csv must contain a 'cohort_id' column.", call. = FALSE)
  }
  if (!(dir_header %in% names(cohorts_df))) {
    stop(
      sprintf(
        "dir_header '%s' not found in config/cohorts.csv. Available columns: %s",
        dir_header, paste(names(cohorts_df), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # ---- locate cohort ----
  cohort_id <- as.character(cohort_id)
  idx <- which(as.character(cohorts_df$cohort_id) == cohort_id)
  
  if (length(idx) == 0L) {
    stop(sprintf("cohort_id '%s' not found in config/cohorts.csv.", cohort_id), call. = FALSE)
  }
  if (length(idx) > 1L) {
    stop(sprintf("cohort_id '%s' appears more than once in config/cohorts.csv.", cohort_id), call. = FALSE)
  }
  
  # ---- return configured value ----
  value <- cohorts_df[[dir_header]][idx]
  value <- as.character(value)
  
  if (!nzchar(value) || is.na(value)) {
    stop(
      sprintf(
        "Empty/NA value for cohort_id '%s' in column '%s' (config/cohorts.csv).",
        cohort_id, dir_header
      ),
      call. = FALSE
    )
  }
  
  # normalize slashes for consistency (does not touch filesystem)
  value <- gsub("\\\\", "/", value)
  
  value
}

# Save an object to the configured cohort artifact path
# Uses cohort_config_path(); does NOT create directories
save_cohort_artifact <- function(
    cohort_id,
    artifact_header,
    object,
    format = c("rda", "csv"),
    cohorts_csv = "config/cohorts.csv"
) {
  
  format <- match.arg(format)
  
  # Resolve target path
  out_path <- cohort_config_path(
    cohort_id   = cohort_id,
    dir_header = artifact_header,
    cohorts_csv = cohorts_csv
  )
  
  # Ensure parent directory exists
  out_dir <- dirname(out_path)
  if (!dir.exists(out_dir)) {
    stop(sprintf("Output directory does not exist: %s", out_dir), call. = FALSE)
  }
  
  # Write object
  if (format == "rda") {
    save(object, file = out_path)
  } else if (format == "csv") {
    if (!is.data.frame(object)) {
      stop("CSV output requires a data.frame.", call. = FALSE)
    }
    write.csv(object, out_path, row.names = FALSE)
  }
  
  invisible(out_path)
}


# ==============================================================================
# Identity & alignment
# ==============================================================================

# Compare expression sample_id to clinical patient_id
# Enforce one-sample-per-patient and compute overlap metrics
compare_expression_clinical_ids <- function(
    expr_df,
    clin_df,
    sample_id_col  = "sample_id",
    patient_id_col = "bcr_patient_barcode"
) {
}

cohort_matching_subgroup <- function(
    sub_clin_df,
    expr_df,
    cohort_id,
    patient_id_col = "bcr_patient_barcode",
    sample_id_col  = "sample_id"
) {
  if (!is.data.frame(sub_clin_df)) stop("sub_clin_df must be a data.frame", call. = FALSE)
  if (!is.data.frame(expr_df)) stop("expr_df must be a data.frame", call. = FALSE)
  
  if (missing(cohort_id) || length(cohort_id) != 1L || is.na(cohort_id) || !nzchar(cohort_id)) {
    stop("cohort_id must be a single non-empty string (e.g., 'HNSC').", call. = FALSE)
  }
  
  if (!patient_id_col %in% names(sub_clin_df)) {
    stop(sprintf("Missing patient_id_col '%s' in sub_clin_df", patient_id_col), call. = FALSE)
  }
  if (!sample_id_col %in% names(expr_df)) {
    stop(sprintf("Missing sample_id_col '%s' in expr_df", sample_id_col), call. = FALSE)
  }
  
  clinical_dir  <- file.path("data", "cohorts", cohort_id, "clinical")
  clinical_file <- get_single_file_in_dir(clinical_dir)
  
  ext <- tolower(tools::file_ext(clinical_file))
  if (ext == "csv") {
    full_clin_df <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    full_clin_df <- read.delim(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop(sprintf("Unsupported clinical file extension: %s", clinical_file), call. = FALSE)
  }
  
  if (!patient_id_col %in% names(full_clin_df)) {
    stop(sprintf("Missing patient_id_col '%s' in full_clin_df (%s)", patient_id_col, clinical_file), call. = FALSE)
  }
  
  sub_patients  <- unique(as.character(sub_clin_df[[patient_id_col]]))
  full_patients <- unique(as.character(full_clin_df[[patient_id_col]]))
  
  sub_patients  <- sub_patients[!is.na(sub_patients) & nzchar(sub_patients)]
  full_patients <- full_patients[!is.na(full_patients) & nzchar(full_patients)]
  
  if (length(sub_patients) >= length(full_patients)) {
    stop(
      paste0(
        "sub_clin_df must contain fewer patient ids than full cohort clinical.\n",
        "cohort_id: ", cohort_id, "\n",
        "n_sub_patients: ", length(sub_patients), "\n",
        "n_full_patients: ", length(full_patients)
      ),
      call. = FALSE
    )
  }
  
  # ---- sample_id and derived patient_id ----
  sample_ids <- as.character(expr_df[[sample_id_col]])
  sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
  if (length(sample_ids) == 0L) stop("expr_df has no non-empty sample_id values.", call. = FALSE)
  
  # patient_id = sample_id without the last 3 chars ("-NN")
  derived_patient_ids <- substr(sample_ids, 1, nchar(sample_ids) - 3)
  
  # keep only samples that belong to subgroup patients
  keep <- derived_patient_ids %in% sub_patients
  sub_expr_df <- expr_df[keep, , drop = FALSE]
  
  # attach patient_id (do NOT drop sample_id yet)
  sub_expr_df$patient_id <- derived_patient_ids[keep]
  
  # ---- resolve duplicates: pick smallest terminal "-NN" per patient ----
  pick_one_by_terminal_code <- function(ids) {
    ids <- as.character(ids)
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (length(ids) <= 1L) return(ids)
    
    ok <- grepl("-\\d{2}$", ids)
    if (!all(ok)) {
      bad <- ids[!ok]
      stop(
        paste0(
          "Invalid sample_id format. Expected '<prefix>-NN' (two digits at end).\n",
          paste0(" - ", head(bad, 10), collapse = "\n")
        ),
        call. = FALSE
      )
    }
    
    code <- as.integer(sub(".*-(\\d{2})$", "\\1", ids))
    ids[order(code, ids)][1]
  }
  
  if (any(duplicated(sub_expr_df$patient_id))) {
    chosen <- tapply(sub_expr_df[[sample_id_col]], sub_expr_df$patient_id, pick_one_by_terminal_code)
    sub_expr_df <- sub_expr_df[sub_expr_df[[sample_id_col]] %in% unname(chosen), , drop = FALSE]
  }
  
  # hard guarantee
  if (any(duplicated(sub_expr_df$patient_id))) {
    stop("BUG: still more than one expression sample per patient after selection.", call. = FALSE)
  }
  
  # now drop sample_id if your pipeline expects patient_id as key
  sub_expr_df[[sample_id_col]] <- NULL
  
  # overlap metrics vs subgroup clinical
  expr_patients <- unique(as.character(sub_expr_df$patient_id))
  expr_patients <- expr_patients[!is.na(expr_patients) & nzchar(expr_patients)]
  
  overlap_n   <- length(intersect(expr_patients, sub_patients))
  overlap_pct <- if (length(sub_patients) > 0L) overlap_n / length(sub_patients) else NA_real_
  
  if (!is.na(overlap_pct) && overlap_pct < 1) {
    warning(
      sprintf("Patient_id overlap is not 100%% for cohort_id=%s (overlap=%0.2f%%).",
              cohort_id, overlap_pct * 100),
      call. = FALSE
    )
  }
  
  CP_overlap <- if (is.na(overlap_pct)) "NA" else sprintf("%0.2f%%", overlap_pct * 100)
  
  list(
    sub_expr_df = sub_expr_df,
    CP_overlap  = CP_overlap
  )
}


# ==============================================================================
# Clinical hygiene (pre-statistics)
# ==============================================================================

# Recursively remove non-dead patients with low follow-up
# Ensure low-survival subgroup contains dead patients only
remove_low_survival_non_dead <- function(
    clin_df,
    patient_id_col = "bcr_patient_barcode",
    follow_up_col  = "days_to_last_follow_up",
    status_col     = "vital_status",
    dead_value     = "Dead",
    CP_patients_filtered = 0L
) {
  
  # ---- required headers (canonical + configurable) ----
  if (!is.data.frame(clin_df)) {
    stop("clin_df must be a data.frame", call. = FALSE)
  }
  req <- c(patient_id_col, follow_up_col, status_col)
  missing_cols <- setdiff(req, names(clin_df))
  if (length(missing_cols) > 0L) {
    stop(
      paste0(
        "clin_df is missing required column(s): ",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # ---- saves number of patient_id in integer patient_num ----
  patient_ids <- unique(as.character(clin_df[[patient_id_col]]))
  patient_ids <- patient_ids[nzchar(patient_ids)]
  patient_num <- as.integer(length(patient_ids))
  
  # ---- calculate average last_day_to_follow_up ----
  follow_up <- suppressWarnings(as.numeric(clin_df[[follow_up_col]]))
  avg_follow_up <- mean(follow_up, na.rm = TRUE)
  
  if (is.na(avg_follow_up)) {
    stop(
      sprintf("Cannot compute average %s (all values are NA or non-numeric).", follow_up_col),
      call. = FALSE
    )
  }
  
  # ---- remove from clin_df all non dead under average last_day_to_follow_up ----
  is_dead <- as.character(clin_df[[status_col]]) == dead_value
  keep <- is_dead | (!is_dead & (follow_up >= avg_follow_up | is.na(follow_up)))
  new_clin_df <- clin_df[keep, , drop = FALSE]
  
  # ---- save number of patient_id in integer new_patient_num ----
  new_patient_ids <- unique(as.character(new_clin_df[[patient_id_col]]))
  new_patient_ids <- new_patient_ids[nzchar(new_patient_ids)]
  new_patient_num <- as.integer(length(new_patient_ids))
  
  # ---- int CP_patients_filtered += (patient_num - new_patient_num) ----
  CP_patients_filtered <- as.integer(CP_patients_filtered + (patient_num - new_patient_num))
  
  # Guard against non-progress recursion
  if (new_patient_num == 0L) {
    stop("Filtering removed all patients; cannot form a low-survival subgroup.", call. = FALSE)
  }
  
  # ---- integer sub_group is equal to a rounded number of 1/3 times new_patient_num ----
  sub_group <- as.integer(round((1 / 3) * new_patient_num))
  if (sub_group < 1L) sub_group <- 1L
  
  # ---- low_clin_df contain sub_group number of patient_id with lowest last_day_to_follow_up ----
  # Sort by follow-up ascending; keep first sub_group rows (not unique patient-level aggregation).
  # Assumption: one row per patient; if not, upstream must standardize clinical granularity.
  ord <- order(suppressWarnings(as.numeric(new_clin_df[[follow_up_col]])), na.last = TRUE)
  low_clin_df <- new_clin_df[ord, , drop = FALSE]
  if (nrow(low_clin_df) > sub_group) {
    low_clin_df <- low_clin_df[seq_len(sub_group), , drop = FALSE]
  }
  
  # ---- checks for non dead in low_clin_df if it greater than zero ----
  low_is_dead <- as.character(low_clin_df[[status_col]]) == dead_value
  non_dead_in_low <- sum(!low_is_dead, na.rm = TRUE)
  
  # ---- if it does call remove_low_survival_non_dead(clin_df) ----
  if (non_dead_in_low > 0L) {
    
    # Prevent infinite recursion if no further filtering is possible
    if (new_patient_num == patient_num) {
      stop(
        paste0(
          "Cannot remove additional non-dead patients from low-survival subgroup; ",
          "no change after filtering step."
        ),
        call. = FALSE
      )
    }
    
    return(
      remove_low_survival_non_dead(
        clin_df = new_clin_df,
        patient_id_col = patient_id_col,
        follow_up_col = follow_up_col,
        status_col = status_col,
        dead_value = dead_value,
        CP_patients_filtered = CP_patients_filtered
      )
    )
  }
  
  # ---- else return clin_df, sub_group, CP_patients_filtered ----
  list(
    clin_df = new_clin_df,
    sub_group = sub_group,
    CP_patients_filtered = as.integer(CP_patients_filtered)
  )
}


# ==============================================================================
# Generic filtering (object-agnostic)
# ==============================================================================

# Generic allowlist-based row filter for data.frames
filter_df_by_allowlist <- function(
    df,
    allow_df,
    df_key_col,
    allow_key_col,
    keep_na = FALSE
) {
}

# ==============================================================================
# BioMart-based filtering
# ==============================================================================
# remove_chr_prefix: helper to normalize chromosome names by stripping leading "chr"
remove_chr_prefix <- function(x) {
  # Remove leading "chr" (case-insensitive). Works on vectors and preserves NAs.
  if (is.null(x)) return(x)
  x_chr <- as.character(x)
  sub("^chr", "", x_chr, ignore.case = TRUE)
}

# memory-friendly exon overlap: per-chromosome processing + interval merging
filter_df_by_biomart_exon_overlap_chunked <- function(
  df,
  biomart_exon_df,
  sample_id_col = "sample_id",
  exon_coord_col = "exon_coord",
  ignore_chr_prefix = TRUE,
  ignore_strand = FALSE,
  merge_exons = TRUE,
  merge_gap = 1L,
  dt_threads = 1L
) {
  if (!is.data.frame(df)) stop("df must be a data.frame")
  if (!is.data.frame(biomart_exon_df)) stop("biomart_exon_df must be a data.frame")
  if (!sample_id_col %in% names(df)) stop("sample_id column not found in df")
  if (!exon_coord_col %in% names(biomart_exon_df)) stop("exon_coord column not found in biomart_exon_df")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required; install it first")

  # limit data.table threads to reduce memory pressure (optional)
  old_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)
  data.table::setDTthreads(as.integer(dt_threads))

  # helper parser (reuse your project's parser if available)
  parse_coord_local <- function(coords) {
    coords <- as.character(coords)
    coords <- trimws(coords)
    coords <- sub("^chr", "", coords, ignore.case = TRUE) # normalize here to simplify
    rx <- "^([^:]+):(\\d+)-(\\d+)(?::([+-]))?$"
    m <- regmatches(coords, regexec(rx, coords, perl = TRUE))
    out <- do.call(rbind, lapply(m, function(mm) {
      if (length(mm) == 0) return(c(chr = NA_character_, start = NA_integer_, end = NA_integer_, strand = NA_character_))
      c(chr = mm[2], start = mm[3], end = mm[4], strand = ifelse(length(mm) >= 5, mm[5], NA_character_))
    }))
    out_df <- as.data.frame(out, stringsAsFactors = FALSE)
    out_df$start <- as.integer(out_df$start); out_df$end <- as.integer(out_df$end)
    out_df
  }

  # parse coordinates
  feat_cols <- setdiff(names(df), sample_id_col)
  feat_tbl <- if (exists("parse_coord_global", mode = "function")) parse_coord_global(feat_cols) else parse_coord_local(feat_cols)
  ex_tbl   <- if (exists("parse_coord_global", mode = "function")) parse_coord_global(as.character(biomart_exon_df[[exon_coord_col]])) else parse_coord_local(as.character(biomart_exon_df[[exon_coord_col]]))

  # optional extra normalization (strip chr only if asked)
  if (ignore_chr_prefix) {
    feat_tbl$chr <- sub("^chr", "", feat_tbl$chr, ignore.case = TRUE)
    ex_tbl$chr   <- sub("^chr", "", ex_tbl$chr, ignore.case = TRUE)
  }
  if (ignore_strand) {
    feat_tbl$strand <- NA_character_
    ex_tbl$strand <- NA_character_
  }

  ok_feat_idx <- which(!is.na(feat_tbl$chr) & !is.na(feat_tbl$start) & !is.na(feat_tbl$end))
  ok_ex_idx   <- which(!is.na(ex_tbl$chr) & !is.na(ex_tbl$start) & !is.na(ex_tbl$end))
  if (length(ok_feat_idx) == 0L) stop("No valid features after parsing")
  if (length(ok_ex_idx) == 0L) stop("No valid exons after parsing")

  feat_tbl_ok <- feat_tbl[ok_feat_idx, , drop = FALSE]
  feat_orig_map <- feat_cols[ok_feat_idx]     # map rows -> original feature names (subset)
  ex_tbl_ok   <- ex_tbl[ok_ex_idx, , drop = FALSE]

  # convert to data.table and key by chr,start,end when needed per-chr
  data.table::setDT(feat_tbl_ok)
  data.table::setDT(ex_tbl_ok)

  # chromosomes to iterate over = intersection of both sets (reduces wasted work)
  chrs <- intersect(unique(feat_tbl_ok$chr), unique(ex_tbl_ok$chr))
  if (length(chrs) == 0L) stop("No chromosomes in common between features and biomart exons after normalization")

  kept_features <- character(0)

  for (chr in chrs) {
    # subset per-chr
    fsub <- feat_tbl_ok[feat_tbl_ok$chr == chr]
    esub <- ex_tbl_ok[ex_tbl_ok$chr == chr]

    # merge exons for this chr if requested (coalesce overlapping/adjacent intervals)
    if (merge_exons && nrow(esub) > 0L) {
      setorder(esub, start, end)
      merged_starts <- integer(0); merged_ends <- integer(0)
      cur_s <- esub$start[1]; cur_e <- esub$end[1]
      if (nrow(esub) > 1L) {
        for (i in 2:nrow(esub)) {
          s <- esub$start[i]; e <- esub$end[i]
          if (s <= cur_e + merge_gap) {
            cur_e <- max(cur_e, e)
          } else {
            merged_starts <- c(merged_starts, cur_s); merged_ends <- c(merged_ends, cur_e)
            cur_s <- s; cur_e <- e
          }
        }
      }
      merged_starts <- c(merged_starts, cur_s); merged_ends <- c(merged_ends, cur_e)
      esub <- data.table::data.table(chr = chr, start = merged_starts, end = merged_ends)
    }

    # small sanity skip
    if (nrow(esub) == 0L || nrow(fsub) == 0L) {
      rm(fsub, esub); gc(); next
    }

    # prepare for foverlaps: need 'start' and 'end' columns and keys
    setkey(esub, start, end)
    setkey(fsub, start, end)

    # run foverlaps: as we only need feature ids, keep memory lower by subsetting minimal columns
    ov <- data.table::foverlaps(fsub[, .(start, end, .I)], esub, nomatch = 0L)
    if (nrow(ov) > 0L) {
      matched_idx <- unique(ov$.I)              # indices in fsub (these are row numbers relative to feat_tbl_ok)
      # map back to original feature names
      # fsub is a subset of feat_tbl_ok - compute global indices
      global_idx <- as.integer(rownames(fsub))   # careful: rownames on data.table are characters of the integer index
      if (is.null(global_idx) || length(global_idx) == 0L) {
        # fallback: compute using matches on start/end (safer but slower)
        matched_global <- which(feat_tbl_ok$chr == chr & feat_tbl_ok$start %in% fsub$start[matched_idx] & feat_tbl_ok$end %in% fsub$end[matched_idx])
        kept_features <- unique(c(kept_features, feat_orig_map[matched_global]))
      } else {
        kept_features <- unique(c(kept_features, feat_orig_map[as.integer(global_idx)[matched_idx]]))
      }
    }

    # cleanup per-chr objects & free memory
    rm(fsub, esub, ov); gc()
    # give simple progress message
    message("Processed chr=", chr, " ; kept so far: ", length(kept_features))
  }

  # preserve original order of features
  final_kept <- feat_cols[feat_cols %in% kept_features]
  if (length(final_kept) == 0L) stop("No overlapping feature columns found (after chunked processing)")

  # return sample_id + kept features
  res <- df[, c(sample_id_col, final_kept), drop = FALSE]

  # restore DT threads (on.exit handles it); return
  res
}
# Filter expression feature columns by overlap with BioMart exons (NO regex)
# Expects:
#   df: sample-rows expression table with a "sample_id" column
#       and feature columns named as "chr:start-end:strand"
#   biomart_exon_df: data.frame with exon_coord "chr:start-end:strand"
# Returns:
#   df subset with sample_id + only overlapping feature columns
# filter_df_by_biomart_exon_overlap: fast, robust filter of feature columns by exon overlap
filter_df_by_biomart_exon_overlap <- function(
    df,
    biomart_exon_df,
    sample_id_col = "sample_id",
    exon_coord_col = "exon_coord",
    ignore_chr_prefix = TRUE,
    ignore_strand = FALSE,
    use_foverlaps = TRUE
) {
  # Filters feature columns in `df` by overlap with exons in `biomart_exon_df`.
  # Returns df with sample_id_col + only overlapping feature columns kept.
  #
  if (!is.data.frame(df)) stop("df must be a data.frame", call. = FALSE)
  if (!is.data.frame(biomart_exon_df)) stop("biomart_exon_df must be a data.frame", call. = FALSE)
  if (!sample_id_col %in% names(df)) stop(sprintf("sample_id column '%s' not found in df", sample_id_col), call. = FALSE)
  if (!exon_coord_col %in% names(biomart_exon_df)) stop(sprintf("exon_coord column '%s' not found in biomart_exon_df", exon_coord_col), call. = FALSE)
  
  # feature columns (preserve their original order and names)
  drop_cols <- c(sample_id_col, "patient_id_tmp", "sample_type_tmp")
  feat_cols <- setdiff(names(df), drop_cols)
  if (length(feat_cols) == 0L) stop("No feature columns found in df (only sample id columns present).", call. = FALSE)
  
  # parse coordinates (use parse_coord_global if available; otherwise fallback parser below)
  parse_coord_fallback <- function(coords) {
    coords <- as.character(coords)
    coords_trim <- trimws(coords)
    coords_norm <- remove_chr_prefix(coords_trim)
    rx <- "^([^:]+):(\\d+)-(\\d+)(?::([+-]))?$"
    m <- regmatches(coords_norm, regexec(rx, coords_norm, perl = TRUE))
    out <- lapply(m, function(mm) {
      if (length(mm) == 0) return(list(chr = NA_character_, start = NA_integer_, end = NA_integer_, strand = NA_character_))
      list(chr = mm[2], start = as.integer(mm[3]), end = as.integer(mm[4]), strand = ifelse(length(mm) >= 5, mm[5], NA_character_))
    })
    dfp <- as.data.frame(do.call(rbind, lapply(out, function(x) c(x$chr, x$start, x$end, x$strand))), stringsAsFactors = FALSE)
    names(dfp) <- c("chr", "start", "end", "strand")
    dfp$start <- as.integer(dfp$start); dfp$end <- as.integer(dfp$end)
    dfp
  }
  
  if (exists("parse_coord_global", mode = "function")) {
    feat_tbl <- parse_coord_global(feat_cols)
    ex_tbl   <- parse_coord_global(as.character(biomart_exon_df[[exon_coord_col]]))
    # `parse_coord_global` may not strip "chr" — do so below if requested
  } else {
    feat_tbl <- parse_coord_fallback(feat_cols)
    ex_tbl   <- parse_coord_fallback(as.character(biomart_exon_df[[exon_coord_col]]))
  }
  
  # normalize chromosome names
  if (ignore_chr_prefix) {
    feat_tbl$chr <- remove_chr_prefix(feat_tbl$chr)
    ex_tbl$chr   <- remove_chr_prefix(ex_tbl$chr)
  }
  
  # optionally strip strand (not used for overlaps)
  if (ignore_strand) {
    feat_tbl$strand <- NA_character_
    ex_tbl$strand <- NA_character_
  }
  
  # remove invalid parsed rows
  ok_ex <- !is.na(ex_tbl$chr) & !is.na(ex_tbl$start) & !is.na(ex_tbl$end)
  ex_tbl_ok <- ex_tbl[ok_ex, , drop = FALSE]
  if (nrow(ex_tbl_ok) == 0L) stop("No valid exons after parsing exon_coord.", call. = FALSE)
  
  ok_feat <- !is.na(feat_tbl$chr) & !is.na(feat_tbl$start) & !is.na(feat_tbl$end)
  feat_tbl_ok <- feat_tbl[ok_feat, , drop = FALSE]
  if (nrow(feat_tbl_ok) == 0L) stop("No valid features after parsing feature coordinates.", call. = FALSE)
  
  # Build data.table for fast overlap if requested
  if (use_foverlaps) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("data.table is required for fast overlaps. Please install data.table or set use_foverlaps = FALSE.", call. = FALSE)
    }
    library(data.table)
    # keep mapping to original feature names
    feat_dt <- as.data.table(feat_tbl_ok)
    feat_dt[, feature := feat_cols[which(ok_feat)] ]
    bm_dt <- as.data.table(ex_tbl_ok)
    bm_dt[, bm_id := .I]
    
    # set keys for foverlaps (data.table wants 'start' and 'end')
    setkey(bm_dt, chr, start, end)
    setkey(feat_dt, chr, start, end)
    
    ov <- foverlaps(feat_dt, bm_dt, nomatch = 0L)
    if (nrow(ov) == 0L) stop("filter_df_by_biomart_exon_overlap: no overlapping feature columns found after normalization.", call. = FALSE)
    
    kept_features <- unique(ov$feature)
    # preserve original feature column order as in feat_cols
    kept_features <- feat_cols[feat_cols %in% kept_features]
  } else {
    # fallback to the original loop logic (slower). Build ex_by_chr first.
    ex_by_chr <- split(data.frame(start = ex_tbl_ok$start, end = ex_tbl_ok$end), ex_tbl_ok$chr)
    
    keep_feat <- logical(nrow(feat_tbl_ok))
    for (i in seq_len(nrow(feat_tbl_ok))) {
      chr_i <- feat_tbl_ok$chr[i]
      exons <- ex_by_chr[[chr_i]]
      if (is.null(exons)) { keep_feat[i] <- FALSE; next }
      s <- feat_tbl_ok$start[i]
      e <- feat_tbl_ok$end[i]
      keep_feat[i] <- any(s <= exons$end & e >= exons$start)
    }
    if (!any(keep_feat)) stop("filter_df_by_biomart_exon_overlap: no overlapping feature columns found after normalization.", call. = FALSE)
    kept_features <- feat_cols[which(ok_feat)][keep_feat]
  }
  
  # return sample_id + kept features
  df[, c(sample_id_col, kept_features), drop = FALSE]
}
# Filter expression feature columns by overlap with BioMart exons (NO regex)
# Expects:
#   df: sample-rows expression table with a "sample_id" column
#       and feature columns named as "chr:start-end:strand"
#   biomart_exon_df: data.frame with exon_coord "chr:start-end:strand"
# Returns:
#   df subset with sample_id + only overlapping feature columns
# Filter expression feature columns by overlap with BioMart exons (NO regex)
# + keep one TCGA sample per patient: prefer sample_type "01", else min(sample_type) with warning


filter_df_by_biomart_exon_overlap_chunked_GPT <- function(
    df,
    biomart_exon_df,
    sample_id_col = "sample_id",
    exon_coord_col = "exon_coord",
    ignore_chr_prefix = TRUE,
    ignore_strand = FALSE,
    merge_exons = TRUE,
    merge_gap = 1L,
    dt_threads = 1L
) {
  if (!is.data.frame(df)) stop("df must be a data.frame", call. = FALSE)
  if (!is.data.frame(biomart_exon_df)) stop("biomart_exon_df must be a data.frame", call. = FALSE)
  if (!sample_id_col %in% names(df)) stop("sample_id column not found in df", call. = FALSE)
  if (!exon_coord_col %in% names(biomart_exon_df)) stop("exon_coord column not found in biomart_exon_df", call. = FALSE)
  
  if (!exists("parse_coord_global", mode = "function")) {
    stop("parse_coord_global() must be loaded before calling filter_df_by_biomart_exon_overlap_chunked().", call. = FALSE)
  }
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table required; install it first", call. = FALSE)
  }
  
  remove_chr_prefix_local <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
  
  old_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)
  data.table::setDTthreads(as.integer(dt_threads))
  
  # parse coordinates using GLOBAL parser only
  feat_cols <- setdiff(names(df), sample_id_col)
  feat_tbl <- parse_coord_global(feat_cols)
  ex_tbl   <- parse_coord_global(as.character(biomart_exon_df[[exon_coord_col]]))
  
  if (ignore_chr_prefix) {
    feat_tbl$chr <- remove_chr_prefix_local(feat_tbl$chr)
    ex_tbl$chr   <- remove_chr_prefix_local(ex_tbl$chr)
  }
  if (ignore_strand) {
    feat_tbl$strand <- NA_character_
    ex_tbl$strand   <- NA_character_
  }
  
  ok_feat_idx <- which(!is.na(feat_tbl$chr) & !is.na(feat_tbl$start) & !is.na(feat_tbl$end))
  ok_ex_idx   <- which(!is.na(ex_tbl$chr)   & !is.na(ex_tbl$start)   & !is.na(ex_tbl$end))
  if (length(ok_feat_idx) == 0L) stop("No valid features after parsing", call. = FALSE)
  if (length(ok_ex_idx) == 0L)   stop("No valid exons after parsing", call. = FALSE)
  
  feat_tbl_ok   <- feat_tbl[ok_feat_idx, , drop = FALSE]
  feat_orig_map <- feat_cols[ok_feat_idx]
  ex_tbl_ok     <- ex_tbl[ok_ex_idx, , drop = FALSE]
  
  data.table::setDT(feat_tbl_ok)
  data.table::setDT(ex_tbl_ok)
  
  chrs <- intersect(unique(feat_tbl_ok$chr), unique(ex_tbl_ok$chr))
  if (length(chrs) == 0L) stop("No chromosomes in common between features and biomart exons after normalization", call. = FALSE)
  
  kept_features <- character(0)
  
  for (chr in chrs) {
    fsub <- feat_tbl_ok[feat_tbl_ok$chr == chr]
    esub <- ex_tbl_ok[ex_tbl_ok$chr == chr]
    
    if (merge_exons && nrow(esub) > 0L) {
      data.table::setorder(esub, start, end)  # <-- fixed
      
      merged_starts <- integer(0); merged_ends <- integer(0)
      cur_s <- esub$start[1]; cur_e <- esub$end[1]
      if (nrow(esub) > 1L) {
        for (i in 2:nrow(esub)) {
          s <- esub$start[i]; e <- esub$end[i]
          if (s <= cur_e + merge_gap) {
            cur_e <- max(cur_e, e)
          } else {
            merged_starts <- c(merged_starts, cur_s); merged_ends <- c(merged_ends, cur_e)
            cur_s <- s; cur_e <- e
          }
        }
      }
      merged_starts <- c(merged_starts, cur_s); merged_ends <- c(merged_ends, cur_e)
      esub <- data.table::data.table(chr = chr, start = merged_starts, end = merged_ends)
    }
    
    if (nrow(esub) == 0L || nrow(fsub) == 0L) {
      rm(fsub, esub); gc(); next
    }
    
    data.table::setkey(esub, start, end)
    data.table::setkey(fsub, start, end)
    
    ov <- data.table::foverlaps(fsub[, .(start, end, .I)], esub, nomatch = 0L)
    if (nrow(ov) > 0L) {
      matched_idx <- unique(ov$.I)
      
      # fsub keeps original row indices from feat_tbl_ok in its rownames
      global_idx <- as.integer(rownames(fsub))
      if (is.null(global_idx) || length(global_idx) == 0L) {
        matched_global <- which(feat_tbl_ok$chr == chr &
                                  feat_tbl_ok$start %in% fsub$start[matched_idx] &
                                  feat_tbl_ok$end   %in% fsub$end[matched_idx])
        kept_features <- unique(c(kept_features, feat_orig_map[matched_global]))
      } else {
        kept_features <- unique(c(kept_features, feat_orig_map[global_idx[matched_idx]]))
      }
    }
    
    rm(fsub, esub, ov); gc()
    message("Processed chr=", chr, " ; kept so far: ", length(kept_features))
  }
  
  final_kept <- feat_cols[feat_cols %in% kept_features]
  if (length(final_kept) == 0L) stop("No overlapping feature columns found (after chunked processing)", call. = FALSE)
  
  df[, c(sample_id_col, final_kept), drop = FALSE]
}

filter_df_by_biomart_exon_overlap <- function(
    df,
    biomart_exon_df,
    sample_id_col = "sample_id",
    exon_coord_col = "exon_coord",
    prefer_sample_type = "01",
    keep_smallest_if_missing = TRUE,
    warn_if_missing_preferred = TRUE
) {
  
  if (!is.data.frame(df)) stop("df must be a data.frame", call. = FALSE)
  if (!is.data.frame(biomart_exon_df)) stop("biomart_exon_df must be a data.frame", call. = FALSE)
  if (!(sample_id_col %in% names(df))) stop("Missing df column: ", sample_id_col, call. = FALSE)
  if (!(exon_coord_col %in% names(biomart_exon_df))) stop("Missing biomart_exon_df column: ", exon_coord_col, call. = FALSE)
  if (!exists("parse_coord_global", mode = "function")) {
    stop("parse_coord_global() must be loaded before calling filter_df_by_biomart_exon_overlap().", call. = FALSE)
  }
  
  # ---- sample selection: one row per patient (prefer 01) ----
  ids <- as.character(df[[sample_id_col]])
  parts <- strsplit(ids, "-", fixed = TRUE)
  
  patient_id <- vapply(parts, function(p) {
    if (length(p) >= 3) paste(p[1], p[2], p[3], sep = "-") else NA_character_
  }, character(1))
  
  sample_type <- vapply(parts, function(p) {
    if (length(p) >= 4) substr(p[4], 1, 2) else NA_character_
  }, character(1))
  
  df$patient_id_tmp <- patient_id
  df$sample_type_tmp <- sample_type
  
  # split by patient, keep preferred or fallback
  spl <- split(df, df$patient_id_tmp)
  kept_list <- lapply(spl, function(d) {
    
    st <- unique(as.character(d$sample_type_tmp))
    st <- st[!is.na(st) & nzchar(st)]
    st_sorted <- sort(st)
    
    if (length(st_sorted) == 0L) {
      # no sample_type parsed -> keep first row (avoid crash)
      return(d[1, , drop = FALSE])
    }
    
    if (prefer_sample_type %in% st_sorted) {
      d2 <- d[d$sample_type_tmp == prefer_sample_type, , drop = FALSE]
      return(d2)
    }
    
    # fallback
    if (!keep_smallest_if_missing) {
      # keep first row only
      if (warn_if_missing_preferred) {
        warning(sprintf("No sample_type %s for patient %s; keeping first available row.",
                        prefer_sample_type, unique(d$patient_id_tmp)), call. = FALSE)
      }
      return(d[1, , drop = FALSE])
    }
    
    fallback <- st_sorted[1]
    if (warn_if_missing_preferred) {
      warning(sprintf("No sample_type %s for patient %s; using %s instead.",
                      prefer_sample_type, unique(d$patient_id_tmp), fallback), call. = FALSE)
    }
    d[d$sample_type_tmp == fallback, , drop = FALSE]
  })
  
  df <- do.call(rbind, kept_list)
  rownames(df) <- NULL
  
  # ---- parse exons ----
  ex_tbl <- parse_coord_global(biomart_exon_df[[exon_coord_col]])
  ok_ex <- !is.na(ex_tbl$chr) & !is.na(ex_tbl$start) & !is.na(ex_tbl$end)
  ex_tbl <- ex_tbl[ok_ex, , drop = FALSE]
  if (nrow(ex_tbl) == 0L) stop("No valid exons after parsing exon_coord.", call. = FALSE)
  
  ex_by_chr <- split(data.frame(start = ex_tbl$start, end = ex_tbl$end), ex_tbl$chr)
  
  # ---- features are column names (excluding ids + tmp helper cols) ----
  drop_cols <- c(sample_id_col, "patient_id_tmp", "sample_type_tmp")
  feat_cols <- setdiff(names(df), drop_cols)
  if (length(feat_cols) == 0L) stop("df has no feature columns to filter.", call. = FALSE)
  
  feat_tbl <- parse_coord_global(feat_cols)
  ok_feat <- !is.na(feat_tbl$chr) & !is.na(feat_tbl$start) & !is.na(feat_tbl$end)
  
  feat_cols <- feat_cols[ok_feat]
  feat_tbl  <- feat_tbl[ok_feat, , drop = FALSE]
  if (length(feat_cols) == 0L) stop("No valid feature coords in column names.", call. = FALSE)
  
  keep_feat <- logical(length(feat_cols))
  for (i in seq_along(feat_cols)) {
    chr_i <- feat_tbl$chr[i]
    exons <- ex_by_chr[[chr_i]]
    if (is.null(exons)) { keep_feat[i] <- FALSE; next }
    s <- feat_tbl$start[i]
    e <- feat_tbl$end[i]
    keep_feat[i] <- any(s <= exons$end & e >= exons$start)
  }
  
  # return: keep sample_id + filtered features (drop tmp cols)
  out <- df[, c(sample_id_col, feat_cols[keep_feat]), drop = FALSE]
  out
}

# ==============================================================================
# files
# ==============================================================================
# Read a cohort expression TSV (or similar) and save as CSV under intermediate/
# - cohort_id: cohort name used to build input/output paths
# - input_path: optional full path to the input file (overrides default location)
# - overwrite: whether to overwrite existing output CSV (default FALSE)
# Returns the data.frame (invisibly) and writes CSV to:
#   intermediate/pre_analytics_expressions/<COHORT>/<COHORT>_expr.csv
read_and_save_cohort_expression <- function(cohort_id,
                                            input_path = NULL,
                                            overwrite = FALSE,
                                            save_root = "intermediate/pre_analytics_expressions",
                                            sep = "\t",
                                            stringsAsFactors = FALSE,
                                            check.names = FALSE) {
  if (missing(cohort_id) || !nzchar(cohort_id)) stop("cohort_id must be provided", call. = FALSE)
  
  # Build candidate input paths if input_path not provided
  if (is.null(input_path)) {
    base <- file.path("GDCdata", "cohorts", cohort_id, "expression", paste0(cohort_id, "_expr"))
    candidates <- c(base,
                    paste0(base, ".tsv"),
                    paste0(base, ".txt"),
                    paste0(base, ".tsv.gz"),
                    paste0(base, ".txt.gz"),
                    paste0(base, ".gz"))
    input_path <- NULL
    for (p in candidates) {
      if (file.exists(p)) { input_path <- p; break }
    }
    if (is.null(input_path)) {
      stop("No input file found for cohort. Tried:\n", paste(candidates, collapse = "\n"), call. = FALSE)
    }
  } else {
    if (!file.exists(input_path)) stop("Provided input_path does not exist: ", input_path, call. = FALSE)
  }
  
  message("Reading expression file: ", input_path)
  
  # Prefer data.table::fread for speed when available
  df <- NULL
  if (requireNamespace("data.table", quietly = TRUE)) {
    # use fread; returns data.frame when data.table = FALSE
    df <- data.table::fread(input_path, sep = sep, header = TRUE, stringsAsFactors = FALSE, data.table = FALSE, showProgress = FALSE)
    # fread generally preserves names; enforce check.names if requested
    if (identical(check.names, TRUE)) names(df) <- make.names(names(df))
  } else {
    # fallback to base read.delim
    df <- utils::read.delim(input_path, sep = sep, stringsAsFactors = stringsAsFactors, check.names = check.names)
  }
  
  # Ensure it's a data.frame
  if (!is.data.frame(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  # Prepare save path
  out_dir <- file.path(save_root, cohort_id)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0(cohort_id, "_expr.csv"))
  
  if (file.exists(out_file) && !overwrite) {
    message("Output file already exists and overwrite = FALSE: ", out_file, "  (skipping write).")
  } else {
    message("Writing CSV to: ", out_file)
    if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::fwrite(df, out_file)
    } else {
      utils::write.csv(df, out_file, row.names = FALSE)
    }
  }
  
  invisible(df)
}

filter_expr_cols_by_exon_endpoint_inside_chunked <- function(
    expr_df,
    biomart_exon_df,
    sample_id_col = "sample_id",
    exon_coord_col = "exon_coord",
    ignore_chr_prefix = TRUE,
    ignore_strand = TRUE,
    chunk_size = 5000L
) {
  if (!is.data.frame(expr_df)) stop("expr_df must be a data.frame", call. = FALSE)
  if (!is.data.frame(biomart_exon_df)) stop("biomart_exon_df must be a data.frame", call. = FALSE)
  if (!sample_id_col %in% names(expr_df)) stop("Missing sample_id_col in expr_df", call. = FALSE)
  if (!exon_coord_col %in% names(biomart_exon_df)) stop("Missing exon_coord_col in biomart_exon_df", call. = FALSE)
  if (!exists("parse_coord_global", mode = "function")) stop("parse_coord_global() must be loaded.", call. = FALSE)
  
  rm_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
  
  # preserve sample_id exactly
  sample_id_vec <- expr_df[[sample_id_col]]
  
  feat_cols <- setdiff(names(expr_df), sample_id_col)
  if (length(feat_cols) == 0L) stop("No feature columns to filter.", call. = FALSE)
  
  exon_vec <- as.character(biomart_exon_df[[exon_coord_col]])
  if (length(exon_vec) == 0L) stop("biomart_exon_df exon_coord is empty.", call. = FALSE)
  
  # parse exons once
  ex_tbl <- parse_coord_global(exon_vec)
  if (ignore_chr_prefix) ex_tbl$chr <- rm_chr(ex_tbl$chr)
  if (ignore_strand) ex_tbl$strand <- NA_character_
  
  ex_ok <- !is.na(ex_tbl$chr) & is.finite(ex_tbl$start) & is.finite(ex_tbl$end)
  if (!any(ex_ok)) stop("No parsable exon_coord values after parsing.", call. = FALSE)
  
  ex_chr <- as.character(ex_tbl$chr[ex_ok])
  ex_st  <- as.integer(ex_tbl$start[ex_ok])
  ex_en  <- as.integer(ex_tbl$end[ex_ok])
  
  # keep raw exon intervals (no exon merging): evaluate overlap per feature vs single exons
  ex_by_chr <- split(data.frame(start = ex_st, end = ex_en), ex_chr)
  ex_index <- lapply(ex_by_chr, function(d) {
    st <- as.integer(d$start)
    en <- as.integer(d$end)
    o <- order(st, en)
    list(start = st[o], end = en[o])
  })
  
  n <- length(feat_cols)
  chunk_size <- as.integer(chunk_size)
  if (!is.finite(chunk_size) || chunk_size < 1L) chunk_size <- 5000L
  
  kept <- logical(n)
  
  for (from in seq.int(1L, n, by = chunk_size)) {
    to <- min(n, from + chunk_size - 1L)
    cols_chunk <- feat_cols[from:to]
    
    feat_tbl <- parse_coord_global(cols_chunk)
    if (ignore_chr_prefix) feat_tbl$chr <- rm_chr(feat_tbl$chr)
    if (ignore_strand) feat_tbl$strand <- NA_character_
    
    feat_ok <- !is.na(feat_tbl$chr) & is.finite(feat_tbl$start) & is.finite(feat_tbl$end)
    keep_chunk <- logical(length(cols_chunk))
    
    if (any(feat_ok)) {
      chr_vec <- as.character(feat_tbl$chr)
      st_vec  <- as.integer(feat_tbl$start)
      en_vec  <- as.integer(feat_tbl$end)
      
      chrs <- intersect(unique(chr_vec[feat_ok]), names(ex_index))
      for (chr in chrs) {
        idx <- which(feat_ok & chr_vec == chr)
        if (length(idx) == 0L) next

        ex_start <- ex_index[[chr]]$start
        ex_end <- ex_index[[chr]]$end
        if (length(ex_start) == 0L) next

        svals <- st_vec[idx]
        evals <- en_vec[idx]

        # strict single-exon overlap: keep feature if ANY original exon overlaps it.
        # overlap condition for one pair: exon_start <= feature_end && exon_end >= feature_start
        has_overlap <- vapply(seq_along(idx), function(i) {
          any(ex_start <= evals[i] & ex_end >= svals[i])
        }, logical(1))

        keep_chunk[idx] <- has_overlap
      }
    }
    
    kept[from:to] <- keep_chunk
    
    rm(cols_chunk, feat_tbl, feat_ok, keep_chunk)
    gc()
  }
  
  kept_features <- feat_cols[kept]
  
  out <- expr_df[, c(sample_id_col, kept_features), drop = FALSE]
  
  # hard guarantee: sample_id preserved
  if (!sample_id_col %in% names(out)) stop("BUG: sample_id column was lost.", call. = FALSE)
  if (!identical(out[[sample_id_col]], sample_id_vec)) stop("BUG: sample_id values changed.", call. = FALSE)
  
  out
}


# Example:
# expr_df_filt <- filter_expr_cols_by_exoncoord_exact_or_endpoint_inside_chunked(expr_df, biomart_exon_df)

# Post-checks (optional but recommended):
# stopifnot("sample_id" %in% names(expr_df_filt))
# stopifnot(identical(expr_df_filt$sample_id, expr_df$sample_id

# Example usage:
# df <- read_and_save_cohort_expression("PAAD")
# df <- read_and_save_cohort_expression("BRCA", input_path = "GDCdata/cohorts/BRCA/expression/BRCA_expr.tsv.gz", overwrite = TRUE)

swap_expr_df <- function(expr_df,
                         sample_id_col = 1,
                         sample_id_name = "sample_id",
                         make_numeric = TRUE,
                         check_names = FALSE) {
  if (!is.data.frame(expr_df)) stop("expr_df must be a data.frame", call. = FALSE)
  
  # resolve sample_id_col to index
  if (is.character(sample_id_col)) {
    if (!sample_id_col %in% names(expr_df)) stop("sample_id_col name not found in expr_df", call. = FALSE)
    sid_idx <- which(names(expr_df) == sample_id_col)[1]
  } else if (is.numeric(sample_id_col) && sample_id_col >= 1 && sample_id_col <= ncol(expr_df)) {
    sid_idx <- as.integer(sample_id_col)
  } else {
    stop("sample_id_col must be a column name or a valid column index", call. = FALSE)
  }
  
  # extract sample ids and data portion
  sample_ids <- as.character(expr_df[[sid_idx]])
  if (any(is.na(sample_ids))) stop("sample_id column contains NA values", call. = FALSE)
  if (any(duplicated(sample_ids))) stop("sample_id values must be unique", call. = FALSE)
  
  data_df <- expr_df[ , -sid_idx, drop = FALSE]
  if (ncol(data_df) == 0L) stop("expr_df contains no feature columns to transpose", call. = FALSE)
  
  # build matrix and coerce to numeric if requested
  mat <- as.matrix(data_df)
  if (make_numeric) {
    # attempt numeric coercion; keep NAs if coercion fails for some entries
    storage.mode(mat) <- "numeric"
  }
  
  # set rownames from sample ids and transpose
  rownames(mat) <- sample_ids
  tmat <- t(mat)
  
  # construct output data.frame: leading sample_id column followed by transposed matrix
  out <- data.frame(sample_id = rownames(tmat), tmat, row.names = NULL,
                    check.names = check_names, stringsAsFactors = FALSE)
  # rename the sample id column if requested
  if (!identical(sample_id_name, "sample_id")) names(out)[1] <- sample_id_name
  
  out
}

convert_rda_to_csv <- function(rda_path, out_csv, obj_name = NULL, row_names = FALSE) {
  if (!file.exists(rda_path)) stop("RDA file not found: ", rda_path)
  e <- new.env(parent = emptyenv())
  nm <- load(rda_path, envir = e)
  if (length(nm) == 0) stop("No objects found in RDA: ", rda_path)
  if (is.null(obj_name)) {
    obj_name <- nm[[1]]  # default to first object
  } else if (!obj_name %in% nm) {
    stop("Object '", obj_name, "' not found in RDA. Available objects: ", paste(nm, collapse = ", "))
  }
  obj <- e[[obj_name]]
  # attempt to coerce matrices to data.frame; otherwise fail early
  if (!is.data.frame(obj)) {
    if (is.matrix(obj)) obj <- as.data.frame(obj)
    else stop("Object '", obj_name, "' is not a data.frame or matrix; cannot convert to CSV.")
  }
  write.csv(obj, out_csv, row.names = row_names)
  invisible(out_csv)
}

# ==============================================================================
# End of utility
# ==============================================================================
