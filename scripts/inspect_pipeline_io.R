#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) < 1L || is.na(x)) y else x
}

parse_args <- function(args) {
  opts <- list(
    cohort_id = NULL,
    normalization_gene = "asitself",
    out_file = NULL,
    run_pipeline = FALSE,
    max_p_value = 0.05
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--cohort-id")) {
      i <- i + 1L
      opts$cohort_id <- args[[i]]
    } else if (identical(arg, "--normalization-gene")) {
      i <- i + 1L
      opts$normalization_gene <- args[[i]]
    } else if (identical(arg, "--out-file")) {
      i <- i + 1L
      opts$out_file <- args[[i]]
    } else if (identical(arg, "--run-pipeline")) {
      i <- i + 1L
      value <- tolower(args[[i]])
      opts$run_pipeline <- value %in% c("1", "true", "yes", "y")
    } else if (identical(arg, "--max-p-value")) {
      i <- i + 1L
      opts$max_p_value <- as.numeric(args[[i]])
    } else if (identical(arg, "--help")) {
      cat(
        paste(
          "Usage:",
          "  Rscript scripts/inspect_pipeline_io.R --cohort-id <COHORT|pancan> [--normalization-gene <gene|asitself>] [--out-file path] [--run-pipeline true|false] [--max-p-value 0.05]",
          "",
          "Purpose:",
          "  Inspect input/output contracts and failure points for prepare_isoform_normalization.R and single_main.R.",
          "  Always writes a report, even when one or more stages fail.",
          sep = "\n"
        )
      )
      quit(status = 0L)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
    i <- i + 1L
  }

  if (is.null(opts$cohort_id) || !nzchar(opts$cohort_id)) {
    stop("--cohort-id is required.", call. = FALSE)
  }
  if (!is.finite(opts$max_p_value) || opts$max_p_value <= 0 || opts$max_p_value > 1) {
    stop("--max-p-value must be a number in (0, 1].", call. = FALSE)
  }

  opts
}

resolve_repo_root <- function(default_script_relpath) {
  script_path_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_path <- if (length(script_path_arg) > 0L) {
    sub("^--file=", "", script_path_arg[[1L]])
  } else {
    default_script_relpath
  }

  candidate_roots <- unique(c(
    tryCatch(normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE), error = function(...) NA_character_),
    tryCatch(normalizePath(Sys.getenv("SGE_O_WORKDIR", unset = ""), winslash = "/", mustWork = TRUE), error = function(...) NA_character_),
    tryCatch(normalizePath(getwd(), winslash = "/", mustWork = TRUE), error = function(...) NA_character_)
  ))
  candidate_roots <- candidate_roots[!is.na(candidate_roots) & nzchar(candidate_roots)]

  for (candidate_root in candidate_roots) {
    if (file.exists(file.path(candidate_root, "R", "utility.R"))) {
      return(candidate_root)
    }
  }

  stop(
    "Could not resolve repository root. Tried:\n",
    paste(candidate_roots, collapse = "\n"),
    call. = FALSE
  )
}

fmt_value <- function(x) {
  if (length(x) < 1L) return("NA")
  if (is.logical(x) && length(x) == 1L) return(if (isTRUE(x)) "TRUE" else "FALSE")
  if (is.numeric(x) && length(x) == 1L && is.finite(x)) return(formatC(x, digits = 8L, format = "fg", flag = "#"))
  if (length(x) == 1L && !is.na(x)) return(as.character(x))
  paste(as.character(x), collapse = ", ")
}

new_report <- function() {
  env <- new.env(parent = emptyenv())
  env$lines <- character(0)
  env
}

append_line <- function(report_env, ...) {
  report_env$lines <- c(report_env$lines, paste0(...))
}

append_named_list <- function(report_env, x, prefix = "- ") {
  for (nm in names(x)) {
    append_line(report_env, prefix, "`", nm, "`: ", fmt_value(x[[nm]]))
  }
}

safe_stage <- function(name, expr) {
  started_at <- Sys.time()
  tryCatch(
    {
      value <- force(expr)
      list(
        name = name,
        ok = TRUE,
        error = NULL,
        value = value,
        started_at = started_at,
        finished_at = Sys.time()
      )
    },
    error = function(e) {
      list(
        name = name,
        ok = FALSE,
        error = conditionMessage(e),
        value = NULL,
        started_at = started_at,
        finished_at = Sys.time()
      )
    }
  )
}

open_text_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
}

read_table_header <- function(path) {
  sep <- detect_table_sep(path)
  con <- open_text_connection(path)
  on.exit(close(con), add = TRUE)
  first_line <- readLines(con, n = 1L, warn = FALSE)
  if (length(first_line) < 1L) {
    stop("Input table is empty: ", path, call. = FALSE)
  }
  header_fields <- strsplit(first_line[[1L]], split = sep, fixed = TRUE)[[1L]]
  list(sep = sep, fields = header_fields)
}

read_table_preview <- function(path, nrows = 3L) {
  can_use_fread <- requireNamespace("data.table", quietly = TRUE) &&
    (!grepl("\\.gz$", path, ignore.case = TRUE) || requireNamespace("R.utils", quietly = TRUE))

  if (can_use_fread) {
    preview <- data.table::fread(
      path,
      nrows = nrows,
      data.table = FALSE,
      check.names = FALSE,
      showProgress = FALSE
    )
  } else {
    sep <- detect_table_sep(path)
    con <- open_text_connection(path)
    on.exit(close(con), add = TRUE)
    preview <- utils::read.table(
      con,
      sep = sep,
      header = TRUE,
      nrows = nrows,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = "",
      quote = "\""
    )
  }

  if (!is.data.frame(preview)) {
    preview <- as.data.frame(preview, stringsAsFactors = FALSE, check.names = FALSE)
  }
  preview
}

summarize_sample_ids <- function(sample_ids_raw) {
  sample_ids_raw <- as.character(sample_ids_raw)
  sample_ids_canonical <- canonicalize_matrix_sample_ids(sample_ids_raw)
  duplicate_raw <- unique(sample_ids_raw[duplicated(sample_ids_raw)])
  duplicate_canonical <- unique(sample_ids_canonical[duplicated(sample_ids_canonical)])

  list(
    sample_count = length(sample_ids_raw),
    duplicate_raw_count = length(duplicate_raw),
    duplicate_raw_examples = head(duplicate_raw, 10L),
    duplicate_canonical_count = length(duplicate_canonical),
    duplicate_canonical_examples = head(duplicate_canonical, 10L),
    raw_examples = head(sample_ids_raw, 10L),
    canonical_examples = head(sample_ids_canonical, 10L)
  )
}

profile_expression_source <- function(path, cohort_label) {
  header <- read_table_header(path)
  sample_ids_raw <- header$fields[-1L]
  sample_summary <- summarize_sample_ids(sample_ids_raw)
  preview <- read_table_preview(path, nrows = 3L)
  feature_preview <- if (ncol(preview) > 0L) head(as.character(preview[[1L]]), 10L) else character(0)

  tcga_primary_keep <- rep(FALSE, length(sample_ids_raw))
  if (tolower(cohort_label) == "pancan") {
    tcga_primary_keep <- grepl("^TCGA-[^-]+-[^-]+-[^-]+-01", sample_ids_raw, perl = TRUE)
  }

  pancan_counts <- list(
    tcga_primary_tumor_columns = sum(grepl("^TCGA-[^-]+-[^-]+-[^-]+-01", sample_ids_raw, perl = TRUE)),
    other_tcga_columns = sum(grepl("^TCGA-", sample_ids_raw) & !grepl("^TCGA-[^-]+-[^-]+-[^-]+-01", sample_ids_raw, perl = TRUE)),
    gtex_columns = sum(grepl("^GTEX-", sample_ids_raw)),
    other_columns = sum(!grepl("^(TCGA|GTEX)-", sample_ids_raw))
  )

  post_filter_sample_summary <- if (tolower(cohort_label) == "pancan") {
    summarize_sample_ids(sample_ids_raw[tcga_primary_keep])
  } else {
    NULL
  }

  list(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    separator = header$sep,
    first_column_name = header$fields[[1L]],
    preview_feature_ids = feature_preview,
    sample_summary = sample_summary,
    pancan_counts = pancan_counts,
    tcga_primary_filter_keeps = sum(tcga_primary_keep),
    post_filter_sample_summary = post_filter_sample_summary
  )
}

profile_gene_file <- function(path, target_sample_ids_canonical, needed_gene_ids) {
  stage <- safe_stage(
    paste0("profile_gene_file:", basename(path)),
    standardize_gene_expression_matrix(
      gene_file = path,
      target_sample_ids_canonical = target_sample_ids_canonical,
      needed_gene_ids = needed_gene_ids
    )
  )

  preview <- safe_stage(paste0("gene_preview:", basename(path)), read_table_preview(path, nrows = 3L))
  header <- safe_stage(paste0("gene_header:", basename(path)), read_table_header(path))

  out <- list(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    stage_ok = stage$ok,
    stage_error = stage$error
  )

  if (header$ok) {
    out$header_first_column <- header$value$fields[[1L]]
    out$header_column_count <- length(header$value$fields)
    out$header_sample_examples <- head(header$value$fields[-1L], 10L)
  }
  if (preview$ok && ncol(preview$value) > 0L) {
    out$preview_first_values <- head(as.character(preview$value[[1L]]), 10L)
  }
  if (stage$ok) {
    out$orientation <- stage$value$orientation
    out$matched_gene_rows <- nrow(stage$value$gene_matrix)
    out$sample_count <- ncol(stage$value$gene_matrix)
    out$sample_summary <- summarize_sample_ids(stage$value$sample_ids_raw)
  }

  out
}

build_split_context <- function(cohort_label) {
  if (tolower(cohort_label) == "pancan") {
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

  low_survival_clin <- split_res$low_survival_clin
  high_survival_clin <- split_res$high_survival_clin
  T1_ids <- patient_ids_to_tier_samples(low_survival_clin$bcr_patient_barcode)
  T3_ids <- patient_ids_to_tier_samples(high_survival_clin$bcr_patient_barcode)

  list(
    clinical_file = clinical_file,
    clinical_df = clin_df,
    low_survival_clin = low_survival_clin,
    high_survival_clin = high_survival_clin,
    T1_ids = T1_ids,
    T3_ids = T3_ids
  )
}

append_stage_result <- function(report_env, stage_result, detail_fun = NULL) {
  append_line(report_env, "### ", stage_result$name)
  append_line(report_env, "- `status`: ", if (stage_result$ok) "ok" else "error")
  append_line(report_env, "- `started_at`: ", format(stage_result$started_at, "%Y-%m-%d %H:%M:%S %Z"))
  append_line(report_env, "- `finished_at`: ", format(stage_result$finished_at, "%Y-%m-%d %H:%M:%S %Z"))
  if (!stage_result$ok) {
    append_line(report_env, "- `error`: ", stage_result$error)
  }
  if (is.function(detail_fun)) {
    detail_fun(stage_result)
  }
  append_line(report_env, "")
}

args <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(args)

repo_root <- resolve_repo_root(file.path("scripts", "inspect_pipeline_io.R"))
setwd(repo_root)

source(file.path(repo_root, "scripts", "single_main.R"))

cohort_label <- if (tolower(opts$cohort_id) == "pancan") "pancan" else toupper(opts$cohort_id)
normalization_label <- normalization_label_from_id(opts$normalization_gene)

out_file <- opts$out_file %||% file.path(
  repo_root,
  "reports",
  "pipeline_diagnostics",
  paste0(cohort_label, "_", normalization_label, "_pipeline_io_report.md")
)
out_file <- path.expand(out_file)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
stage_csv <- sub("\\.[^.]+$", "_stages.csv", out_file)

report_env <- new_report()
stage_rows <- list()
remember_stage <- function(stage_result) {
  stage_rows[[length(stage_rows) + 1L]] <<- data.frame(
    stage = stage_result$name,
    status = if (stage_result$ok) "ok" else "error",
    error = stage_result$error %||% "",
    started_at = format(stage_result$started_at, "%Y-%m-%d %H:%M:%S %Z"),
    finished_at = format(stage_result$finished_at, "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  stage_result
}

append_line(report_env, "# Pipeline IO Diagnostic Report")
append_line(report_env, "")
append_line(report_env, "## Scope")
append_line(report_env, "- `cohort_id`: ", cohort_label)
append_line(report_env, "- `normalization_gene`: ", normalization_label)
append_line(report_env, "- `run_pipeline`: ", if (isTRUE(opts$run_pipeline)) "TRUE" else "FALSE")
append_line(report_env, "- `generated_at`: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
append_line(report_env, "")

append_line(report_env, "## Environment")
append_named_list(
  report_env,
  list(
    repo_root = repo_root,
    cwd = getwd(),
    SGE_O_WORKDIR = Sys.getenv("SGE_O_WORKDIR", unset = ""),
    PANCAN_TRANSCRIPT_TABLE = Sys.getenv("PANCAN_TRANSCRIPT_TABLE", unset = ""),
    PANCAN_CLINICAL_FILE = Sys.getenv("PANCAN_CLINICAL_FILE", unset = ""),
    CONDA_PREFIX = Sys.getenv("CONDA_PREFIX", unset = ""),
    Rscript = Sys.which("Rscript")
  )
)
append_line(report_env, "")

append_line(report_env, "## Expected Output Paths")
append_named_list(
  report_env,
  list(
    prepare_isoform_tpm_dir = file.path(repo_root, "intermediate", "isoform_tpm", cohort_label),
    prepare_normalization_dir = file.path(repo_root, "intermediate", "normalization", normalization_label, cohort_label),
    presto_results_dir = file.path(repo_root, "results", "presto", normalization_label, cohort_label),
    split_output_dir = file.path(repo_root, "intermediate", "splitted_cohorts", "clin_based", cohort_label)
  )
)
append_line(report_env, "")

genes_stage <- remember_stage(safe_stage("read_genes_config", read_genes_config("config/genes.csv")))
append_stage_result(report_env, genes_stage, function(stage_result) {
  if (stage_result$ok) {
    genes_df <- stage_result$value
    append_line(report_env, "- `row_count`: ", nrow(genes_df))
    append_line(report_env, "- `type_counts`:")
    type_counts <- sort(table(genes_df$type), decreasing = TRUE)
    for (nm in names(type_counts)) {
      append_line(report_env, "  - `", nm, "`: ", type_counts[[nm]])
    }
  }
})

isoform_path_stage <- remember_stage(safe_stage("resolve_isoform_matrix_file", resolve_isoform_matrix_file(cohort_label)))
append_stage_result(report_env, isoform_path_stage, function(stage_result) {
  if (stage_result$ok) {
    profile <- profile_expression_source(stage_result$value, cohort_label)
    append_line(report_env, "- `path`: ", profile$path)
    append_line(report_env, "- `separator`: ", if (identical(profile$separator, "\t")) "\\t" else profile$separator)
    append_line(report_env, "- `first_column_name`: ", profile$first_column_name)
    append_line(report_env, "- `preview_feature_ids`: ", fmt_value(profile$preview_feature_ids))
    append_line(report_env, "- `sample_count`: ", profile$sample_summary$sample_count)
    append_line(report_env, "- `duplicate_raw_count`: ", profile$sample_summary$duplicate_raw_count)
    append_line(report_env, "- `duplicate_raw_examples`: ", fmt_value(profile$sample_summary$duplicate_raw_examples))
    append_line(report_env, "- `duplicate_canonical_count`: ", profile$sample_summary$duplicate_canonical_count)
    append_line(report_env, "- `duplicate_canonical_examples`: ", fmt_value(profile$sample_summary$duplicate_canonical_examples))
    append_line(report_env, "- `raw_sample_examples`: ", fmt_value(profile$sample_summary$raw_examples))
    if (tolower(cohort_label) == "pancan") {
      append_line(report_env, "- `tcga_primary_tumor_columns`: ", profile$pancan_counts$tcga_primary_tumor_columns)
      append_line(report_env, "- `other_tcga_columns`: ", profile$pancan_counts$other_tcga_columns)
      append_line(report_env, "- `gtex_columns`: ", profile$pancan_counts$gtex_columns)
      append_line(report_env, "- `other_columns`: ", profile$pancan_counts$other_columns)
      append_line(report_env, "- `post_filter_sample_count`: ", profile$post_filter_sample_summary$sample_count)
      append_line(report_env, "- `post_filter_duplicate_canonical_count`: ", profile$post_filter_sample_summary$duplicate_canonical_count)
      append_line(report_env, "- `post_filter_duplicate_canonical_examples`: ", fmt_value(profile$post_filter_sample_summary$duplicate_canonical_examples))
    }
  }
})

biomart_path_stage <- remember_stage(safe_stage("resolve_biomart_isoform_info_file", resolve_biomart_isoform_info_file()))
append_stage_result(report_env, biomart_path_stage, function(stage_result) {
  if (stage_result$ok) {
    append_line(report_env, "- `path`: ", normalizePath(stage_result$value, winslash = "/", mustWork = TRUE))
  }
})

biomart_stage <- remember_stage(safe_stage("load_biomart_isoform_info", load_biomart_isoform_info()))
append_stage_result(report_env, biomart_stage, function(stage_result) {
  if (stage_result$ok) {
    iso_df <- stage_result$value
    append_line(report_env, "- `row_count`: ", nrow(iso_df))
    append_line(report_env, "- `unique_gene_count`: ", length(unique(iso_df$gene_id)))
    append_line(report_env, "- `unique_isoform_count`: ", length(unique(iso_df$isoform_id)))
  }
})

clinical_stage <- if (tolower(cohort_label) == "pancan") {
  remember_stage(safe_stage("load_pancan_clinical_df", load_pancan_clinical_df()))
} else {
  remember_stage(safe_stage("resolve_clinical_file", resolve_clinical_file(cohort_label)))
}
append_stage_result(report_env, clinical_stage, function(stage_result) {
  if (stage_result$ok) {
    if (tolower(cohort_label) == "pancan") {
      append_line(report_env, "- `clinical_file`: ", normalizePath(stage_result$value$clinical_file, winslash = "/", mustWork = TRUE))
      append_line(report_env, "- `clinical_rows`: ", nrow(stage_result$value$clinical_df))
      append_line(report_env, "- `cohort_count`: ", length(unique(stage_result$value$clinical_df$cohort_id)))
    } else {
      append_line(report_env, "- `clinical_file`: ", normalizePath(stage_result$value, winslash = "/", mustWork = TRUE))
    }
  }
})

load_filter_stage <- remember_stage(safe_stage("load_and_filter_isoform_expression", load_and_filter_isoform_expression(cohort_label, genes_config = "config/genes.csv")))
append_stage_result(report_env, load_filter_stage, function(stage_result) {
  if (stage_result$ok) {
    value <- stage_result$value
    append_line(report_env, "- `isoform_matrix_file`: ", normalizePath(value$isoform_matrix_file, winslash = "/", mustWork = TRUE))
    append_line(report_env, "- `filtered_feature_count`: ", nrow(value$filtered_map_df))
    append_line(report_env, "- `sample_count`: ", ncol(value$isoform_expr_df) - 1L)
    append_line(report_env, "- `mapped_gene_count`: ", length(unique(value$filtered_map_df$gene_id)))
    append_line(report_env, "- `mapped_gene_examples`: ", fmt_value(head(unique(value$filtered_map_df$gene_id), 10L)))
  }
})

converted_stage <- if (load_filter_stage$ok) {
  remember_stage(safe_stage("convert_log2_tpm_df_to_tpm", convert_log2_tpm_df_to_tpm(load_filter_stage$value$isoform_expr_df)))
} else {
  remember_stage(list(name = "convert_log2_tpm_df_to_tpm", ok = FALSE, error = "Skipped because load_and_filter_isoform_expression failed.", value = NULL, started_at = Sys.time(), finished_at = Sys.time()))
}
append_stage_result(report_env, converted_stage, function(stage_result) {
  if (stage_result$ok) {
    append_line(report_env, "- `feature_count`: ", nrow(stage_result$value$expr_df))
    append_line(report_env, "- `sample_count`: ", ncol(stage_result$value$expr_df) - 1L)
    append_line(report_env, "- `value_min`: ", fmt_value(min(stage_result$value$value_matrix, na.rm = TRUE)))
    append_line(report_env, "- `value_max`: ", fmt_value(max(stage_result$value$value_matrix, na.rm = TRUE)))
  }
})

needed_gene_ids <- if (load_filter_stage$ok) {
  required_gene_ids_for_normalization(load_filter_stage$value$filtered_map_df, normalization_label)
} else {
  character(0)
}

gene_file_stage <- if (tolower(cohort_label) == "pancan") {
  remember_stage(safe_stage("discover_pancan_gene_expression_files", discover_pancan_gene_expression_files()))
} else {
  remember_stage(safe_stage("resolve_gene_matrix_file", resolve_gene_matrix_file(cohort_label)))
}
append_stage_result(report_env, gene_file_stage, function(stage_result) {
  if (stage_result$ok) {
    if (is.character(stage_result$value) && length(stage_result$value) > 1L) {
      append_line(report_env, "- `file_count`: ", length(stage_result$value))
      append_line(report_env, "- `first_files`: ", fmt_value(head(normalizePath(stage_result$value, winslash = "/", mustWork = TRUE), 10L)))
    } else {
      append_line(report_env, "- `file`: ", normalizePath(stage_result$value, winslash = "/", mustWork = TRUE))
    }
  }
})

if (gene_file_stage$ok && load_filter_stage$ok) {
  append_line(report_env, "### Gene File Profiles")
  gene_files <- if (is.character(gene_file_stage$value) && length(gene_file_stage$value) > 1L) gene_file_stage$value else c(gene_file_stage$value)
  for (gene_file in gene_files) {
    profile <- profile_gene_file(gene_file, canonicalize_matrix_sample_ids(names(load_filter_stage$value$isoform_expr_df)[-1L]), needed_gene_ids)
    append_line(report_env, "- `", basename(gene_file), "`")
    append_line(report_env, "  - `path`: ", profile$path)
    append_line(report_env, "  - `stage_ok`: ", if (isTRUE(profile$stage_ok)) "TRUE" else "FALSE")
    append_line(report_env, "  - `stage_error`: ", profile$stage_error %||% "")
    if (!is.null(profile$orientation)) append_line(report_env, "  - `orientation`: ", profile$orientation)
    if (!is.null(profile$matched_gene_rows)) append_line(report_env, "  - `matched_gene_rows`: ", profile$matched_gene_rows)
    if (!is.null(profile$sample_count)) append_line(report_env, "  - `sample_count`: ", profile$sample_count)
    if (!is.null(profile$sample_summary)) append_line(report_env, "  - `duplicate_canonical_count`: ", profile$sample_summary$duplicate_canonical_count)
  }
  append_line(report_env, "")
}

gene_matrix_stage <- if (load_filter_stage$ok) {
  remember_stage(safe_stage(
    "load_required_gene_expression_matrix",
    load_required_gene_expression_matrix(
      cohort_label = cohort_label,
      isoform_sample_ids_canonical = canonicalize_matrix_sample_ids(names(load_filter_stage$value$isoform_expr_df)[-1L]),
      needed_gene_ids = needed_gene_ids
    )
  ))
} else {
  remember_stage(list(name = "load_required_gene_expression_matrix", ok = FALSE, error = "Skipped because load_and_filter_isoform_expression failed.", value = NULL, started_at = Sys.time(), finished_at = Sys.time()))
}
append_stage_result(report_env, gene_matrix_stage, function(stage_result) {
  if (stage_result$ok) {
    append_line(report_env, "- `source_file_count`: ", length(stage_result$value$source_files))
    append_line(report_env, "- `gene_matrix_rows`: ", nrow(stage_result$value$gene_matrix))
    append_line(report_env, "- `gene_matrix_cols`: ", ncol(stage_result$value$gene_matrix))
    append_line(report_env, "- `source_files`: ", fmt_value(normalizePath(stage_result$value$source_files, winslash = "/", mustWork = TRUE)))
  }
})

aligned_stage <- if (gene_matrix_stage$ok && load_filter_stage$ok) {
  remember_stage(safe_stage(
    "align_gene_expression_matrix_to_isoforms",
    align_gene_expression_matrix_to_isoforms(
      gene_expr = gene_matrix_stage$value,
      isoform_sample_ids_raw = names(load_filter_stage$value$isoform_expr_df)[-1L]
    )
  ))
} else {
  remember_stage(list(name = "align_gene_expression_matrix_to_isoforms", ok = FALSE, error = "Skipped because gene matrix loading failed.", value = NULL, started_at = Sys.time(), finished_at = Sys.time()))
}
append_stage_result(report_env, aligned_stage, function(stage_result) {
  if (stage_result$ok) {
    append_line(report_env, "- `matched_sample_count`: ", stage_result$value$matched_sample_count)
    append_line(report_env, "- `unmatched_sample_count`: ", length(stage_result$value$unmatched_sample_ids))
    append_line(report_env, "- `unmatched_sample_examples`: ", fmt_value(head(stage_result$value$unmatched_sample_ids, 10L)))
  }
})

normalize_stage <- if (converted_stage$ok && aligned_stage$ok && load_filter_stage$ok) {
  if (normalization_label == "asitself") {
    remember_stage(safe_stage(
      "normalize_as_itself",
      normalize_as_itself(
        isoform_tpm_matrix = converted_stage$value$value_matrix,
        filtered_map_df = load_filter_stage$value$filtered_map_df,
        aligned_gene_matrix = aligned_stage$value$aligned_gene_matrix
      )
    ))
  } else {
    remember_stage(safe_stage(
      "normalize_by_single_gene",
      normalize_by_single_gene(
        isoform_tpm_matrix = converted_stage$value$value_matrix,
        aligned_gene_matrix = aligned_stage$value$aligned_gene_matrix,
        normalization_gene_label = normalization_label
      )
    ))
  }
} else {
  remember_stage(list(name = if (normalization_label == "asitself") "normalize_as_itself" else "normalize_by_single_gene", ok = FALSE, error = "Skipped because prerequisite stages failed.", value = NULL, started_at = Sys.time(), finished_at = Sys.time()))
}
append_stage_result(report_env, normalize_stage, function(stage_result) {
  if (stage_result$ok) {
    append_line(report_env, "- `normalized_rows`: ", nrow(stage_result$value$normalized_matrix))
    append_line(report_env, "- `normalized_cols`: ", ncol(stage_result$value$normalized_matrix))
    if (!is.null(stage_result$value$missing_parent_gene_ids)) {
      append_line(report_env, "- `missing_parent_gene_ids`: ", fmt_value(stage_result$value$missing_parent_gene_ids))
    }
    if (!is.null(stage_result$value$missing_denominator_samples)) {
      append_line(report_env, "- `missing_denominator_samples`: ", fmt_value(head(stage_result$value$missing_denominator_samples, 10L)))
    }
  }
})

prepare_full_stage <- remember_stage(safe_stage(
  "prepare_isoform_normalization",
  prepare_isoform_normalization(cohort_id = cohort_label, normalization_gene = normalization_label)
))
append_stage_result(report_env, prepare_full_stage, function(stage_result) {
  if (stage_result$ok) {
    append_line(report_env, "- `isoform_tpm_output_file`: ", normalizePath(stage_result$value$isoform_tpm_output_file, winslash = "/", mustWork = TRUE))
    append_line(report_env, "- `normalization_output_file`: ", normalizePath(stage_result$value$normalization_output_file, winslash = "/", mustWork = TRUE))
    append_line(report_env, "- `normalization_manifest`: ", normalizePath(stage_result$value$normalization_manifest, winslash = "/", mustWork = TRUE))
  }
})

split_stage <- remember_stage(safe_stage("build_split_context", build_split_context(cohort_label)))
append_stage_result(report_env, split_stage, function(stage_result) {
  if (stage_result$ok) {
    append_line(report_env, "- `clinical_file`: ", normalizePath(stage_result$value$clinical_file, winslash = "/", mustWork = TRUE))
    append_line(report_env, "- `clinical_rows`: ", nrow(stage_result$value$clinical_df))
    append_line(report_env, "- `T1_count`: ", length(stage_result$value$T1_ids))
    append_line(report_env, "- `T3_count`: ", length(stage_result$value$T3_ids))
  }
})

presto_input_stage <- if (prepare_full_stage$ok && split_stage$ok) {
  remember_stage(safe_stage(
    "prepare_presto_inputs",
    prepare_presto_inputs(
      expr_df = prepare_full_stage$value$normalized_expr_df,
      T1_ids = split_stage$value$T1_ids,
      T3_ids = split_stage$value$T3_ids
    )
  ))
} else {
  remember_stage(list(name = "prepare_presto_inputs", ok = FALSE, error = "Skipped because normalized matrix or split context was unavailable.", value = NULL, started_at = Sys.time(), finished_at = Sys.time()))
}
append_stage_result(report_env, presto_input_stage, function(stage_result) {
  if (stage_result$ok) {
    append_line(report_env, "- `matched_sample_count`: ", length(stage_result$value$y))
    append_line(report_env, "- `feature_count`: ", nrow(stage_result$value$X))
    append_line(report_env, "- `T1_count`: ", sum(stage_result$value$y == "T1"))
    append_line(report_env, "- `T3_count`: ", sum(stage_result$value$y == "T3"))
  }
})

presto_pkg_stage <- remember_stage(safe_stage("requireNamespace(presto)", requireNamespace("presto", quietly = TRUE)))
append_stage_result(report_env, presto_pkg_stage, function(stage_result) {
  append_line(report_env, "- `presto_available`: ", if (stage_result$ok && isTRUE(stage_result$value)) "TRUE" else "FALSE")
})

pipeline_stage <- if (isTRUE(opts$run_pipeline)) {
  remember_stage(safe_stage(
    "run_significants_pipeline",
    run_significants_pipeline(
      cohort_id = cohort_label,
      normalization_gene = normalization_label,
      max_p_value = opts$max_p_value
    )
  ))
} else {
  remember_stage(list(name = "run_significants_pipeline", ok = FALSE, error = "Skipped because --run-pipeline was not enabled.", value = NULL, started_at = Sys.time(), finished_at = Sys.time()))
}
append_stage_result(report_env, pipeline_stage, function(stage_result) {
  if (stage_result$ok) {
    append_line(report_env, "- `results_dir`: ", file.path(repo_root, "results", "presto", normalization_label, cohort_label))
  }
})

append_line(report_env, "## Current Output Existence")
expected_paths <- c(
  file.path("intermediate", "isoform_tpm", cohort_label, paste0(cohort_label, "_isoform_tpm_filtered.rds")),
  file.path("intermediate", "normalization", normalization_label, cohort_label, paste0(cohort_label, "_normalized_isoform_matrix.rds")),
  file.path("results", "presto", normalization_label, cohort_label, paste0(cohort_label, "_T1_vs_T3_presto.csv")),
  file.path("results", "presto", normalization_label, cohort_label, paste0(cohort_label, "_T1_vs_T3_presto_padj_le_", opts$max_p_value, ".csv"))
)
for (path in expected_paths) {
  append_line(report_env, "- `", path, "`: ", if (file.exists(path)) "exists" else "missing")
}
append_line(report_env, "")

stage_df <- do.call(rbind, stage_rows)
utils::write.csv(stage_df, stage_csv, row.names = FALSE)
writeLines(report_env$lines, out_file)

failed_stages <- stage_df$stage[stage_df$status == "error" & !grepl("^run_significants_pipeline$", stage_df$stage)]
if (length(failed_stages) > 0L) {
  message("Diagnostic report written to: ", out_file)
  message("Stage summary written to: ", stage_csv)
  quit(status = 1L)
}

message("Diagnostic report written to: ", out_file)
message("Stage summary written to: ", stage_csv)
