#!/usr/bin/env Rscript

script_path_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_path_arg) < 1L) {
  stop("Could not determine script path from commandArgs().", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", script_path_arg[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: Rscript scripts/run_normalized_single_main.R <cohort_id> <normalization_gene>",
    call. = FALSE
  )
}

cohort_id <- trimws(as.character(args[[1L]]))
normalization_gene <- trimws(as.character(args[[2L]]))

if (!nzchar(cohort_id)) {
  stop("cohort_id must be a non-empty string.", call. = FALSE)
}
if (!nzchar(normalization_gene)) {
  stop("normalization_gene must be a non-empty string.", call. = FALSE)
}
if (tolower(cohort_id) == "pancan") {
  stop("Normalized single_main runs currently support cohort-level runs only, not pancan.", call. = FALSE)
}

cohort_label <- toupper(cohort_id)
normalized_matrix_file <- file.path(
  repo_root,
  "intermediate",
  "cohorts",
  cohort_label,
  "normalization",
  normalization_gene,
  paste0(cohort_label, "_normalized_isoform_matrix.RDA")
)
if (!file.exists(normalized_matrix_file)) {
  legacy_csv <- file.path(
    repo_root,
    "intermediate",
    "cohorts",
    cohort_label,
    "normalization",
    normalization_gene,
    paste0(cohort_label, "_normalized_isoform_matrix.csv")
  )
  if (file.exists(legacy_csv)) {
    normalized_matrix_file <- legacy_csv
  } else {
    stop(
      "Normalized matrix file does not exist. Tried:\n",
      normalized_matrix_file, "\n", legacy_csv,
      call. = FALSE
    )
  }
}
matrix_size <- file.info(normalized_matrix_file)$size
if (is.na(matrix_size) || matrix_size <= 0) {
  stop(
    "Normalized matrix file is empty: ",
    normalized_matrix_file,
    call. = FALSE
  )
}
normalized_matrix_file <- normalizePath(normalized_matrix_file, winslash = "/", mustWork = TRUE)

run_root <- file.path(
  repo_root,
  "results",
  "normalized_runs",
  cohort_label,
  normalization_gene
)
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)

original_wd <- getwd()
setwd(repo_root)
on.exit(setwd(original_wd), add = TRUE)

source(file.path(repo_root, "scripts", "single_main.R"))

resolve_matrix_file <- function(
    cohort_id,
    cohort_data_roots = default_cohort_data_roots(),
    pancan_matrix_file = Sys.getenv(
      "PANCAN_TRANSCRIPT_TABLE",
      unset = "~/tcga_transcripts/TcgaTargetGtex_rsem_isoform_tpm.gz"
    )
) {
  normalized_matrix_file
}

cohort_data_roots <- nh_default_cohort_data_roots()
pancan_data_roots <- c(
  cohort_data_roots,
  path.expand("~/data/cohorts")
)
max_p_value <- suppressWarnings(as.numeric(Sys.getenv("MAX_P_VALUE", unset = "0.05")))
if (!is.finite(max_p_value) || max_p_value <= 0) {
  stop("MAX_P_VALUE must be a positive finite numeric value.", call. = FALSE)
}

setwd(run_root)

message("Running normalized single_main for cohort_id=", cohort_label, ", normalization_gene=", normalization_gene)
message("Using normalized matrix: ", normalized_matrix_file)
message("Writing outputs under: ", run_root)

res <- run_significants_one_cohort(
  cohort_id = cohort_label,
  max_p_value = max_p_value,
  cohort_data_roots = cohort_data_roots,
  pancan_data_roots = pancan_data_roots
)

presto_final_dir <- file.path(run_root, "presto")
split_final_dir <- file.path(run_root, "split")
meta_final_dir <- file.path(run_root, "meta")
dir.create(presto_final_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(split_final_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_final_dir, recursive = TRUE, showWarnings = FALSE)

presto_tmp_dir <- file.path(run_root, "results", "presto", cohort_label)
split_tmp_dir <- file.path(run_root, "intermediate", "splitted_cohorts", "clin_based", cohort_label)

move_artifact <- function(from_path, to_path) {
  if (!file.exists(from_path)) {
    stop("Expected output artifact not found: ", from_path, call. = FALSE)
  }
  if (file.exists(to_path)) {
    unlink(to_path, recursive = TRUE, force = TRUE)
  }
  ok <- file.rename(from_path, to_path)
  if (!ok) {
    ok <- file.copy(from_path, to_path, overwrite = TRUE)
    if (!ok) {
      stop("Failed to move artifact from ", from_path, " to ", to_path, call. = FALSE)
    }
    unlink(from_path, recursive = TRUE, force = TRUE)
  }
}

padj_suffix <- as.character(max_p_value)
move_artifact(
  file.path(presto_tmp_dir, paste0(cohort_label, "_T1_vs_T3_presto.csv")),
  file.path(presto_final_dir, "T1_vs_T3_presto.csv")
)
move_artifact(
  file.path(presto_tmp_dir, paste0(cohort_label, "_T1_vs_T3_presto_padj_le_", padj_suffix, ".csv")),
  file.path(presto_final_dir, paste0("T1_vs_T3_presto_padj_le_", padj_suffix, ".csv"))
)
move_artifact(
  file.path(split_tmp_dir, paste0(cohort_label, "_low_survival_clin.RDA")),
  file.path(split_final_dir, "low_survival_clin.RDA")
)
move_artifact(
  file.path(split_tmp_dir, paste0(cohort_label, "_high_survival_clin.RDA")),
  file.path(split_final_dir, "high_survival_clin.RDA")
)
move_artifact(
  file.path(split_tmp_dir, paste0(cohort_label, "_T1_T3_lists.rds")),
  file.path(split_final_dir, "T1_T3_lists.rds")
)
move_artifact(
  file.path(split_tmp_dir, paste0(cohort_label, "_T1_samples.txt")),
  file.path(split_final_dir, "T1_samples.txt")
)
move_artifact(
  file.path(split_tmp_dir, paste0(cohort_label, "_T3_samples.txt")),
  file.path(split_final_dir, "T3_samples.txt")
)

run_info_lines <- c(
  paste0("cohort_id=", cohort_label),
  paste0("normalization_gene=", normalization_gene),
  paste0("normalized_matrix_file=", normalized_matrix_file),
  paste0("clinical_file=", normalizePath(res$clinical_file, winslash = "/", mustWork = TRUE)),
  paste0("max_p_value=", padj_suffix),
  paste0("run_root=", normalizePath(run_root, winslash = "/", mustWork = TRUE)),
  paste0("timestamp=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)
writeLines(run_info_lines, con = file.path(meta_final_dir, "run_info.txt"))

unlink(file.path(run_root, "results"), recursive = TRUE, force = TRUE)
unlink(file.path(run_root, "intermediate"), recursive = TRUE, force = TRUE)

invisible(res)
