#!/usr/bin/env Rscript

source("R/normalization_helper.R")

read_next_run_cohorts <- function(cohorts_csv = "config/cohorts.csv") {
  if (!file.exists(cohorts_csv)) {
    stop("Missing cohorts config file: ", cohorts_csv, call. = FALSE)
  }

  cohorts_df <- utils::read.csv(cohorts_csv, stringsAsFactors = FALSE, check.names = FALSE)
  assert_required_columns(cohorts_df, c("cohort_id", "next_run"), df_name = "cohorts_config")

  keep <- trimws(toupper(as.character(cohorts_df$next_run))) == "TRUE"
  cohort_ids <- trimws(as.character(cohorts_df$cohort_id[keep]))
  cohort_ids <- cohort_ids[nzchar(cohort_ids)]

  if (length(cohort_ids) < 1L) {
    stop("No cohorts are marked TRUE in config/cohorts.csv column 'next_run'.", call. = FALSE)
  }

  unique(cohort_ids)
}

run_normalization_one_cohort <- function(
    cohort_id,
    normalization_gene,
    cohort_data_roots = nh_default_cohort_data_roots(),
    calculated_tpm_root = file.path("intermediate", "calculated_TPM", "cohort"),
    genes_config = "config/genes.csv",
    biomart_isoform_rda = file.path("intermediate", "biomart", "NK_isoform_info.RDA"),
    normalization_root = file.path("intermediate", "cohorts")
) {
  prepared_tpm_file <- nh_get_calculated_tpm_file(
    cohort_id = cohort_id,
    calculated_tpm_root = calculated_tpm_root
  )
  prepared_tpm_file_legacy <- nh_get_calculated_tpm_file_legacy_csv(
    cohort_id = cohort_id,
    calculated_tpm_root = calculated_tpm_root
  )

  if (file.exists(prepared_tpm_file) || file.exists(prepared_tpm_file_legacy)) {
    existing_prepared_tpm <- if (file.exists(prepared_tpm_file)) prepared_tpm_file else prepared_tpm_file_legacy
    message("Prepared TPM artifact already exists. Skipping TPM recalculation: ", existing_prepared_tpm)
    tpm_df <- nh_read_prepared_tpm_matrix(
      cohort_id = cohort_id,
      calculated_tpm_root = calculated_tpm_root
    )
    attr(tpm_df, "output_file") <- prepared_tpm_file
  } else {
    tpm_df <- convert_log2_tpm_to_tpm(
      cohort_id = cohort_id,
      cohort_data_roots = cohort_data_roots,
      calculated_tpm_root = calculated_tpm_root
    )
  }

  normalized_df <- if (tolower(normalization_gene) == "asitself") {
    normalize_isoform_asitself(
      cohort_id = cohort_id,
      cohort_data_roots = cohort_data_roots,
      calculated_tpm_root = calculated_tpm_root,
      genes_config = genes_config,
      biomart_isoform_rda = biomart_isoform_rda,
      normalization_root = normalization_root
    )
  } else {
    normalize_isoform_by_gene(
      cohort_id = cohort_id,
      normalization_gene = normalization_gene,
      cohort_data_roots = cohort_data_roots,
      calculated_tpm_root = calculated_tpm_root,
      genes_config = genes_config,
      biomart_isoform_rda = biomart_isoform_rda,
      normalization_root = normalization_root
    )
  }

  invisible(list(
    cohort_id = toupper(cohort_id),
    normalization_gene = normalization_gene,
    calculated_tpm_file = attr(tpm_df, "output_file"),
    normalized_file = attr(normalized_df, "output_file")
  ))
}

run_normalization_request <- function(
    cohort_id,
    normalization_gene,
    cohorts_csv = "config/cohorts.csv",
    cohort_data_roots = nh_default_cohort_data_roots(),
    calculated_tpm_root = file.path("intermediate", "calculated_TPM", "cohort"),
    genes_config = "config/genes.csv",
    biomart_isoform_rda = file.path("intermediate", "biomart", "NK_isoform_info.RDA"),
    normalization_root = file.path("intermediate", "cohorts")
) {
  cohort_id <- trimws(as.character(cohort_id))
  normalization_gene <- trimws(as.character(normalization_gene))

  if (!nzchar(cohort_id)) {
    stop("cohort_id must be provided.", call. = FALSE)
  }
  if (!nzchar(normalization_gene)) {
    stop("normalization_gene must be provided.", call. = FALSE)
  }

  cohort_ids <- if (tolower(cohort_id) == "pancan") {
    read_next_run_cohorts(cohorts_csv = cohorts_csv)
  } else {
    toupper(cohort_id)
  }

  results <- vector("list", length(cohort_ids))
  failures <- character(0)

  for (idx in seq_along(cohort_ids)) {
    current_cohort <- cohort_ids[[idx]]
    message("Preparing normalization for cohort_id=", current_cohort, ", normalization_gene=", normalization_gene)

    res <- tryCatch(
      run_normalization_one_cohort(
        cohort_id = current_cohort,
        normalization_gene = normalization_gene,
        cohort_data_roots = cohort_data_roots,
        calculated_tpm_root = calculated_tpm_root,
        genes_config = genes_config,
        biomart_isoform_rda = biomart_isoform_rda,
        normalization_root = normalization_root
      ),
      error = function(err) err
    )

    if (inherits(res, "error")) {
      warn_msg <- paste0(
        "Normalization failed for cohort_id=", current_cohort,
        ", normalization_gene=", normalization_gene,
        ": ", conditionMessage(res)
      )
      warning(warn_msg, call. = FALSE, immediate. = TRUE)
      failures <- c(failures, warn_msg)
      next
    }

    results[[idx]] <- res
  }

  results <- Filter(Negate(is.null), results)

  if (tolower(cohort_id) == "pancan" && length(failures) > 0L) {
    warning(
      "Completed pancan normalization with ", length(failures),
      " cohort warning(s). Successful cohorts: ", length(results), ".",
      call. = FALSE,
      immediate. = TRUE
    )
  }

  if (length(results) < 1L) {
    stop("No cohorts completed normalization successfully.", call. = FALSE)
  }

  invisible(results)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop(
      "Usage: Rscript scripts/prepare_isoform_normalization.R <cohort_id> <normalization_gene>",
      call. = FALSE
    )
  }

  run_normalization_request(
    cohort_id = args[[1L]],
    normalization_gene = args[[2L]]
  )
}

main()
