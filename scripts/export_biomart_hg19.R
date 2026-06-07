#!/usr/bin/env Rscript

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

repo_root <- resolve_repo_root(file.path("scripts", "export_biomart_hg19.R"))
setwd(repo_root)

source(file.path(repo_root, "R", "utility.R"))
source(file.path(repo_root, "R", "biomart_hg19_info.R"))

args <- commandArgs(trailingOnly = TRUE)
genes_config <- if (length(args) >= 1L && nzchar(args[[1L]])) args[[1L]] else "config/genes.csv"
out_dir <- if (length(args) >= 2L && nzchar(args[[2L]])) args[[2L]] else file.path("intermediate", "biomart")

result <- export_biomart_tables_hg19(
  genes_csv_path = genes_config,
  out_dir = out_dir
)

expected_files <- c(
  file.path(out_dir, "NK_genes_info.RDA"),
  file.path(out_dir, "NK_exons_info.RDA"),
  file.path(out_dir, "NK_isoform_info.RDA"),
  file.path(out_dir, "mapped_exons_to_isoforms.RDA")
)

write_manifest_file(
  folder_path = out_dir,
  manifest_entries = list(
    genes_config = normalizePath(genes_config, winslash = "/", mustWork = TRUE),
    output_files = normalizePath(expected_files, winslash = "/", mustWork = TRUE),
    supplemented_gene_symbols = result$supplemented$gene_id,
    supplemental_sources = result$supplemental_sources,
    missing_gene_symbols = result$missing$gene_id,
    transformation = "Exported BioMart GRCh37 gene, exon, isoform, and exon-to-isoform annotation tables.",
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  ),
  file_name = "manifest.txt"
)

message("Biomart annotation written to: ", out_dir)
