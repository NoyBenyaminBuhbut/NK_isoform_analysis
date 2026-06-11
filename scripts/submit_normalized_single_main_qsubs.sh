#!/bin/bash

set -euo pipefail

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo_root"

cohorts_csv="config/cohorts.csv"
genes_csv="config/genes.csv"
runner_script="scripts/run_normalized_single_main.R"
log_dir="logs/qsub_normalized_single_main"

if [[ ! -f "$cohorts_csv" ]]; then
  echo "Missing cohorts config: $cohorts_csv" >&2
  exit 1
fi
if [[ ! -f "$genes_csv" ]]; then
  echo "Missing genes config: $genes_csv" >&2
  exit 1
fi
if [[ ! -f "$runner_script" ]]; then
  echo "Missing runner script: $runner_script" >&2
  exit 1
fi
if [[ -z "${CONDA_PREFIX:-}" || ! -x "${CONDA_PREFIX}/bin/Rscript" ]]; then
  echo "Activate the conda environment first so CONDA_PREFIX/bin/Rscript exists." >&2
  exit 1
fi

mkdir -p "$log_dir"

mapfile -t cohorts < <(
  awk -F, '
    NR == 1 {
      for (j = 1; j <= NF; j++) sub(/\r$/, "", $j)
      for (i = 1; i <= NF; i++) {
        if ($i == "cohort_id") cohort_col = i
        if ($i == "next_run") next_run_col = i
      }
      next
    }
    {
      for (j = 1; j <= NF; j++) sub(/\r$/, "", $j)
    }
    cohort_col > 0 && next_run_col > 0 && toupper($next_run_col) == "TRUE" {
      print toupper($cohort_col)
    }
  ' "$cohorts_csv"
)

mapfile -t norm_genes < <(
  awk -F, '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "gene_id") gene_col = i
        if ($i == "type") type_col = i
      }
      next
    }
    gene_col > 0 && type_col > 0 && $type_col == "normalization" {
      print $gene_col
    }
  ' "$genes_csv"
)
norm_genes+=("asitself")

submitted=0
skipped_existing=0
skipped_missing_matrix=0

for cohort in "${cohorts[@]}"; do
  for gene in "${norm_genes[@]}"; do
    normalized_matrix="intermediate/cohorts/${cohort}/normalization/${gene}/${cohort}_normalized_isoform_matrix.csv"
    result_file="results/normalized_runs/${cohort}/${gene}/presto/T1_vs_T3_presto.csv"
    job_name="nsm_${cohort}_${gene}"
    log_file="${log_dir}/${job_name}.log"

    if [[ ! -f "$normalized_matrix" ]]; then
      echo "Skipping missing normalized matrix: ${normalized_matrix}"
      skipped_missing_matrix=$((skipped_missing_matrix + 1))
      continue
    fi

    if [[ -f "$result_file" ]]; then
      echo "Skipping existing normalized result: ${result_file}"
      skipped_existing=$((skipped_existing + 1))
      continue
    fi

    echo "Submitting normalized single_main: cohort=${cohort} normalization_gene=${gene}"
    qsub -b y -cwd -j y -N "$job_name" -o "$log_file" \
      "${CONDA_PREFIX}/bin/Rscript" "$runner_script" "$cohort" "$gene"
    submitted=$((submitted + 1))
  done
done

echo "Submitted ${submitted} jobs."
echo "Skipped ${skipped_existing} existing results."
echo "Skipped ${skipped_missing_matrix} missing normalized matrices."
