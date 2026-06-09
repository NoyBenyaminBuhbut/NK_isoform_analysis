#!/bin/bash

set -euo pipefail

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo_root"

cohorts_csv="config/cohorts.csv"
genes_csv="config/genes.csv"
runner="scripts/sge_run_normalization.sh"

if [[ ! -f "$cohorts_csv" ]]; then
  echo "Missing cohorts config: $cohorts_csv" >&2
  exit 1
fi
if [[ ! -f "$genes_csv" ]]; then
  echo "Missing genes config: $genes_csv" >&2
  exit 1
fi
if [[ ! -f "$runner" ]]; then
  echo "Missing runner script: $runner" >&2
  exit 1
fi

exclude_cohort="${1:-LUAD}"

mapfile -t cohorts < <(
  awk -F, -v exclude="$exclude_cohort" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "cohort_id") cohort_col = i
        if ($i == "next_run") next_run_col = i
      }
      next
    }
    cohort_col > 0 && next_run_col > 0 && toupper($next_run_col) == "TRUE" && $cohort_col != exclude {
      print $cohort_col
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

if [[ "${#cohorts[@]}" -lt 1 ]]; then
  echo "No eligible cohorts found in $cohorts_csv after excluding $exclude_cohort." >&2
  exit 1
fi
if [[ "${#norm_genes[@]}" -lt 1 ]]; then
  echo "No normalization genes found in $genes_csv." >&2
  exit 1
fi

for cohort in "${cohorts[@]}"; do
  for gene in "${norm_genes[@]}"; do
    job_name="norm_${cohort}_${gene}"
    echo "Submitting cohort=${cohort} normalization_gene=${gene}"
    qsub \
      -N "$job_name" \
      -v "COHORT_ID=${cohort},NORMALIZATION_GENE=${gene},CONDA_ENV_PATH=${CONDA_ENV_PATH:-},CONDA_ENV_NAME=${CONDA_ENV_NAME:-}" \
      "$runner"
  done
done
