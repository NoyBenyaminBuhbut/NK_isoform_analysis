#!/bin/bash

set -euo pipefail

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo_root"

cohorts_csv="config/cohorts.csv"
genes_csv="config/genes.csv"
log_dir="logs/qsub_norm"

if [[ ! -f "$cohorts_csv" ]]; then
  echo "Missing cohorts config: $cohorts_csv" >&2
  exit 1
fi
if [[ ! -f "$genes_csv" ]]; then
  echo "Missing genes config: $genes_csv" >&2
  exit 1
fi
if [[ -z "${CONDA_PREFIX:-}" || ! -x "${CONDA_PREFIX}/bin/Rscript" ]]; then
  echo "Activate the normalization conda environment first so CONDA_PREFIX/bin/Rscript exists." >&2
  exit 1
fi

mkdir -p "$log_dir"

declare -A exclude_map=()
for cohort in "$@"; do
  exclude_map["$cohort"]=1
done

mapfile -t cohorts < <(
  awk -F, '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "cohort_id") cohort_col = i
        if ($i == "next_run") next_run_col = i
      }
      next
    }
    cohort_col > 0 && next_run_col > 0 && toupper($next_run_col) == "TRUE" {
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

submitted=0
for cohort in "${cohorts[@]}"; do
  if [[ -n "${exclude_map[$cohort]:-}" ]]; then
    echo "Skipping excluded cohort=${cohort}"
    continue
  fi
  for gene in "${norm_genes[@]}"; do
    job_name="norm_${cohort}_${gene}"
    log_file="${log_dir}/${job_name}.log"
    echo "Submitting: ${job_name}"
    qsub -b y -cwd -j y -N "$job_name" -o "$log_file" \
      "${CONDA_PREFIX}/bin/Rscript" scripts/prepare_isoform_normalization.R "$cohort" "$gene"
    submitted=$((submitted + 1))
  done
done

if [[ "$submitted" -eq 0 ]]; then
  echo "No jobs submitted." >&2
  exit 1
fi

echo "Submitted ${submitted} normalization jobs."
