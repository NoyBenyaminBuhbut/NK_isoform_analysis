#!/bin/bash

set -euo pipefail

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo_root"

paths=(
  "intermediate/calculated_TPM/cohort"
  "results/normalized_runs"
)

shopt -s nullglob
for dir_path in intermediate/cohorts/*/normalization; do
  paths+=("$dir_path")
done
shopt -u nullglob

if [[ ${#paths[@]} -eq 0 ]]; then
  echo "No normalized artifact paths were found."
  exit 0
fi

echo "Normalized artifact paths:"
for path in "${paths[@]}"; do
  echo " - $path"
done

if [[ "${1:-}" != "--force" ]]; then
  echo
  echo "Dry run only. Re-run with --force to delete these paths."
  exit 0
fi

for path in "${paths[@]}"; do
  if [[ -e "$path" ]]; then
    rm -rf -- "$path"
  fi
done

echo "Deleted normalized intermediates and normalized results."
