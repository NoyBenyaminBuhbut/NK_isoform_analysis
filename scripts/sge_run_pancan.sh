#!/bin/bash
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N nk_iso_pancan
#$ -j y

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "$REPO_ROOT"

export PATH="$HOME/.local/bin:$PATH"
export PANCAN_TRANSCRIPT_TABLE="${PANCAN_TRANSCRIPT_TABLE:-$HOME/tcga_transcripts/TcgaTargetGtex_rsem_isoform_tpm.gz}"
export PANCAN_CLINICAL_FILE="${PANCAN_CLINICAL_FILE:-$HOME/tcga_transcripts/pancan_clinical.csv}"

Rscript -e 'source("scripts/single_main.R"); run_significants_one_cohort("pancan", max_p_value=0.05)'
