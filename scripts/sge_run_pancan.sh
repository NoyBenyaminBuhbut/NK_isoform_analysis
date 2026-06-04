#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N nk_iso_pancan
#$ -j y

set -euo pipefail

cd "${SGE_O_WORKDIR:?}"

export PATH="$HOME/.local/bin:$PATH"
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/x86_64-redhat-linux-gnu-library/4.5}"
mkdir -p "$R_LIBS_USER"
export PANCAN_TRANSCRIPT_TABLE="${PANCAN_TRANSCRIPT_TABLE:-$HOME/tcga_transcripts/TcgaTargetGtex_rsem_isoform_tpm.gz}"
export PANCAN_CLINICAL_FILE="${PANCAN_CLINICAL_FILE:-$HOME/tcga_transcripts/pancan_clinical.csv}"

Rscript -e 'source("scripts/single_main.R"); run_significants_one_cohort("pancan", max_p_value=0.05)'
