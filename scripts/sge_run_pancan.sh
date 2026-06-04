#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N nk_iso_pancan
#$ -j y

set -euo pipefail

cd "${SGE_O_WORKDIR:?}"

USE_CONDA_R=0

if [[ -n "${CONDA_ENV_PATH:-}" || -n "${CONDA_ENV_NAME:-}" ]]; then
  if [[ -n "${CONDA_ENV_PATH:-}" && -x "${CONDA_ENV_PATH}/bin/Rscript" ]]; then
    export PATH="${CONDA_ENV_PATH}/bin:$PATH"
    export CONDA_PREFIX="${CONDA_ENV_PATH}"
    USE_CONDA_R=1
  else
    for conda_bin_dir in \
      "${HOME}/miniconda3/bin" \
      "${HOME}/anaconda3/bin" \
      "${HOME}/mambaforge/bin" \
      "${HOME}/miniforge3/bin"
    do
      if [[ -d "$conda_bin_dir" ]]; then
        export PATH="$conda_bin_dir:$PATH"
      fi
    done

    if command -v conda >/dev/null 2>&1; then
      conda_base="$(conda info --base 2>/dev/null || true)"
      if [[ -n "$conda_base" && -f "$conda_base/etc/profile.d/conda.sh" ]]; then
        source "$conda_base/etc/profile.d/conda.sh"
      else
        eval "$(conda shell.bash hook)"
      fi
    elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
      source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
      source "$HOME/anaconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/mambaforge/etc/profile.d/conda.sh" ]]; then
      source "$HOME/mambaforge/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]]; then
      source "$HOME/miniforge3/etc/profile.d/conda.sh"
    else
      echo "Conda activation requested but conda.sh was not found." >&2
      exit 1
    fi

    if [[ -n "${CONDA_ENV_PATH:-}" ]]; then
      conda activate "$CONDA_ENV_PATH"
    else
      conda activate "$CONDA_ENV_NAME"
    fi
    USE_CONDA_R=1
  fi
fi

export PATH="$HOME/.local/bin:$PATH"
if [[ "$USE_CONDA_R" -eq 1 ]]; then
  unset R_LIBS_USER R_LIBS R_LIBS_SITE
else
  export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/x86_64-redhat-linux-gnu-library/4.5}"
  mkdir -p "$R_LIBS_USER"
fi
export PANCAN_TRANSCRIPT_TABLE="${PANCAN_TRANSCRIPT_TABLE:-$HOME/tcga_transcripts/TcgaTargetGtex_rsem_isoform_tpm.gz}"
export PANCAN_CLINICAL_FILE="${PANCAN_CLINICAL_FILE:-$HOME/tcga_transcripts/pancan_clinical.csv}"

Rscript -e 'source("scripts/single_main.R"); run_significants_one_cohort("pancan", max_p_value=0.05)'
