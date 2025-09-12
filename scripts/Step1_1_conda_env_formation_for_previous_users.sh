#!/bin/bash
set -e

# 1. Detect Mamba/Conda
if command -v mamba >/dev/null 2>&1; then
  CREATE_CMD="mamba"
elif command -v conda >/dev/null 2>&1; then
  CREATE_CMD="conda"
else
  echo "ERROR: Neither 'mamba' nor 'conda' found in PATH. Follow Step1"
  exit 1
fi

# 2. Initialize conda
if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  eval "$(conda shell.bash hook)" || true
fi

# 3. Create a virtual environment
ENV_NAME="enhancer-env"
"${CREATE_CMD}" create -n "${ENV_NAME}" python=3.9 numpy pandas r-base -y

# 7. Launch a new shell with the "enhancer-env1" environment activated
echo "Environment 'enhancer-env' created. Opening new shell with it activated."
exec bash --rcfile <(cat <<'RC'

# initialize shell for mamba/conda in this new shell only
if command -v mamba >/dev/null 2>&1; then
  eval "$(mamba shell hook --shell bash)"
elif command -v conda >/dev/null 2>&1; then
  source <(conda shell.bash hook) >/dev/null 2>&1 || true
fi

# activate env (prefer mamba if available)
if command -v mamba >/dev/null 2>&1; then
  mamba activate enhancer-env1
else
  conda activate enhancer-env1
fi
RC
)
