#!/bin/bash

# Usage: sh post_threads.sh [THREADS]
THREADS=${1:-1}

# Activate the conda environment (if not already active)
eval "$(conda shell.bash hook)"
conda activate grnclust

# Find the active Python executable
PYTHON_EXEC=$(which python)
if [ -z "$PYTHON_EXEC" ]; then
  echo "Error: Python executable not found. Ensure you have activated the Conda environment."
  exit 1
fi

# Find the site-packages dir
SITE_PACKAGES_DIR=$($PYTHON_EXEC -c "import site; print(site.getsitepackages()[0])")
if [ -z "$SITE_PACKAGES_DIR" ]; then
  echo "Error: Could not locate the site-packages directory."
  exit 1
fi

# Locate arboreto/core.py
CORE_PY_PATH=$(find "$SITE_PACKAGES_DIR" -name "core.py" | grep "arboreto")
if [ -z "$CORE_PY_PATH" ]; then
  echo "Error: arboreto core.py file not found in $SITE_PACKAGES_DIR."
  exit 1
fi

# Replace any 'n_jobs': <number> with 'n_jobs': $THREADS (robust to previous value)
sed -i "s/'n_jobs': *[0-9]\+/'n_jobs': $THREADS/g" "$CORE_PY_PATH"

# Print confirmation
echo "Set n_jobs in $CORE_PY_PATH to $THREADS threads."
