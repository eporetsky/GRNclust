#!/bin/bash

# [BUG] TypeError: Must supply at least one delayed object
# https://github.com/aertslab/pySCENIC/issues/592
# This issue is caused by the create_graph function in arboreto.core. 
# To resolve it, comment out the line that throws the error.
# ~/miniconda3/envs/pyscenic/lib/python3.10/site-packages/arboreto/core.py
# all_meta_df = from_delayed(delayed_meta_dfs, meta=_META_SCHEMA)

# Locate the Python executable for the active Conda environment
PYTHON_EXEC=$(which python)

if [ -z "$PYTHON_EXEC" ]; then
  echo "Error: Python executable not found. Ensure you have activated the Conda environment."
  exit 1
fi

# Locate the site-packages directory for the active Conda environment
SITE_PACKAGES_DIR=$($PYTHON_EXEC -c "import site; print(site.getsitepackages()[0])")

if [ -z "$SITE_PACKAGES_DIR" ]; then
  echo "Error: Could not locate the site-packages directory."
  exit 1
fi

# Locate the arboreto core.py file
CORE_PY_PATH=$(find "$SITE_PACKAGES_DIR" -name "core.py" | grep "arboreto")

if [ -z "$CORE_PY_PATH" ]; then
  echo "Error: arboreto core.py file not found in $SITE_PACKAGES_DIR."
  exit 1
fi

# Comment out the problematic line in core.py
sed -i 's/^.*all_meta_df = from_delayed(delayed_meta_dfs, meta=_META_SCHEMA).*$/# &/' "$CORE_PY_PATH"

echo "Successfully commented out the problematic line in arboreto core.py."