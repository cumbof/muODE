#!/usr/bin/env bash
# Snakemake post-deploy hook for envs/carvefungi.yaml.
#
# Runs after the conda env is created, with the env activated and CWD set to the
# Snakemake working directory (the repo root).  Installs the local muode package
# editable so workflow/scripts/reconstruct_one.py can import muode here.
#
# This replaces a `pip: ["."]` entry in the env YAML, which was unreliable:
# conda runs its pip step from an unspecified directory, so a relative `.` did
# not resolve to the repo root ("Directory '.' is not installable").
set -euo pipefail

if [[ ! -f pyproject.toml ]]; then
    echo "muode post-deploy: pyproject.toml not found in $(pwd)." >&2
    echo "Run snakemake from the repository root so muode can be installed." >&2
    exit 1
fi

pip install -e .
