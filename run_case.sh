#!/usr/bin/env bash
# Build and run one MFiX case.
#
# Usage:  ./run_case.sh <case_dir> <project.mfx>
# Example: ./run_case.sh gasifier_2d gasifier_2d.mfx
#
# Requires an MFiX environment on PATH (e.g. `conda activate mfix-<version>`
# first). If the case directory contains user Fortran (usr*.f), the custom
# solver is built with build_mfixsolver before running; otherwise the stock
# mfixsolver is used (e.g. the cold-flow case).
set -euo pipefail

CASE_DIR="${1:?usage: run_case.sh <case_dir> <project.mfx>}"
MFX="${2:?usage: run_case.sh <case_dir> <project.mfx>}"

cd "$CASE_DIR"

if ls usr*.f >/dev/null 2>&1; then
    echo "[run_case] building custom solver (user Fortran present)..."
    build_mfixsolver
    echo "[run_case] running $MFX ..."
    ./mfixsolver -f "$MFX"
else
    echo "[run_case] running $MFX with the stock solver ..."
    mfixsolver -f "$MFX"
fi
