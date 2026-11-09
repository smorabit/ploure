#!/usr/bin/env bash
# R CMD check on a clean export of the tracked files: R CMD build copies the whole directory before applying
# .Rbuildignore, so building in place would copy the pixi environment (.pixi/, several GB)
set -euo pipefail
root=$(git -C "$(dirname "$0")/.." rev-parse --show-toplevel)
work=$(mktemp -d "${TMPDIR:-/tmp}/ploure-check.XXXXXX")
mkdir -p "$work/ploure"
git -C "$root" ls-files --cached --others --exclude-standard | tar -C "$root" -cf - -T - | tar -C "$work/ploure" -xf -
cd "$work"
R CMD build ploure
R CMD check --no-manual ploure_*.tar.gz
echo "check directory: $work"
