#!/usr/bin/env bash
#
# run_batch.sh
#
# Batch wrapper around usgs_splitter.py. Processes every .tif file in
# inputs/ (ignoring .gitkeep and non-.tif files), then moves each
# processed file into inputs-completed/.
#
# Usage:
#   ./run_batch.sh

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
inputs_dir="$repo_root/inputs"
completed_dir="$repo_root/inputs-completed"

mkdir -p "$completed_dir"

shopt -s nullglob
tif_files=("$inputs_dir"/*.tif)
shopt -u nullglob

if [ ${#tif_files[@]} -eq 0 ]; then
    echo "No .tif files found in $inputs_dir"
    exit 0
fi

for tif in "${tif_files[@]}"; do
    echo "=== Processing $(basename "$tif") ==="
    python3 "$repo_root/usgs_splitter.py" "$tif"
    mv "$tif" "$completed_dir/"
    echo "=== Moved $(basename "$tif") to $completed_dir ==="
    echo
done
