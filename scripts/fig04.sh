#!/usr/bin/env bash
set -euo pipefail

# Ensure script is run from the project root
if [[ ! -f Makefile ]]; then
    echo "Error: this script must be run from the project root directory (Makefile not found)."
    exit 1
fi

# ---------------- Configuration ----------------
CC="${1:-icx}"
TARGET="p_PDF"
DOF=1024
EPS=100
TOTAL_TIME_VALS=(
    100000
    1000000
    10000000
    100000000
)
SEED=1313
# ------------------------------------------------

echo "Compiling ${TARGET} with CC=${CC}"
make CC="${CC}" "${TARGET}"

EXE="bin/${TARGET}.x"
if [[ ! -x "${EXE}" ]]; then
    echo "Error: executable ${EXE} not found. Did compilation fail?"
    exit 1
fi

for TOTAL_TIME in "${TOTAL_TIME_VALS[@]}"; do
    exe_cmd="${EXE} ${DOF} ${EPS} ${TOTAL_TIME} ${SEED}"
    echo "${exe_cmd}"
    ${exe_cmd}
    # ${exe_cmd} & if you want to run all of them at once
done
