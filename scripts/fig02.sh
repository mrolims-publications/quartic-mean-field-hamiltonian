#!/usr/bin/env bash
set -euo pipefail

# Ensure script is run from the project root
if [[ ! -f Makefile ]]; then
    echo "Error: this script must be run from the project root directory (Makefile not found)."
    exit 1
fi

# ---------------- Configuration ----------------
CC="${1:-icx}"
TARGET="std_on_site_energy_history"
DOFS=(
    1024
    2048
    4096
    8192
    16384
    32768
    65536
    131072
    262144
    524288
    1048576
    2097152
    4194304
)
EPS_VALS=(1 10 100 1000)
TOTAL_TIME=100000
SEED=1313
# ------------------------------------------------

echo "Compiling ${TARGET} with CC=${CC}"
make CC="${CC}" "${TARGET}"

EXE="bin/${TARGET}.x"
if [[ ! -x "${EXE}" ]]; then
    echo "Error: executable ${EXE} not found. Did compilation fail?"
    exit 1
fi

for DOF in "${DOFS[@]}"; do
    for EPS in "${EPS_VALS[@]}"; do
        exe_cmd="${EXE} ${DOF} ${EPS} ${TOTAL_TIME} ${SEED}"
        echo "${exe_cmd}"
        ${exe_cmd}
        # ${exe_cmd} & if you want to run all of them at once
    done
done
