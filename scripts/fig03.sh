#!/usr/bin/env bash
set -euo pipefail

# Ensure script is run from the project root
if [[ ! -f Makefile ]]; then
    echo "Error: this script must be run from the project root directory (Makefile not found)."
    exit 1
fi

# ---------------- Configuration ----------------
CC="${1:-icx}"
TARGET="lle_history"
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
)
EPS_VALS=()
while IFS= read -r eps; do
    EPS_VALS+=( "$eps" )
done < <(awk 'BEGIN{
    exp_min=-3; exp_max=3; dexp=0.5;
    for (e=exp_min; e<=exp_max+1e-12; e+=dexp) {
        printf "%.15g\n", 10^e
    }
}')
TOTAL_TIME=100000000
RTOL=1e-7
ATOL=1e-12
NUM_VALS=500
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
        exe_cmd="${EXE} ${DOF} ${EPS} ${TOTAL_TIME} ${RTOL} ${ATOL} ${NUM_VALS} ${SEED}"
        echo "${exe_cmd}"
        ${exe_cmd}
        # ${exe_cmd} & if you want to run all of them at once
    done
done
