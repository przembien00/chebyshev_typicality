#!/bin/bash

## Real-time extension of Data/Random_realtime to Tmax=15.
## Submit all 60 runs with: sbatch array_job_random_realtime_tmax15.sh
## Resubmit a subset with: sbatch --array=<indices> array_job_random_realtime_tmax15.sh
## Tasks 0-29: h_z=0, symm_type=A.  Tasks 30-59: h_z=0.5, symm_type=C.
## The original Tmax=10 data use 100 points; 150 points retain its nominal 0.1 spacing.

#SBATCH --job-name=Chebyshev_Random_rt15
#SBATCH --output=logs/%x-%A-%a.txt
#SBATCH --error=logs/%x-%A-%a.err
#SBATCH --array=0-59
#SBATCH --time=48:00:00
#SBATCH --partition=long
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --mail-user=przemyslaw.bieniek@tu-dortmund.de
#SBATCH --mail-type=NONE

set -euo pipefail

data_file="random_realtime_tmax15_array.txt"
expected_tasks=60
coupling_seed=0

if [[ ! -f "$data_file" ]]; then
    echo "File not found: $data_file" >&2
    exit 1
fi

mapfile -t data < <(sed '/^[[:space:]]*#/d; /^[[:space:]]*$/d' "$data_file")
if [[ "${#data[@]}" -ne "$expected_tasks" ]]; then
    echo "Expected $expected_tasks parameter lines in $data_file, found ${#data[@]}" >&2
    exit 1
fi

task_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set; submit this as a Slurm array job}"
if ! [[ "$task_id" =~ ^[0-9]+$ ]] || (( task_id >= ${#data[@]} )); then
    echo "Invalid array index: $task_id" >&2
    exit 1
fi

read -r N beta num_configs h_z symm_type <<< "${data[$task_id]}"
if [[ -z "${symm_type:-}" ]]; then
    echo "Malformed parameter line for array index $task_id" >&2
    exit 1
fi

echo "Running N=$N beta=$beta numCouplingConfigs=$num_configs h_z=$h_z symm_type=$symm_type seed=$coupling_seed"
mpirun -n 16 ./executable_rc_DOUBLE.out \
    --beta="$beta" \
    --numSpins="$N" \
    --numTimePoints=150 \
    --Tmax=15 \
    --spinmodel=ISO \
    --fulldiag \
    --numCouplingConfigs="$num_configs" \
    --h_z="$h_z" \
    --symm_type="$symm_type" \
    --evol_type=real \
    --project=Random_realtime_tmax15
