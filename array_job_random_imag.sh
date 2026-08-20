#!/bin/bash

## Imaginary-time counterpart of the Random_realtime data set.
## Submit the full data set with:  sbatch array_job_random_imag.sh
## Resubmit selected files with:   sbatch --array=<indices> array_job_random_imag.sh
## Each array task takes one "N beta numConfigs h_z symm_type" line from random_imag_array.txt.
## Tasks 0-29  : h_z = 0  , symm_type = A
## Tasks 30-59 : h_z = 0.5, symm_type = C

## Mandatory:
#SBATCH --job-name=Chebyshev_Random_imag
#SBATCH --output=logs/%x-%A-%a.txt    ## File for stdout & stderr
#SBATCH --error=logs/%x-%A-%a.err
#SBATCH --array=0-59
#SBATCH --time=48:00:00		## maximum runtime; hours:minutes:seconds
#SBATCH --partition=long		## choose queue

#SBATCH --ntasks=16		## number of tasks has to be = 1 for single core jobs
#SBATCH --cpus-per-task=1	## number of cpus per task has to be = 1, too!

#SBATCH --mem=32gb		## give maximum required memory in mb or gb

#SBATCH --mail-user=przemyslaw.bieniek@tu-dortmund.de	## replace mail by personal mail address
#SBATCH --mail-type=NONE	## most relevant options: NONE, BEGIN, END, FAIL

## Optional:
##SBATCH --hint=nomultithread	## deactivate Hyperthreading (recommended); for Hyperthreading comment out this line
##SBATCH --constraint=Haswell	## chose a specific feature, e.g., only nodes with Haswell-architecture
				## Feature-Output by "cat /etc/slurm/slurm.conf | grep Feature"

set -euo pipefail

date

## Use the same disorder realization for every rerun. For a given N this also
## matches the coupling configurations across beta and h_z.
coupling_seed=0
expected_tasks=60
file="random_imag_array.txt"

if [ -f "$file" ]; then
	mapfile -t data < "$file"
else
	echo "File not found: $file" >&2
	exit 1
fi

if [ "${#data[@]}" -ne "$expected_tasks" ]; then
	echo "Expected $expected_tasks parameter lines in $file, found ${#data[@]}" >&2
	exit 1
fi

task_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set; submit this as a Slurm array job}"
if ! [[ "$task_id" =~ ^[0-9]+$ ]] || [ "$task_id" -ge "${#data[@]}" ]; then
	echo "Invalid array index: $task_id" >&2
	exit 1
fi

echo "--- START ---"


line="${data[$task_id]}"
read -r N beta numConfigs h_z symm <<< "$line"

if [ -z "$symm" ]; then
	echo "No parameters for array index $task_id" >&2
	exit 1
fi

echo "RUNNING JOB WITH N = $N, BETA = $beta, numCouplingConfigs = $numConfigs, h_z = $h_z, symm_type = $symm, seed = $coupling_seed"
mpirun -n 16 ./executable_rc_DOUBLE.out \
	--beta="$beta" \
	--numSpins="$N" \
	--numTimePoints=100 \
	--spinmodel=ISO \
	--fulldiag \
	--numCouplingConfigs="$numConfigs" \
	--h_z="$h_z" \
	--symm_type="$symm" \
	--evol_type=imaginary \
	--project=Random_imagtime_new

echo "--- END ---"
date
echo
echo "$(whoami) is leaving from $(hostname) ..."
echo
