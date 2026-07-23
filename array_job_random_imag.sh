#!/bin/bash

## Imaginary-time counterpart of the Random_realtime data set.
## Submit with:  sbatch --array=0-59 array_job_random_imag.sh
## Each array task takes one "N beta numConfigs h_z symm_type" line from random_imag_array.txt.
## Tasks 0-29  : h_z = 0  , symm_type = A
## Tasks 30-59 : h_z = 0.5, symm_type = C

## Mandatory:
#SBATCH --job-name=Chebyshev_Random_imag
#SBATCH --output=logs/%x-%A-%a.txt    ## File for stdout & stderr
#SBATCH --error=logs/%x-%A-%a.err
#SBATCH --time=24:00:00		## maximum runtime; hours:minutes:seconds
#SBATCH --partition=long		## choose queue

#SBATCH --ntasks=16		## number of tasks has to be = 1 for single core jobs
#SBATCH --cpus-per-task=1	## number of cpus per task has to be = 1, too!

#SBATCH --mem=24gb		## give maximum required memory in mb or gb

#SBATCH --mail-user=przemyslaw.bieniek@tu-dortmund.de	## replace mail by personal mail address
#SBATCH --mail-type=NONE	## most relevant options: NONE, BEGIN, END, FAIL

## Optional:
##SBATCH --hint=nomultithread	## deactivate Hyperthreading (recommended); for Hyperthreading comment out this line
##SBATCH --constraint=Haswell	## chose a specific feature, e.g., only nodes with Haswell-architecture
				## Feature-Output by "cat /etc/slurm/slurm.conf | grep Feature"


date

file="random_imag_array.txt"
if [ -e "$file" ]; then
	mapfile -t data < "$file"
else
	echo "File not found: $file"
	exit 0
fi

cd "../chebyshev_typicality"

echo "--- START ---"


line="${data[$SLURM_ARRAY_TASK_ID]}"
read -r N beta numConfigs h_z symm <<< "$line"

if [ -z "$symm" ]; then
	echo "No parameters for array index $SLURM_ARRAY_TASK_ID"
	exit 0
fi

echo "RUNNING JOB WITH N = $N, BETA = $beta, numCouplingConfigs = $numConfigs, h_z = $h_z, symm_type = $symm"
mpirun -n 16 executable_rc_DOUBLE.out --beta=$beta --numSpins=$N --numTimePoints=100 --spinmodel=ISO --fulldiag --numCouplingConfigs=$numConfigs --h_z=$h_z --symm_type=$symm --evol_type=imaginary --project=Random_imagtime

echo "--- END ---"
date
echo
echo "$(whoami) is leaving from $(hostname) ..."
echo
