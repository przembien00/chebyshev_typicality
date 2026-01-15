# Import data from file:
file="temperature_array.txt"
if [ -e "$file" ]; then
	mapfile -t data < "$file"
else
	echo "File not found: $file"
	exit 0
fi

echo "STARTING ARRAY JOB : EXACT DIAGONALIZATION CHANGING T"
len=${#data[@]}
for((i=0;i<$len;i++))
do
	d="${data[$i]}"
	echo "RUNNING JOB WITH BETA = "$d
	mpirun -n 8 executable_rc_DOUBLE.out --spinmodel=ISO --numSpins=13 --numCouplingConfigs=5 --numTimePoints=50 --beta=$d --fulldiag --project="Test" --symm_type=A
wait
done
echo "DONE"