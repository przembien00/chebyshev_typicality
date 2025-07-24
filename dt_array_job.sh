# Import data from file:
file="dt_array.txt"
if [ -e "$file" ]; then
	mapfile -t data < "$file"
else
	echo "File not found: $file"
	exit 0
fi

echo "STARTING ARRAY JOB : EXACT DIAGONALIZATION CHANGING DT"
len=${#data[@]}
for((i=0;i<$len;i++))
do
	d="${data[$i]}"
	echo "RUNNING JOB WITH TIMESTEP = "$d
	mpirun -n 4 executable_DOUBLE.out --spinmodel=ISO --srcfile=Square_NN_PBC_N=16 --beta=2 --numVectorsPerCore=2 --rescale=0.5 --seed=2 --dt=$d --fileext="dt=$d" --project="dt_therm"

wait
done
echo "DONE"
