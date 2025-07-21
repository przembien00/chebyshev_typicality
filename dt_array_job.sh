# Import data from file:
file="dt_array.txt"
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
	echo "RUNNING JOB WITH THERMAL DISCRETIZATION = "$d
	mpirun -n 4 executable_DOUBLE.out --spinmodel=ISO --srcfile=Square_NN_PBC_N=16 --numTimePoints=300 --beta=0.5 --numVectorsPerCore=1 --rescale=0.5 --seed=1 --dtThermalize=$d --fileext="dt=$d" --project="dt_therm"

wait
done
echo "DONE"
