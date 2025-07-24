# Import data from file:
file="rescale_array.txt"
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
	echo "RUNNING JOB WITH RESCALE = "$d
	mpirun -n 4 executable_DOUBLE.out --spinmodel=ISO --srcfile=Square_NN_PBC_N=16 --dt=0.001 --beta=0.1 --numVectorsPerCore=5 --rescale=$d --project=J_constT
wait
done
echo "DONE"
