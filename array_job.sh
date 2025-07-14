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
	mpirun -n 4 executable_FLOAT.out --spinmodel=ISO --srcfile=Square_NN_PBC_N=16 --numTimePoints=201 --beta=$d --numVectorsPerCore=1 --rescale=0.5
wait
done
echo "DONE"
