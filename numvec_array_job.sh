# Import data from file:
file="numvec_array.txt"
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
	n="${data[$i]}"
	echo "RUNNING JOB WITH NUMBER OF VECTORS = "$n
	mpirun -n 4 executable_DOUBLE.out --spinmodel=ISO --srcfile=Square_NN_PBC_N=16 --beta=0.5 --numVectorsPerCore=$n --rescale=0.5 --dt=0.001 --fileext="numVecPerCore=$n" --project="NumVec_DOUBLE"

wait
done
echo "DONE"
