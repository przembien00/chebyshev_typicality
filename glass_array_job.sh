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
	mpirun -n 6 executable_DOUBLE.out --spinmodel=ISO --srcfile=Square_Glass_N=16 --numTimePoints=100 --beta=$d --ChebyshevCutoff=3 --numVectorsPerCore=100
wait
done
echo "DONE"
