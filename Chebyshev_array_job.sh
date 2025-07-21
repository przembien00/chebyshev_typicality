# Import data from file:
file="Chebyshev_array.txt"
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
	echo "RUNNING JOB WITH CHEBYSHEV ORDER = "$n
	mpirun -n 4 executable_FLOAT.out --spinmodel=ISO --srcfile=Square_NN_PBC_N=16 --numTimePoints=300 --beta=5 --numVectorsPerCore=1 --rescale=0.5 --seed=10 --ChebyshevCutoff=$n --fileext="Cheb_cut=$n"

wait
done
echo "DONE"
