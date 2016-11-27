#clean old runs
./omp > omp_runs.txt
# append to file
for((n=0; n < 10; n++));
do
	./omp >> omp_runs.txt
done