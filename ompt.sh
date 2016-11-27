#clean old runs (assuming only run once per change in implementation)
./omp > omp_runs.txt
# append to file
for((n=0; n < 4; n++));
do
	./omp >> omp_runs.txt
done

python ./quick_stats.py omp_runs.txt