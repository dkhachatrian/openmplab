#clean old runs
./seq > seq_runs.txt
# append to file
for((n=0; n < 4; n++));
do
	./seq >> seq_runs.txt
done

python ./quick_stats.py seq_runs.txt