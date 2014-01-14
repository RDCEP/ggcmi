# run after job is completed
cat timestamps_*.txt > timestamps.txt
rm timestamps_*.txt
mv timestamps.txt ..
rm runtask.aggregator.*
rm slurm*
