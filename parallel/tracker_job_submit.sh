#!/bin/bash
#SBATCH --account=SAMPLE_NAME
#SBATCH --time=3:00:00
#SBATCH --job-name=par_tracker
#SBATCH --output=/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/run3/%x-%j.out
#SBATCH --array=1-10
#SBATCH --mem-per-cpu=128M

#-------------------------------------------------------------------------------------------------
echo 'Hello, track maker!'

sleep $SLURM_ARRAY_TASK_ID

#-------------------------------------------------------------------------------------------------

# where am I writting to?
out_dir="/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir//"

# what is the file with the names of the trees to track?
in_names="/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/Parallel/tracker_names.txt"

# where is the tracker executable?
exec_dir="/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/tracker/build/"

# How many series tracker runs per job
MAX_SERIES=10

#-------------------------------------------------------------------------------------------------

track_dir="$out_dir/tracker_data/"
if [ ! -d "$track_dir" ]
then
	mkdir "$track_dir"
fi

# make write directories
day=$(date +"%d_%m_%y")

# check if parent directory (date) exists
day_dir="$track_dir/$day"
if [ ! -d "$day_dir" ]
then
        mkdir "$day_dir"
fi

#make an array job specific dir to write to
ID_dir="$day_dir/$SLURM_ARRAY_TASK_ID"
if [ ! -d "$ID_dir" ]
then
        mkdir "$ID_dir"
fi

# count series runs
n=0

while [ $n -lt $MAX_SERIES ]
do
	line="$(<$in_names cut -d $'\n' -f $(($MAX_SERIES*$SLURM_ARRAY_TASK_ID+$n)))"

	#echo "$line"

	"$exec_dir/tracker" "$line" "$ID_dir"

	mv "$ID_dir/stat0.root" "$ID_dir/stat_${SLURM_ARRAY_TASK_ID}_${n}.root"

	((n++))

	# redundant ... but for safety!
	if [ $n -eq $MAX_SERIES ]
	then
		break
	fi
done

echo 'Finito!'
