#!/bin/bash
#SBATCH --account=SAMPLE_NAME
#SBATCH --time=3:00:00
#SBATCH --job-name=par_analyser
#SBATCH --output=/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/run3/%x-%j.out
#SBATCH --array=1-19
#SBATCH --mem-per-cpu=128M

#-------------------------------------------------------------------------------------------------
echo 'Hello, analysis doer!'

#sleep $SLURM_ARRAY_TASK_ID

#-------------------------------------------------------------------------------------------------

# where am I writting to?
out_dir="/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/run3/"

# what is the file with the names of the trees to analyse?
in_names="/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/Parallel/analysis_names.txt"

# where is the analysis script?
exec_dir="/home/keeganh/GitHub/MATHUSLA-Kalman-Algorithm/analysis/"

# How many series analysis runs per job
MAX_SERIES=10

#-------------------------------------------------------------------------------------------------

analysis_dir="$out_dir/analysis_data/"
if [ ! -d "$analysis_dir" ]
then
	mkdir "$analysis_dir"
fi

# make write directories
day=$(date +"%d_%m_%y")

# check if parent directory (date) exists
day_dir="$analysis_dir/$day"
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

	python "$exec_dir/cutflow.py" "$line" "$ID_dir" "$(($MAX_SERIES*$SLURM_ARRAY_TASK_ID+$n))"

	((n++))

	# redundant ... but for safety!
	if [ $n -eq $MAX_SERIES ]
	then
		break
	fi
done

echo 'Finito!'
