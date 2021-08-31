#!/bin/bash

if [ $# -eq 0 ]
then
    echo "Usage: ./Study <path to write directory>"
    exit 1
fi

echo "Hello user!"

# make write directories
day=$(date +"%d_%m_%y")
time=$(date +"%H_%M_%S")

# check if parent directory (date) exists
if [ ! -d "$1/$day" ]
then
	mkdir "$1/$day"
fi

# create directories for the run
mkdir "$1/$day/$time"
mkdir "$1/$day/$time/trees"
mkdir "$1/$day/$time/prints"
mkdir "$1/$day/$time/plots"

# get script directory (so run.sh can be executed from a different directory)
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
echo "script directory is: $script_dir"

# check for files left from previous run
# (necessary after keyboard interrupt)
Done="$script_dir/done.txt"
if [ -f "$Done" ]
then
	rm "$Done"
	echo "done.txt has been removed"
fi

# move to /path/to/tracker/
cd "$script_dir/../"

# keep source for version control
zip -r include.zip ./include/
zip -r src.zip ./src/
mv include.zip "$1/$day/$time"
mv src.zip "$1/$day/$time"

# file containing paths to input trees
input="./run/DataNames.txt"

n=0

while IFS= read -r line
do
	echo "current dataset: $line"

	c=0

	if [ -f "$line" ]
	then
		while [ ! -f "./run/done.txt" ]
		do
			#modify run parameters
			cd ./run
			python mod_globs.py "$c" "$n"
			cd ../

			#build and run tracker
			cd build
#			make
			./tracker "$line" ../
			cd ../

            # move and rename output files
			mv stat0.root "$1/$day/$time/trees/stat_${n}_${c}.root"
#			mv ./build/print.txt "$1/$day/$time/prints/print_${n}_${c}.txt"

            # run cutflow script on output tree
			python "$script_dir/../../analysis/cutflow.py" "$1/$day/$time/trees/stat_${n}_${c}.root"

			((c++))
		done
	else
		echo "File Doesn't exist, moving on..."
	fi

	((n++))

	Done=./run/done.txt
	if [ -f "$Done" ]
	then
		rm "$Done"
	fi

done < "$input"

# for reference in plots*.zip file!
cp ./run/*.txt "$1/$day/$time/"
cp ./run/*.txt "$1/$day/$time/plots"

# make directories to write plots into
mkdir "$1/$day/$time/plots/NumTracks"
mkdir "$1/$day/$time/plots/NumTracks_k"
mkdir "$1/$day/$time/plots/NumTracks_k_m"
mkdir "$1/$day/$time/plots/Track_numHits"
mkdir "$1/$day/$time/plots/Track_k_numHits"
mkdir "$1/$day/$time/plots/local_chi_f"
mkdir "$1/$day/$time/plots/local_chi_s"
mkdir "$1/$day/$time/plots/Track_k_smooth_chi_sum"
mkdir "$1/$day/$time/plots/vertex_k_f_beta"
mkdir "$1/$day/$time/plots/vertex_k_s_beta"
mkdir "$1/$day/$time/plots/Track_k_beta"
mkdir "$1/$day/$time/plots/Track_k_beta_err"

python ./run/plotting.py "$1/$day/$time/"

cd "$1"

zip -r "$day/$time/plots_$day-$time.zip" "$day/$time/plots/"

# if nohup was used, include output file
nohup="$script_dir/nohup.out"

if [ -f "$nohup" ]
then
        mv "$nohup" "$1/$day/$time/"
fi

echo "write directory: $1/$day/$time/"
echo "All done!"
