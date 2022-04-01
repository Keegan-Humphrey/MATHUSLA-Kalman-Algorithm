#!/bin/bash

# specify a directory with find_dir
# any root files in that directory or it's children will be written into
# the write file (useful for running tracker on parrallel jobs)

#find_dir=/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/run6/data/
#find_dir=/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/run3/tracker_data/16_01_22/
find_dir=/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/run7/tracker_data/31_03_22/

write_file=analysis_names_6sim_7tracker_1e5.txt
#write_file=tracker_names_6.txt

current_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

touch $write_file

#find $find_dir -name '*.root' -exec echo "$current_dir/"{} >> $write_file \;
find $find_dir -name '*.root' -exec echo {} >> $write_file \;

