#!/bin/bash

# specify a directory with find_dir
# any root files in that directory or it's children will be written into
# the write file (useful for running tracker on parrallel jobs)

#find_dir=~/GitHub/MATHUSLA-Kalman-Algorithm/09_12_21/18_42_42/
find_dir=~/scratch
write_file=find_names.txt

current_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

touch $write_file

find $find_dir -name '*.root' -exec echo "$current_dir/"{} >> $write_file \;

