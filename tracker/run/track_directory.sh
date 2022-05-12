#!/bin/bash

# This script takes an input and an output directory, runs tracker on every root file in the input directory, 
# and each time renames the stat0.root to "Trk" + the input file name 
# the ##*/ is just to cut out the path info from the file name if necessary

# Written by https://github.com/serramilli

for file in "${1}"/*.root
do
    ../build/tracker "$file" "${2}"
    mv "${2}/stat0.root" "${2}/Trk${file##*/}"
done
