The contents of this directory, allow series runs of the MATHUSLA tracker across a set of input datasets 
and run parmameters. The directory is intended to be portable to any version of the tracker by simply
copying /run/ into /tracker/ and correctly placing the build directory.

-- To Prepare:

Program assumes that the build directory is /path/to/tracker/build/. If building from scratch:

mkdir /path/to/tracker/build/
cd /path/to/tracker/build/
cmake ../
make

-- To Run:

./run.sh <path to write directory>

-- Other Notes:

DataNames.txt specifies the datasets over which the tracker will be run.
Paramaters are specified in Parameters.txt, they will be run for each dataset in series.

./run.sh will modify ../include/globals.hh by running ./mod_glob.py. 
./run.sh contains a while loop that is broken when ./done.txt is created by ./mod_globs.py. 

output directories are organised by <path to write directory>/<date>/<time>/

./plotting.py is run to generate ROOT plots with each of the output trees for the run. These are stored
in <path to write directory>/<date>/<time>/plots/. This directory is automatically zipped for downloading
and viewing if tracker has been run on a cluster.

output files are organised as:

<file>_i_j.*

where:
- i - indexes the dataset (ordered as in ./DataNames.txt)
- j - indexes the parameter list (ordered as in ./Parameters.txt)

