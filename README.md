# MATHUSLA-Tracker
Tracking software for MATHUSLA Experiment

## Introduction
This repository, including this documentation, is a modification of the following code from: https://github.com/seg188/MATHUSLA-MLTracker 

In addition to the previously existing Maximum Likelihood Linear Tracking Algorithm, this repository also contains a Kalman Filter Tracking Algorithm.

## Installation

Most of the dependencies and necessary libraries for this code are available through the anaconda software. If you do not have anaconda, it can be installed following the directions here: https://docs.anaconda.com/anaconda/install/

Once anaconda is installed, the enviornment can be created using the .yml file available in this repository as follows. From the top directory:

```bash
$ conda env create -f env/environment.yml
$ conda activate tracker
```

Now, the project can be built using cmake:

```bash
$ cd /tracker/
$ mkdir build
$ cd build
$ cmake ../ 
$ make 
```

At this point, the tracker executable is available in the build directory. Note that the /build/ directory MUST be placed in the /tracker/ directory. 


## Running the Tracker

The tracker requires two command line arguments, the path to an input file, and the path to which the output file should be written. The input file should be the output from a MATHUSLA Mu-Simulation run, a Geant4 based simulation of particle passage through the MATHUSLA detector. 

The Mu-Simulation repository can be found here: https://github.com/MATHUSLA/Mu-Simulation

An example command to run the tracker:

```bash
$ ./tracker path_to_input_file path_to_write_output 
```
A script for automating series runs of the tracker, and further documentation about it, is located in the /run/ directory. 



