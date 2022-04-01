#!/bin/bash
#SBATCH --account=SAMPLE_NAME
#SBATCH --time=12:00:00
#SBATCH --job-name=W_sim
#SBATCH --output=/home/keeganh/projects/rrg-mdiamond/keeganh/job_test/W_sample_dir/run7/%x-%j.out
#SBATCH --array=1-8000%1000
#SBATCH --mem-per-cpu=512M

echo 'Hello, world!'

sleep $SLURM_ARRAY_TASK_ID

/home/keeganh/GitHub/Mu-Simulation/simulation -s /home/keeganh/GitHub/Mu-Simulation/studies/box/background/run_pythia_w.mac

echo 'done'
