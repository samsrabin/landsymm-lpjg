#!/bin/bash
#SBATCH --job-name=plum
#SBATCH --time=10:00:00
#SBATCH --output=plum.out
#SBATCH --error=plum.err
#SBATCH --cpus-per-task=1
#SBATCH --partition=broadwl
#SBATCH --mem-per-cpu=20000
module load python/3.6.1+intel-16.0
python plum.py
