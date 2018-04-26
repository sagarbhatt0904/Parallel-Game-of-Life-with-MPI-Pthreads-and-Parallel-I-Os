#!/bin/bash
#
# USAGE: /full/path/to/./q_MC.out [--help]
#                                 [--init dimension [outfile]]
#                                 [--nonstop dimension outfile steps [increment]]
#                                 [infile [outfile] steps [increment]]
#
#SBATCH --job-name=ppcAssignment4
#SBATCH --mail-type=END
#SBATCH --mail-user=bhatts8@rpi.edu
#SBATCH --partition debug
#SBATCH -t 00:05:00
#SBATCH -N 16
#SBATCH -n 1024
#SBATCH --overcommit
#SBATCH -o /gpfs/u/home/PCP7/PCP7sgrb/scratch/assignment3/out.log


srun --runjob-opts="--mapping TEDCBA" /gpfs/u/home/PCP7/PCP7sgrb/scratch/assignment4/assignment4-5
