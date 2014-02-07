#!/bin/sh
#SBATCH --time=03:00:00
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=6144
#SBATCH --partition=sandyb,westmere
module load parallel
srun="srun -N1 -n1"
parallel="parallel -j $SLURM_NTASKS --joblog runtask.downscaler.log --resume"
$parallel "$srun ./runtask.downscaler {1} > runtask.downscaler.{1}" ::: {1..10}
