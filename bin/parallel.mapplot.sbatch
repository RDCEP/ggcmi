#!/bin/sh
#SBATCH --time=01:00:00
#SBATCH --ntasks=40
#SBATCH --mem-per-cpu=2048
#SBATCH --partition=sandyb,westmere
module load parallel
srun="srun -N1 -n1"
parallel="parallel -j $SLURM_NTASKS --joblog runtask.mapplot.log --resume"
$parallel "$srun ./runtask.mapplot {1} > runtask.mapplot.{1}" ::: {1..40}
