#!/bin/sh
#SBATCH --time=05:00:00
#SBATCH --ntasks=128
#SBATCH --mem-per-cpu=6144
#SBATCH --partition=sandyb
module load parallel
srun="srun -N1 -n1"
parallel="parallel -j $SLURM_NTASKS --joblog runtask.aggplot.log --resume"
$parallel "$srun ./runtask.aggplot {1} > runtask.aggplot.{1}" ::: {1..128}
