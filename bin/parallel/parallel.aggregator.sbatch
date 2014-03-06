#!/bin/sh
#SBATCH --time=08:00:00
#SBATCH --ntasks=128
#SBATCH --mem-per-cpu=6144
#SBATCH --partition=sandyb
module load parallel
srun="srun -N1 -n1"
parallel="parallel -j $SLURM_NTASKS --joblog runtask.aggregator.log --resume"
$parallel "$srun ./runtask.aggregator {1} > runtask.aggregator.{1}" ::: {1..128}
