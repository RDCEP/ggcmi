#!/bin/sh
#SBATCH --time=03:00:00
#SBATCH --ntasks=128
#SBATCH --mem-per-cpu=6144
#SBATCH --partition=sandyb,westmere
module load parallel
srun="srun -N1 -n1"
parallel="parallel -j $SLURM_NTASKS --joblog runtask.biascorrect.log --resume"
$parallel "$srun ./runtask.biascorrect {1} > runtask.biascorrect.{1}" ::: {1..128}
