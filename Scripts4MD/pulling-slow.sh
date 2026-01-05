#!/bin/sh
#SBATCH -t 71:10:00
#SBATCH -A naiss2025-22-1484
#SBATCH -n 32



module load GROMACS/2022.2-nsc1-gcc-9.3.0-bare

srun -n 32 gmx_mpi mdrun -ntomp 1 -dd 8 2 2 -s nnqq1-control.tpr -deffnm slow -c slow.gro
