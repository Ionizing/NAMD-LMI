#!/bin/bash
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --ntasks-per-node=48
#SBATCH --partition=part2

module load intel/2020u4
mpirun /data/app/Hefeitest/vasp.6.3.0_wan_soc_vtst_bin/bin/vasp_ncl
