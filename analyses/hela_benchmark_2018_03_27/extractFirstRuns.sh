#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --account=def-hroest
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu 10000M

unzip /scratch/hroest/tims_dia/20180327_firstRuns_data.zip -d .
