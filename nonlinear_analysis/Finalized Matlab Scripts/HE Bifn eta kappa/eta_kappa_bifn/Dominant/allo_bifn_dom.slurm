#!/bin/bash
#SBATCH --job-name=matlab
#SBATCH --account=rgutenk
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=08:00:00

module load matlab/r2023b

matlab -nodisplay -nosplash < new_allo_bifn_data_dominant.m > allo_bifn_dom.txt
