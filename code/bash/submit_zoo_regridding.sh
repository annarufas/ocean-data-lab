#!/bin/bash

#SBATCH --job-name=mesozoo
#SBATCH --output=log_zoo_%A_%a.log # use both job ID and array index in the log
#SBATCH --time=50:00:00
#SBATCH --partition=long
#SBATCH --clusters=htc
#SBATCH --ntasks=1 
#SBATCH --mem=50G
#SBATCH --array=1-5 # for running the same analysis several times with different input files

# Load the necessary module
module load MATLAB/R2022a

# Set MATLABPATH
export MATLABPATH=/data/eart-slam-dunk/wolf4894/Datasets/zooplankton/CMIP6/ 
export MATLABPATH=/data/eart-slam-dunk/wolf4894/Datasets/grid/ 
export MATLABPATH=/data/eart-slam-dunk/wolf4894/Resources/

# Execute MATLAB script with parameters
matlab -nodesktop -nosplash -r \
	"regridZooplanktonConcentrationFromCMIP6('mesozooClimatologyNative_dz${SLURM_ARRAY_TASK_ID}', \
	'mesozooClimatologyRegular_dz${SLURM_ARRAY_TASK_ID}','grid_MITgcm_2p8deg');exit;"
