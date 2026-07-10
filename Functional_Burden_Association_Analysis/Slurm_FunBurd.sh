#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script to run FunBurd Association Test Across Traits and Genelists
# ---------------------------------------------------------------------
 
# 
#SBATCH --job-name=FunBurd
#SBATCH --time=02:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1GB
#SBATCH --account=rrg-jacquese
#SBATCH --output=_log/%x_%A_%a.out
#SBATCH --error=_log/%x_%A_%a.err
#SBATCH --mail-user=sayeh.kazem@umontreal.ca
#SBATCH --mail-type=ALL

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------

source "$HOME/my_firth_env/bin/activate"
#module load StdEnv/2023 scipy-stack/2023b python/3.10.13 ipython-kernel/3.10.13
module load gcc/12.3
module load openmpi/4.1.5
module load glost/0.3.1

#--------------------------------------------------------------------------------------------
# Paramter files
cd links/projects/rrg-jacquese/sayeh94/NewProject_2023/All_Wiener/FunBurd_Pipeline_on_ToyDataset

srun glost_launch Args_File_FunBurd.txt
