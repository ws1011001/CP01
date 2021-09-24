#!/bin/bash

#SBATCH -J CP01
#SBATCH -p skylake
#SBATCH --nodes=1
#SBATCH -A b222
#SBATCH -t 3-12
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH -o ./run_slurm.log/MATLAB_%j.out
#SBATCH -e ./run_slurm.log/MATLAB_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shuai.wang.notes@gmail.com
#SBATCH --mail-type=BEGIN,END

# setup path
idir='/home/swang/simages'        # singularity images directory
mdir='/scratch/swang/agora/CP01'  # the project main directory
scripts='/CP01/scripts'           # scripts folder in Sy

# set an array of scripts
scripts_matlab[0]='ps05_TFA_individual_high_gamma_mia.m'
script_id=0

# run MATLAB script
script_run=${scripts_matlab[$script_id]}
matlab_log="$mdir/scripts/run_MATLAB.log"  # records of running MATLAB scripts with slurm
SECONDS=0
echo -e "========== Start running $script_run with singularity at $(date) =========="
echo -e "$(date) : $script_run" >> $matlab_log
singularity exec --bind $mdir:/CP01 $idir/nidebian-2.1 \
  matlab -nodisplay -nodesktop -r "try;cd('$scripts');run('$script_run');catch ME;fprintf('ERROR: %s\n', ME.message);end;exit"
echo -e "========== Finish $script_run with singularity at $(date) =========="
duration=$SECONDS
echo -e "$(date) : $script_run - $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." >> $matlab_log
