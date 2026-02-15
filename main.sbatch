#!/bin/bash

##SBATCH --job-name=Matlab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=168:00:00
#SBATCH --mail-user=ajk8776@nyu.edu
#SBATCH --mail-type=END,FAIL
##
echo "I ran on:"
echo $SLURM_JOB_NODELIST
#
module purge
module load matlab/2024b

# export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXX)
# export MATLAB_LOG_DIR=$SLURM_JOBTMP
# can include 'parpool('local', 1)' after the cat line below

cd $SLURM_SUBMIT_DIR
cat<<EOF | matlab -nodisplay
main
exit

EOF