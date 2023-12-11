#!/bin/bash -l
#SBATCH -J pyadfjob
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=5GB
#SBATCH --time=01:00:00
#SBATCH -A plgqcembed-cpu
#SBATCH -p plgrid-now
#SBATCH --output="output.out"
#SBATCH --error="error.err"

# adapt this:
scratch=$SCRATCH/myname-space/qm_calcs/pyadf/molname/spinorbitZORA/PBE/TZ2P
mkdir -p $scratch

project=generic

# normally, this should not need to be adapted:
srun /bin/hostname

module purge
module use /net/pr2/projects/plgrid/plggqcembed/groupmodules
module load pyadf-master

cd $SLURM_SUBMIT_DIR
config='/net/pr2/projects/plgrid/plggqcembed/devel/tools/pyadf-jobrunner.conf'
# to save results add -s flag:
# pyadf -s --jobrunnerconf $config $project.pyadf
# otherwise:
pyadf --jobrunnerconf $config $project.pyadf

cd $SLURM_SUBMIT_DIR
