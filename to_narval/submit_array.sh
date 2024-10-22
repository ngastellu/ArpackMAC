#!/bin/bash
#SBATCH --account=ctb-simine
#SBATCH --mem=249G
#SBATCH --time=0-03:00
#SBATCH --array=1-300
#SBATCH --output=slurm-%a.out
#SBATCH --error=slurm-%a.err

# 1. Set up virtualenv
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt

# 2. Tell Julia to use Python version/packages in the virtualenv
julia build_pycall.jl

n=$SLURM_ARRAY_TASK_ID

# 3. Run desired computation
julia GetExtremalMOs.jl $n '40x40'
