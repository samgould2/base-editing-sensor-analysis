#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1                      # Number of tasks (really number of CPU Cores/task). Don't specify more than 16 unless approved by the system admin
#SBATCH --array=1-6    #change this to match up with the number of parallel computing jobs; also generate a config file 
#SBATCH --mail-type=END           # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE 
#SBATCH --mail-user=samgould@mit.edu           # Email to which notifications will be sent. Equivalent to the -M option in SGE. You must replace [] with your email address.
#SBATCH --exclude=c[5-22]
#SBATCH --nice=100000
#############################################

module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit

#create a new virtual environment if needed; this one is fine from previously
conda activate /home/samgould/.conda/envs/sensor_lib_sg

#go to correct directory
cd /net/bmc-lab2/data/lab/sanchezrivera/samgould/240624San

#access the config file
#NEED TO CHANGE THE CONFIG FILE FOR YOUR SPECIFIC RUN
config=./DIEGO_IDR_CONFIG_MISMATCH_PIPELINE.txt

# Extract R1_FILE name for the current $SLURM_ARRAY_TASK_ID
R1_FILE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
R2_FILE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
folder_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)

# Paramters for the script (NO SPACES ALLOWED)
splitby='protospacer'
proto_mismatches_allowed=2
bc_len=15
sensor_len=42

python3 sensor_extraction_guide_counts.py IDR_ABE_library_columns_renamed.csv ${R1_FILE} ${R2_FILE} ${splitby} ${proto_mismatches_allowed} ${bc_len} ${sensor_len} -o ${folder_name}
