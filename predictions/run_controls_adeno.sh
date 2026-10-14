#!/bin/bash
#SBATCH --job-name=classif_all             # Job name
#SBATCH --output=log/3classif_all_%A_%a.out
#SBATCH --error=log/3classif_all_%A_%a.err
#SBATCH --time=06:00:00                # Maximum runtime (hh:mm:ss)
#SBATCH --mem=2G
#SBATCH --ntasks=1                    # Number of tasks (1 task for a single script)
#SBATCH --nodes 1
#SBATCH --cpus-per-task=5             # Number of CPU cores per task


# Load the conda module if necessary
module load Conda

# Activate the Conda environment
conda activate /lisc/data/scratch/ccr/conda_envs/rML_env/

# Compute the line number corresponding to the current array index
# (Assuming SLURM_ARRAY_TASK_ID is 0-indexed; adjust if necessary)
line_number=$(($SLURM_ARRAY_TASK_ID))

# Read the seed from the file using sed (prints the nth line)
seed=$(sed -n "${line_number}p" $1) #"$seeds_file")

echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Using seed: $seed"

# Change directory to the location of your script
cd /lisc/data/scratch/ccr/CR_projects/CRC/scripts
# Run your Python script
python classification_train_test_v2_cases.py --config /lisc/data/scratch/ccr/CR_projects/CRC/Metadata/config_file_controls_adeno.yaml --seed $seed --run_nested_cv True  --subgroup all --with_oligos True --with_additional_features False --prevalence_threshold_min 2.0 --prevalence_threshold_max 98.0   --out_name controls_adeno  --out_dir /lisc/data/scratch/ccr/CR_projects/CRC/controls_adeno --no_additional_train_test_data FALSE --outer_cv_split 5 --inner_cv_split 5 --model_type random-forest --train '{"group_test": ["Controls", "Adenoma"]}'
