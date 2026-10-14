#!/bin/bash
#SBATCH --job-name=shap_plots             # Job name
#SBATCH --output=log/aucsPlot_%A.out
#SBATCH --error=log/aucsPlot_%A.err
#SBATCH --time=01:00:00                # Maximum runtime (hh:mm:ss)
#SBATCH --mem=5G
#SBATCH --ntasks=1                    # Number of tasks (1 task for a single script)
#SBATCH --nodes 1
#SBATCH --cpus-per-task=5             # Number of CPU cores per task


# Load the conda module if necessary
module load Conda

# Activate the Conda environment
conda activate /lisc/data/scratch/ccr/conda_envs/rML_env

results=$1
group1=$2
group2=$3
size1=$4
size2=$5

python plot_nestedCV_ROC-AUCS.py  --joblib-dir $results --group1 $group1 --size1 $size1 --group2 $group2 --size2 $size2 --out-dir figures/ --out-base roc_curve_${group1}_${group2}_nestedCV --color-indiv '#e5b8b0'  --color-mean '#cb7060'  --color-ci '#e5b8b0'  --color-rand gray
