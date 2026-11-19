#!/bin/bash
#SBATCH --job-name=SpecHLA_LOH
#SBATCH --partition=public
#SBATCH --qos=public
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vnikoone@asu.edu
#SBATCH --output=/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/05_Downstream_Analysis/02_80_patients/03_results/07_HLA_LOH_Production/logs/%A_%a.out
#SBATCH --error=/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/05_Downstream_Analysis/02_80_patients/03_results/07_HLA_LOH_Production/logs/%A_%a.err

set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate /scratch/vnikoone/999_envs/SpecHLA

BASE=/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/05_Downstream_Analysis/02_80_patients

python "$BASE/01_scripts/07_HLA_LOH/11_spechla_production.py" \
    --index "$SLURM_ARRAY_TASK_ID"

