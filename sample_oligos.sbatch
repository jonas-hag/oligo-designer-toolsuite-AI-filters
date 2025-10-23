#!/bin/bash

#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --mem=120G
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --ntasks=16

#SBATCH --job-name=sample_oligos
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

source /home/haicu/jonas.hagenberg/.bashrc
conda activate odt-ai-2
python /lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/oligo_designer_toolsuite_ai_filters/hybridization_probability/generate_real_dataset_sample_oligos.py -c /lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/configs/hybridization_probability/generate_real_dataset_sample_oligos.yaml