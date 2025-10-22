#!/bin/bash

#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --ntasks=1

#SBATCH --job-name=test_memory_consumption_oligo_generation
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

source /home/haicu/jonas.hagenberg/.bashrc
conda activate odt-ai-2
python /lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/oligo_designer_toolsuite_ai_filters/hybridization_probability/test_memory_oligos.py --region $1 --n_jobs $2