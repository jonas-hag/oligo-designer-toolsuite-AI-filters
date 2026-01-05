#!/bin/bash

#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --mem=60G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --ntasks=1

#SBATCH --job-name=generate_real_alignment_blast_short
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

source /home/haicu/jonas.hagenberg/.bashrc
conda activate odt-ai-2
python /lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/oligo_designer_toolsuite_ai_filters/hybridization_probability/generate_real_dataset_by_region.py -c /lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/configs/hybridization_probability/generate_real_dataset_alignment_blast_short.yaml