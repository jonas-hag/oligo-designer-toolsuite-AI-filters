#!/bin/bash

#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=30:00
#SBATCH --ntasks=1

#SBATCH --job-name=generate_real_select_alignments
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

source /home/haicu/jonas.hagenberg/.bashrc
conda activate odt-ai-2
python /lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/oligo_designer_toolsuite_ai_filters/hybridization_probability/generate_real_dataset_select_alignments_updated.py -c /lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/configs/hybridization_probability/generate_real_dataset_select_alignments.yaml