#!/bin/bash

#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --mem=10G
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --ntasks=1

#SBATCH --job-name=split_fasta_files
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

source /home/haicu/jonas.hagenberg/.bashrc
conda activate odt-ai-2
python /lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/oligo_designer_toolsuite_ai_filters/hybridization_probability/generate_real_dataset_split_fasta_file.py