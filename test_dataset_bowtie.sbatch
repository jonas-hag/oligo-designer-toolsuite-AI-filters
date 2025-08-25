#!/bin/bash

#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --mem=120G
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --ntasks=1

#SBATCH --job-name=dataset_bowtie_test
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

source /home/haicu/jonas.hagenberg/.bashrc
conda activate odt-ai
python /lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/oligo_designer_toolsuite_ai_filters/hybridization_probability/generate_real_dataset.py -c /lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/configs/hybridization_probability/generate_real_dataset_bowtie_test_4.yaml