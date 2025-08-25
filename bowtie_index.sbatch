#!/bin/bash

#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --mem=120G
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1

#SBATCH --job-name=bowtie_index
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

source /home/haicu/jonas.hagenberg/.bashrc
conda activate odt-ai
python /home/haicu/jonas.hagenberg/projects/odt-ai/02_code/generate_bowtie_index.py -c /home/haicu/jonas.hagenberg/projects/odt-ai/03_config/generate_bowtie_index.yaml