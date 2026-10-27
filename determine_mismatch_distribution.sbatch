#!/bin/bash

#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --ntasks=1

#SBATCH --job-name=mismatch_distribution
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

source /home/haicu/jonas.hagenberg/.bashrc
conda activate odt-ai-2

# stop the job if one of the regions fails, i.e. don't summarise incomplete results
set -e

SCRIPT_DIR=/lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/oligo_designer_toolsuite_ai_filters/hybridization_probability
ALIGNMENT_DIR=/lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/data/hybridization_probability/real_dataset_sampled_oligos/blast/selected_oligos
MISMATCH_DIR=/lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/data/hybridization_probability/real_dataset_sampled_oligos/blast/selected_oligos/mismatches
MISMATCH_OUT_DIR=/lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/data/hybridization_probability/real_dataset_sampled_oligos/blast/selected_oligos/mismatch_distribution

mkdir -p $MISMATCH_DIR

# one process per region file, $SLURM_CPUS_PER_TASK files at a time
ls $ALIGNMENT_DIR/*.csv | xargs -P $SLURM_CPUS_PER_TASK -I {} python $SCRIPT_DIR/determine_mismatch_distribution.py -i {} -o $MISMATCH_DIR

python $SCRIPT_DIR/summarise_mismatch_distribution.py -i $MISMATCH_DIR -o $MISMATCH_OUT_DIR
