#!/bin/bash
#SBATCH -t 30:00:00
#SBATCH -c 12
#SBATCH --mem=128

module load miniconda/3-2023-09

source $(conda info --base)/etc/profile.d/conda.sh
source activate /home/***[seu_nome]***/miniconda3/envs/phyluce-1.7.3

phyluce_assembly_assemblo_spades \
  --output assembly_2 \
  --cores 12 \
  --memory 128  --log-path log \
  --config assembly.conf
