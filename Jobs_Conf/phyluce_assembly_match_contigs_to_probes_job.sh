#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -c 1


module load miniconda/3-2023-09

source $(conda info --base)/etc/profile.d/conda.sh
source activate /home/tiagobelintani/miniconda3/envs/phyluce-1.7.3

phyluce_assembly_match_contigs_to_probes \
    --contigs /home/tiagobelintani/uce-treinamento/assembly/contigs \
    --probes /home/tiagobelintani/uce-treinamento/probes/probes.fasta \
    --min-coverage 80 \
    --min-identity 80 \
     --output uce-resultados-busca \
    --log-path log 


