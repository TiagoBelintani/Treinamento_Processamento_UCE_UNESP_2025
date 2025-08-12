#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -c 1


module load miniconda/3-2023-09

source $(conda info --base)/etc/profile.d/conda.sh
source activate /home/tiagobelintani/miniconda3/envs/phyluce-1.7.3


phyluce_assembly_get_fastas_from_match_counts \
    --contigs /home/tiagobelintani/uce-treinamento/assembly_2/contigs/ \
    --locus-db /home/tiagobelintani/uce-treinamento/uce-resultados-busca/probe.matches.sqlite \
    --match-count-output /home/tiagobelintani/uce-treinamento/taxon-set/all/all-taxa-incomplete.conf \
    --output  all-taxa-incomplete.fasta \
    --incomplete-matrix /home/tiagobelintani/uce-treinamento/taxon-set/all/all-taxa-incomplete.incomplete \
    --log-path log
