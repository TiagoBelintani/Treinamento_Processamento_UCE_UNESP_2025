#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G

# Ativar ambiente corretamente
module load miniconda/3-2023-09
source activate /home/tiagobelintani/miniconda3/envs/phyluce-1.7.3

# Diretórios
BASE_DIR="$HOME/uce-treinamento/raw-fastq"
OUT_BASE_DIR="$HOME/uce-treinamento/clean-fastq"

mkdir -p "$OUT_BASE_DIR"

# Loop pelos arquivos R1
for r1 in "$BASE_DIR"/*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    r2="${BASE_DIR}/${sample}_R2.fastq.gz"

    echo "🔍 Verificando: $sample"
    echo "  R1: $r1"
    echo "  R2: $r2"

        if [[ -f "$r2" ]]; then
        output_dir="${OUT_BASE_DIR}/${sample}/split-adapter-quality-trimmed"
        mkdir -p "$output_dir"

        echo "✂️ Rodando Trim Galore para $sample..."
        trim_galore \
            --paired "$r1" "$r2" \
            --cores 12 \
            --output_dir "$output_dir" \
            --gzip
    else
        echo "⚠️ Arquivo R2 não encontrado para $sample — pulando."
    fi
done

echo "✅ Trim Galore finalizado com sucesso."
