
# UCE Phylogenomics: Tutorial de Processamento de Dados de Aranhas Mygalomorphae

Este tutorial visa descrever, de forma metodológica e reprodutível, o processamento de dados de enriquecimento de UCEs (Ultra-Conserved Elements) provenientes de amostras de aranhas da subordem Mygalomorphae.

## Referência Bibliográfica

Para detalhes adicionais sobre a extração, construção das bibliotecas e sequenciamento dos dados, consultar o artigo recentemente publicado:

**DOI:** [https://doi.org/10.1016/j.ympev.2025.108323](https://doi.org/10.1016/j.ympev.2025.108323)

## Acesso aos Dados

Os dados utilizados neste tutorial encontram-se depositados no repositório público do NCBI:

**BioProject:** PRJNA1161786  
**Link:** [https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP561602](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP561602)

### Conjunto de Acessions

O conjunto de acessions completo está disponível na URL acima, sendo composto por 24 amostras, cobrindo uma ampla diversidade geográfica e taxonômica.

---

## Estruturação do Ambiente de Trabalho

Inicialmente, é necessário criar a seguinte estrutura de diretórios:

```bash
mkdir -p uce-treinamento/raw-fastq uce-treinamento/clean-fastq uce-treinamento/log
chmod -R 777 uce-treinamento/
```

Verifique a criação correta com:

```bash
ls
```

---

## Transferência dos Dados

Acesse o diretório `raw-fastq`:

```bash
cd uce-treinamento/raw-fastq
```

Transfira os dados de leitura emparelhada (`*_R1.fastq.gz` e `*_R2.fastq.gz`) do seguinte repositório:

**Link:** [https://drive.google.com/drive/folders/10BTIOdb93iEyvOi_YnrTJFKo2wRVNGtS](https://drive.google.com/drive/folders/10BTIOdb93iEyvOi_YnrTJFKo2wRVNGtS)

**Link:** [git clone https://github.com/TiagoBelintani/Treinamento_Processamento_UCE_UNESP_2025.git](https://github.com/TiagoBelintani/Treinamento_Processamento_UCE_UNESP_2025.git)

---

## Verificação de Integridade e Quantificação

Para contar o número de reads em cada arquivo `.fastq.gz`:

```bash
for i in *_R1.fastq.gz; do
  echo "$i"
  gunzip -c "$i" | wc -l | awk '{print $1/4}'
done
```

---

## Renomeação dos Arquivos

### 1. Criar Tabela de Mapeamento

Crie um arquivo `rename_map.tsv` contendo o mapeamento entre os acessions e os nomes dos táxons:

```bash
nano rename_map.tsv
```

Insira as linhas com os pares `SRRxxx<TAB>Genus_species`.

### 2. Script de Renomeação

Salve o seguinte script como `rename_fastq.sh`:

```bash
#!/bin/bash

MAP_FILE="rename_map.tsv"

while IFS=$'	' read -r srr species; do
    for pair in R1 R2; do
        old="${srr}_${pair}.fastq.gz"
        new="${species}_${pair}.fastq.gz"

        if [[ -f "$old" ]]; then
            echo "Renomeando $old → $new"
            mv "$old" "$new"
        else
            echo "Arquivo $old não encontrado — pulando."
        fi
    done
done < "$MAP_FILE"
```

Execute o script:

```bash
bash rename_fastq.sh
```

---

## Limpeza dos Dados com Trim Galore

A limpeza de dados de sequenciamento de nova geração (NGS) é um passo crucial para garantir a qualidade e a confiabilidade das análises subsequentes. Ferramentas como o (Trim Galore)[https://github.com/FelixKrueger/TrimGalore], (Trimmomatic) [https://github.com/timflutre/trimmomatic] atuam removendo sequências adaptadoras e filtrando leituras de baixa qualidade, que podem introduzir ruído ou enviesar resultados. Durante o processo de sequenciamento, é comum que resíduos técnicos — como adaptadores não removidos ou bases com qualidade deteriorada nas extremidades — se acumulem nas leituras. Esses artefatos, se não tratados, podem levar a alinhamentos incorretos, montagem de genomas incompleta e interpretações equivocadas dos dados biológicos.

### 1. Ativar o Ambiente

```bash
conda activate phyluce-1.7.3
```

Instalar o Trim Galore:

```bash
conda install trim-galore
```

Verificar a instalação:

```bash
trim_galore --version
```

### 2. Criar o Script de Trimming

Salve como `trimgalore.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G

module load miniconda/3-2023-09
source activate /home/tiagobelintani/miniconda3/envs/phyluce-1.7.3

BASE_DIR="$HOME/uce-treinamento/raw-fastq"
OUT_BASE_DIR="$HOME/uce-treinamento/clean-fastq"

mkdir -p "$OUT_BASE_DIR"

for r1 in "$BASE_DIR"/*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    r2="${BASE_DIR}/${sample}_R2.fastq.gz"

    if [[ -f "$r2" ]]; then
        output_dir="${OUT_BASE_DIR}/${sample}/split-adapter-quality-trimmed"
        mkdir -p "$output_dir"

        trim_galore             --paired "$r1" "$r2"             --cores 12             --output_dir "$output_dir"             --gzip
    else
        echo "Arquivo R2 não encontrado para $sample — pulando."
    fi
done

echo "Trim Galore finalizado com sucesso."
```

---

## Montagem dos Dados com SPAdes

### 1. Preparar o Arquivo de Configuração

Utilize a listagem dos diretórios `split-adapter-quality-trimmed` para montar o arquivo `assembly.conf`.

Exemplo do conteúdo:

```ini
[samples]
Arbanitis_rapax:/home/tiagobelintani/uce-treinamento/clean-fastq/Arbanitis_rapax/split-adapter-quality-trimmed
...
```

Remova espaços desnecessários com `sed`:

```bash
sed -E 's/[[:space:]]*:[[:space:]]*/:/g' tabela.txt > assembly.conf
```

Assembly.conf model

```bash
[samples]
Arbanitis_rapax:/home/tiagobelintani/uce-treinamento/clean-fastq/Arbanitis_rapax/split-adapter-quality-trimmed
Cteniza_sp.:/home/tiagobelintani/uce-treinamento/clean-fastq/Cteniza_sp./split-adapter-quality-trimmed
Ctenolophus_sp.:/home/tiagobelintani/uce-treinamento/clean-fastq/Ctenolophus_sp./split-adapter-quality-trimmed
Galeosoma_sp.:/home/tiagobelintani/uce-treinamento/clean-fastq/Galeosoma_sp./split-adapter-quality-trimmed
Gorgyrella_namaquensis:/home/tiagobelintani/uce-treinamento/clean-fastq/Gorgyrella_namaquensis/split-adapter-quality-trimm$
Heligmomerus_sp.:/home/tiagobelintani/uce-treinamento/clean-fastq/Heligmomerus_sp./split-adapter-quality-trimmed
Idiops_camelus:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_camelus/split-adapter-quality-trimmed
Idiops_carajas:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_carajas/split-adapter-quality-trimmed
Idiops_clarus:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_clarus/split-adapter-quality-trimmed
Idiops_fryi:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_fryi/split-adapter-quality-trimmed
Idiops_germaini:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_germaini/split-adapter-quality-trimmed
Idiops_guri:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_guri/split-adapter-quality-trimmed
Idiops_kanonganus:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_kanonganus/split-adapter-quality-trimmed
Idiops_petiti:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_petiti/split-adapter-quality-trimmed
Idiops_pirassununguensis:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_pirassununguensis/split-adapter-quality-t$
Idiops_pretoriae:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_pretoriae/split-adapter-quality-trimmed
Idiops_rastratus:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_rastratus/split-adapter-quality-trimmed
Idiops_rohdei:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_rohdei/split-adapter-quality-trimmed
Idiops_sp2_RF2025:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_sp2_RF2025/split-adapter-quality-trimmed
Idiops_sp3_RF2025:/home/tiagobelintani/uce-treinamento/clean-fastq/Idiops_sp3_RF2025/split-adapter-quality-trimmed
Moggridgea_crudeni:/home/tiagobelintani/uce-treinamento/clean-fastq/Moggridgea_crudeni/split-adapter-quality-trimmed
Neocteniza_toba:/home/tiagobelintani/uce-treinamento/clean-fastq/Neocteniza_toba/split-adapter-quality-trimmed
Segregara_transvaalensis:/home/tiagobelintani/uce-treinamento/clean-fastq/Segregara_transvaalensis/split-adapter-quality-t$
Titanidiops_sp.:/home/tiagobelintani/uce-treinamento/clean-fastq/Titanidiops_sp./split-adapter-quality-trimmed
```


### 2. Script de Montagem

Salve como `spades_job.sh`:

```bash
#!/bin/bash
#SBATCH -t 30:00:00
#SBATCH -c 12
#SBATCH --mem=128

module load miniconda/3-2023-09

source $(conda info --base)/etc/profile.d/conda.sh
source activate /home/seu_nome/miniconda3/envs/phyluce-1.7.3 #alterar o caminho para seu ambiente

phyluce_assembly_assemblo_spades \
  --output assembly \
  --cores 12 \
  --memory 64 \
  --log-path log \
  --config assembly.conf
```

Submeta o job via SLURM:

```bash
sbatch spades_job.sh
```

Verifique o status dos jobs:

```bash
squeue -u tiagobelintani
```

---
## Encontrando os Locis UCE (Finding UCE)
















