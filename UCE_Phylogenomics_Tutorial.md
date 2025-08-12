
# UCE Phylogenomics: Tutorial de Processamento de Dados de Aranhas Mygalomorphae

<div align="justify">
Este tutorial visa descrever, de forma metodológica e reprodutível, o processamento de dados de enriquecimento de UCEs (Ultra-Conserved Elements) provenientes de amostras de aranhas da subordem Mygalomorphae.
</div>

## Referência Bibliográfica

Para detalhes adicionais sobre a extração, construção das bibliotecas e sequenciamento dos dados, consultar o artigo recentemente publicado:

**DOI:** [https://doi.org/10.1016/j.ympev.2025.108323](https://doi.org/10.1016/j.ympev.2025.108323)

## Acesso aos Dados

Os dados utilizados neste tutorial encontram-se depositados no repositório público do NCBI:

**BioProject:** PRJNA1161786  
**Link:** [SRA GenBank](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP561602)

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

<div align="justify">
A limpeza de dados de sequenciamento de nova geração (NGS) é um passo crucial para garantir a qualidade e a confiabilidade das análises subsequentes. Ferramentas como o [Trim Galore](https://github.com/FelixKrueger/TrimGalore), [Trimmomatic](https://github.com/timflutre/trimmomatic) atuam removendo sequências adaptadoras e filtrando leituras de baixa qualidade, que podem introduzir ruído ou enviesar resultados. Durante o processo de sequenciamento, é comum que resíduos técnicos, como adaptadores não removidos ou bases com qualidade deteriorada nas extremidades — se acumulem nas leituras. Esses artefatos, se não tratados, podem levar a alinhamentos incorretos, montagem de genomas incompleta e interpretações equivocadas dos dados biológicos.
</div>
  
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

<div align="justify">
No fluxo de análise de dados no PHYLUCE, a etapa de montagem é responsável por reconstruir sequências contíguas (contigs) a partir das leituras limpas de sequenciamento. 
Para isso, o PHYLUCE oferece suporte a diferentes programas de montagem, todos integrados por meio de scripts próprios, mantendo um padrão de entrada e saída.

Os montadores disponíveis no [PHYLUCE] (https://phyluce.readthedocs.io/en/latest/) incluem:

[SPAdes](https://github.com/ablab/spades): geralmente a opção recomendada. É fácil de instalar, produz resultados consistentes e costuma apresentar melhor desempenho na maioria dos conjuntos de dados processados com PHYLUCE.

[Velvet](https://github.com/dzerbino/velvet): indicado para montagens de genomas menores ou dados com boa cobertura, sendo eficiente e rápido em cenários menos complexos.

[ABySS](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411771/): voltado para conjuntos de dados maiores ou genomas mais complexos, capaz de lidar com grandes volumes de leituras.

Para utilizar o SPAdes dentro do PHYLUCE, o comando típico é:
</div>

```bash
phyluce_assembly_assemblo_spades \ 
  --conf assembly.conf \ %+
  --output spades-assemblies \
  --cores 12 \
  --memory 64
```

Explicação dos parâmetros:

--conf assembly.conf → Arquivo de configuração que lista as amostras, os caminhos para os arquivos FASTQ e parâmetros opcionais de montagem.

--output spades-assemblies → Pasta onde os resultados da montagem serão salvos. Cada amostra terá seu próprio diretório com os contigs.

--cores 12 → Número de núcleos de CPU a serem usados, acelerando o processamento.

--memory 64 → Quantidade de memória RAM (em GB) disponível para a execução do SPAdes.

Após a execução, a pasta de saída conterá os contigs prontos para as próximas etapas, como identificação e extração dos loci alvo.


# Passos Práticos

## 1. Preparar o Arquivo de Configuração

O arquivo `assembly.conf` é essencial para que o **PHYLUCE** saiba onde encontrar os arquivos de leituras já limpas de cada amostra.  

Ele deve conter uma lista com o nome da amostra e o caminho para o diretório `split-adapter-quality-trimmed` correspondente.

Pode ser facilmente construido usando um editor de tabelas (excel ou outro pacote).

Exemplo de listagem simples:

```bash
[samples]
Arbanitis_rapax:/home/tiagobelintani/uce-treinamento/clean-fastq/Arbanitis_rapax/split-adapter-quality-trimmed
```

Atenção:

Os nomes das amostras não devem conter espaços.

É importante que o caminho seja exato e que o diretório exista.

Use nomes consistentes (iguais aos usados nos arquivos de leitura) para evitar erros.

Se os caminhos tiverem espaços antes ou depois dos dois-pontos (:), remova-os usando sed:


Agora um comando utilizando "sed" para retirar espaços indesejáveis 

```bash
sed -E 's/[[:space:]]*:[[:space:]]*/:/g' tabela.txt > assembly.conf
```

#Modelo de assembly.conf para várias amostras:

*Este arquivo é propenso a erros; verifique cuidadosamente para evitar problemas causados por espaços extras, caminhos inválidos ou erros de digitação.*

[Arquivo modelo - assembly.conf](https://github.com/TiagoBelintani/Treinamento_Processamento_UCE_UNESP_2025/tree/main/Jobs_Conf) *precisa de edição


```bash
[samples]
Arbanitis_rapax:/home/tiagobelintani/uce-treinamento/clean-fastq/Arbanitis_rapax/split-adapter-quality-trimmed
Cteniza_sp.:/home/tiagobelintani/uce-treinamento/clean-fastq/Cteniza_sp./split-adapter-quality-trimmed
Ctenolophus_sp.:/home/tiagobelintani/uce-treinamento/clean-fastq/Ctenolophus_sp./split-adapter-quality-trimmed
Galeosoma_sp.:/home/tiagobelintani/uce-treinamento/clean-fastq/Galeosoma_sp./split-adapter-quality-trimmed
Gorgyrella_namaquensis:/home/tiagobelintani/uce-treinamento/clean-fastq/Gorgyrella_namaquensis/split-adapter-quality-trimmed
...
Titanidiops_sp.:/home/tiagobelintani/uce-treinamento/clean-fastq/Titanidiops_sp./split-adapter-quality-trimmed
```

2. Criar o Script de Montagem
   
#Crie um arquivo chamado spades_job.sh com o seguinte conteúdo:

```bash
#!/bin/bash
#SBATCH -t 30:00:00           # Tempo máximo de execução (30h)
#SBATCH -c 12                 # Número de CPUs
#SBATCH --mem=128             # Memória total disponível (GB)

# Carregar miniconda e o ambiente do PHYLUCE
module load miniconda/3-2023-09
source $(conda info --base)/etc/profile.d/conda.sh
source activate /home/seu_nome/miniconda3/envs/phyluce-1.7.3  # Ajuste para o caminho do seu ambiente

# Executar o SPAdes via PHYLUCE
phyluce_assembly_assemblo_spades \
  --output assembly \         # Pasta de saída
  --cores 12 \                 # Número de núcleos
  --memory 64 \                 # Memória para o SPAdes (GB)
  --log-path log \              # Pasta para arquivos de log
  --config assembly.conf        # Arquivo de configuração
```

#3. Submeter e Monitorar o Job
Envie o job para execução no SLURM:

```bash
sbatch spades_job.sh
Acompanhe o status:
```

```bash
squeue -u tiagobelintani
```

Possíveis Problemas e Como Evitar

| Problema comum                       | Possível causa                                    | Solução                                                               |
| ------------------------------------ | ------------------------------------------------- | --------------------------------------------------------------------- |
| **Erro "No such file or directory"** | Caminho no `assembly.conf` incorreto              | Verificar caminhos e corrigir.                                        |
| **Script falha logo no início**      | Ambiente conda não ativado corretamente           | Confirmar o caminho do ambiente `phyluce` no script.                  |
| **Montagem muito lenta**             | Poucos núcleos ou memória insuficiente            | Ajustar `--cores` e `--memory` conforme a disponibilidade do cluster. |
| **Diretório de saída incompleto**    | Execução interrompida ou falta de espaço em disco | Checar logs em `log/` e espaço disponível.                            |


#Estrutura de diretorios atuais 

*não é obrigatorio tal estrutura, mas pode facilitar muito a organização do fluxo de trabalho*

```bash
/home/tiagobelintani/uce-treinamento/
├── assembly
│   ├── Arbanitis_rapax_spades
│   ├── contigs
│   ├── ... 
├── clean-fastq
│   ├── Arbanitis_rapax
│   ├── ...
├── config
│   └── assembly.conf
├── job
├── lista.txt
├── log
└── raw-fastq
```

#Encontrar UCE loci (Fingding UCE loci)


<div align="justify">
Após a montagem das leituras em contigs, o próximo passo no PHYLUCE é identificar quais desses contigs contêm loci UCE (Ultra-Conserved Elements).
Esse processo é importante porque, embora a montagem contenha todas as sequências resultantes do sequenciamento, apenas uma parte delas corresponde aos loci-alvo definidos pela sonda de captura utilizada no experimento.
</div>



<div align="justify">
O PHYLUCE realiza essa identificação comparando os contigs montados com um banco de dados de loci UCEs de referência, geralmente fornecido em formato .fasta. Essa comparação é feita usando algoritmos de alinhamento rápido, como lastz, que detectam regiões de alta similaridade.
</div>



#Passos para encontrar loci UCE no PHYLUCE


<div align="justify">
Organizar o diretório de montagem
Certifique-se de que todas as pastas de montagem (por amostra) estão reunidas em um único diretório.
Cada pasta deve conter o arquivo contigs.fasta gerado pelo montador.
</div>

#Preparar o banco de sondas UCE

Baixe ou utilize o conjunto de sondas específico para o seu grupo de estudo (por exemplo, insetos, aves, aracnídeos).
Esse arquivo .fasta será usado como referência.

Executar o alinhamento com lastz
O comando típico no PHYLUCE para essa etapa é:

```bash
phyluce_assembly_match_contigs_to_probes \
  --contigs /caminho/para/assemblies \ #diretorio com os contigs (symlinks)
  --probes uce-probes.fasta \ #diretorio com as probes
  --output /caminho/para/uces_finder_output #arquivos tabulados com squile3
```


Gerar a tabela de loci encontrados

O resultado será um conjunto de arquivos que listam quais loci foram encontrados em cada amostra.
Esses arquivos são usados nas etapas seguintes de filtragem e extração.

Possíveis problemas

```bash
| Problema                          | Causa comum                                 | Como corrigir                                                                       |
| --------------------------------- | ------------------------------------------- | ----------------------------------------------------------------------------------- |
| Nenhum locus encontrado           | Sonda incompatível com seu grupo taxonômico | Verificar se está usando o conjunto de sondas correto.                              |
| Alinhamento muito lento           | Muitas amostras e/ou sonda muito grande     | Usar mais *cores* (`--cores`) e otimizar o cluster.                                 |
| Arquivo de contigs não encontrado | Estrutura de diretórios incorreta           | Conferir se o caminho passado em `--contigs` está correto e contém `contigs.fasta`. |
```


# Passos Práticos

## Baixar as probes (**iscas**)


#criar dir para receber as probes

```bash
mkdir probes
```
#Opções de downloud

[Probes] (https://figshare.com/articles/dataset/Arachnida_14_799_baits_targeting_1_120_conserved_loci_Arachnida-UCE-1_1K-v1_fasta_/3856146)

[Probes] (https://drive.google.com/file/d/1BTGtLKJQHw1uxE7X2kTqiSm0Gmnt8KwM/view?usp=drive_link)

[Probes] (https://github.com/TiagoBelintani/Treinamento_Processamento_UCE_UNESP_2025/tree/main/Probes)

## 2. Preparar o Arquivo de Configuração



<div align="justify">
Agora, no diretório inicial das suas análises, você deverá criar um job e executá-lo no terminal.
</div>

```bash
pwd
/home/seu_nome/uce-treinamento/
```
#criar o arquivo config

```bash
nano phyluce_assembly_match_contigs_to_probes_job.sh
```
#(arquivo config)[

```bash
#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -c 1


module load miniconda/3-2023-09

source $(conda info --base)/etc/profile.d/conda.sh
source activate /home/tiagobelintani/miniconda3/envs/phyluce-1.7.3

phyluce_assembly_match_contigs_to_probes \
    --contigs /home/tiagobelintani/uce-treinamento/assembly/contigs \
    --probes uce-5k-probes.fasta \
    --output uce-search-results
```













