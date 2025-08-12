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
A limpeza de dados de sequenciamento de nova geração (NGS) é um passo crucial para garantir a qualidade e a confiabilidade das análises subsequentes. Ferramentas como o <a href="https://github.com/FelixKrueger/TrimGalore">Trim Galore</a> e o <a href="https://github.com/timflutre/trimmomatic">Trimmomatic</a> atuam removendo sequências adaptadoras e filtrando leituras de baixa qualidade, que podem introduzir ruído ou enviesar resultados. Durante o processo de sequenciamento, é comum que resíduos técnicos, como adaptadores não removidos ou bases com qualidade deteriorada nas extremidades, se acumulem nas leituras. Esses artefatos, se não tratados, podem levar a alinhamentos incorretos, montagem de genomas incompleta e interpretações equivocadas dos dados biológicos.
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

Os montadores disponíveis no <a href="https://phyluce.readthedocs.io/en/latest/">Phyluce</a> incluem:

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

**Explicação dos parâmetros:**

```bash
--conf assembly.conf → Arquivo de configuração que lista as amostras, os caminhos para os arquivos FASTQ e parâmetros opcionais de montagem.

--output spades-assemblies → Pasta onde os resultados da montagem serão salvos. Cada amostra terá seu próprio diretório com os contigs.

--cores 12 → Número de núcleos de CPU a serem usados, acelerando o processamento.

--memory 64 → Quantidade de memória RAM (em GB) disponível para a execução do SPAdes.

Após a execução, a pasta de saída conterá os contigs prontos para as próximas etapas, como identificação e extração dos loci alvo.
```

# Passos Práticos

## 1. Preparar o Arquivo de Configuração

<div align="justify"> 
O arquivo `assembly.conf` é essencial para que o <a href="https://phyluce.readthedocs.io/en/latest/">Phyluce</a> saiba onde encontrar os arquivos de leituras já limpas de cada amostra.  

Ele deve conter uma lista com o nome da amostra e o caminho para o diretório `split-adapter-quality-trimmed` correspondente.

Pode ser facilmente construido usando um editor de tabelas (excel ou outro pacote).
</div>

*Exemplo de listagem simples:*

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

#Resultados 

<div align="justify">  
Nesta prática, duas montagens ficaram incompletas devido à baixa qualidade inicial dos dados de sequenciamento. Nesse caso específico, observou-se um número reduzido de leituras e uma proporção significativa de sequências de baixa qualidade.
</p>
  
<div align="justify"> 
Esse problema poderia ser, em parte, solucionado por meio de novas tentativas utilizando outros montadores (assemblers) ou ajustando parâmetros críticos, como diferentes valores de k-mer. Ainda assim, é importante ressaltar que dados de qualidade insatisfatória impõem limitações severas às análises filogenômicas, comprometendo a recuperação de loci, a completude das montagens e, consequentemente, a robustez das inferências evolutivas.
</p>
  
<div align="justify">   
Portanto, investir em uma etapa de sequenciamento bem planejada — garantindo cobertura adequada, qualidade de leitura elevada e estratégias de limpeza eficientes — é essencial para minimizar perdas de informação e maximizar o sucesso das etapas subsequentes de montagem e análise.
#Estrutura de diretorios atuais 
</p>
  
*não é obrigatorio tal estrutura, mas pode facilitar muito a organização do fluxo de trabalho*

#Vamos continuar verificando a estrutura dos diretórios.

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
---

# Encontrar UCE loci (Finding UCE loci)

<p align="justify">
Encontrando os Loci UCE (Finding UCE)
Após a montagem das leituras em contigs, o próximo passo no PHYLUCE é identificar quais desses contigs contêm loci UCE (Ultra-Conserved Elements).
Esse processo é importante porque, embora a montagem contenha todas as sequências resultantes do sequenciamento, apenas uma parte delas corresponde aos loci-alvo definidos pela sonda de captura utilizada no experimento.
</p>

<p align="justify">
O PHYLUCE realiza essa identificação comparando os contigs montados com um banco de dados de loci UCEs de referência, geralmente fornecido em formato .fasta. Essa comparação é feita usando algoritmos de alinhamento rápido, como lastz, que detectam regiões de alta similaridade.
</p>

#**Passos para encontrar loci UCE no <a href="https://phyluce.readthedocs.io/en/latest/">Phyluce</a>**

<p align="justify">
#Organizar o diretório de montagem
Certifique-se de que todas as pastas de montagem (por amostra) estão reunidas em um único diretório.
Cada pasta deve conter o arquivo contigs.fasta gerado pelo montador.
</p>

#Preparar o banco de sondas UCE
 
Baixe ou utilize o conjunto de sondas específico para o seu grupo de estudo (por exemplo, insetos, aves, aracnídeos).
Esse arquivo .fasta será usado como referência.
 
#Executar o alinhamento com lastz

O comando típico no PHYLUCE para essa etapa é:

```bash
phyluce_assembly_match_contigs_to_probes \
  --contigs /caminho/para/assemblies \ #diretorio com os contigs (symlinks)
  --probes uce-probes.fasta \ #diretorio com as probes
  --output /caminho/para/uces_finder_output #arquivos tabulados com squile3
```


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

----

# Passos Práticos

## Baixar as probes (**iscas**)


#criar dir para receber as probes

```bash
mkdir probes
```
#Opções de downloud

[Probes – Figshare](https://figshare.com/articles/dataset/Arachnida_14_799_baits_targeting_1_120_conserved_loci_Arachnida-UCE-1_1K-v1_fasta_/3856146)  
[Probes – Google Drive](https://drive.google.com/file/d/1BTGtLKJQHw1uxE7X2kTqiSm0Gmnt8KwM/view?usp=drive_link)  
[Probes – GitHub](https://github.com/TiagoBelintani/Treinamento_Processamento_UCE_UNESP_2025/tree/main/Probes)

#renomear probes

```bash
mv * probes.fasta
```
*Estrutura diretorio atual

```bash
/home/tiagobelintani/uce-treinamento/probes/
└── probes.fasta
```

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
#[Arquivo de execução - job slurm](https://github.com/TiagoBelintani/Treinamento_Processamento_UCE_UNESP_2025/blob/main/Jobs_Conf/phyluce_assembly_match_contigs_to_probes_job.sh)

```bash
#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -c 1


module load miniconda/3-2023-09

source $(conda info --base)/etc/profile.d/conda.sh
source activate /home/tiagobelintani/miniconda3/envs/phyluce-1.7.3

phyluce_assembly_match_contigs_to_probes \
    --contigs /home/tiagobelintani/uce-treinamento/assembly/contigs \
    --probes /home/tiagobelintani/uce-treinamento/probes/probes.fasta \
    --output uce-resultados-busca
```

#Deverá ver um resultado semelhante ao seguinte [também armazenado em log](https://github.com/TiagoBelintani/Treinamento_Processamento_UCE_UNESP_2025/blob/main/LOGS/phyluce_assembly_match_contigs_to_probes.log)

```bash
2025-08-12 16:45:05,072 - phyluce_assembly_match_contigs_to_probes - INFO - ======= Starting phyluce_assembly_match_contigs_to_probes =======
2025-08-12 16:45:05,072 - phyluce_assembly_match_contigs_to_probes - INFO - Version: 1.7.3
2025-08-12 16:45:05,072 - phyluce_assembly_match_contigs_to_probes - INFO - Commit: None
2025-08-12 16:45:05,072 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --contigs: /home/tiagobelintani/uce-treinamento/assembly_2/contigs
2025-08-12 16:45:05,072 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --csv: None
2025-08-12 16:45:05,072 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --dupefile: None
2025-08-12 16:45:05,072 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --keep_duplicates: None
2025-08-12 16:45:05,072 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --log_path: None
2025-08-12 16:45:05,072 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --min_coverage: 80
2025-08-12 16:45:05,072 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --min_identity: 80
2025-08-12 16:45:05,072 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --output: /home/tiagobelintani/uce-treinamento/uce-resultados-busca
2025-08-12 16:45:05,072 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --probes: /home/tiagobelintani/uce-treinamento/probes/probes.fasta
2025-08-12 16:45:05,073 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --regex: ^(uce-\d+)(?:_p\d+.*)
2025-08-12 16:45:05,073 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --verbosity: INFO
2025-08-12 16:45:05,142 - phyluce_assembly_match_contigs_to_probes - INFO - Creating the UCE-match database
2025-08-12 16:45:05,733 - phyluce_assembly_match_contigs_to_probes - INFO - Processing contig data
2025-08-12 16:45:05,734 - phyluce_assembly_match_contigs_to_probes - INFO - -----------------------------------------------------------------
2025-08-12 16:45:08,138 - phyluce_assembly_match_contigs_to_probes - INFO - Cteniza_sp: 120 (4.93%) uniques of 2434 contigs, 0 dupe probe matches, 5 UCE loci removed for matching multiple contigs, 18 contigs removed for matching multiple UCE loci
2025-08-12 16:45:10,460 - phyluce_assembly_match_contigs_to_probes - INFO - Ctenolophus_sp: 110 (4.58%) uniques of 2401 contigs, 0 dupe probe matches, 3 UCE loci removed for matching multiple contigs, 14 contigs removed for matching multiple UCE loci
2025-08-12 16:45:11,597 - phyluce_assembly_match_contigs_to_probes - INFO - Gorgyrella_namaquensis: 41 (3.31%) uniques of 1238 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 9 contigs removed for matching multiple UCE loci
2025-08-12 16:45:12,827 - phyluce_assembly_match_contigs_to_probes - INFO - Heligmomerus_sp: 51 (3.69%) uniques of 1383 contigs, 0 dupe probe matches, 1 UCE loci removed for matching multiple contigs, 8 contigs removed for matching multiple UCE loci
2025-08-12 16:45:13,605 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_camelus: 0 (0.00%) uniques of 1033 contigs, 0 dupe probe matches, 1 UCE loci removed for matching multiple contigs, 0 contigs removed for matching multiple UCE loci
2025-08-12 16:45:14,456 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_carajas: 4 (0.10%) uniques of 3860 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 1 contigs removed for matching multiple UCE loci
2025-08-12 16:45:15,279 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_clarus: 4 (0.08%) uniques of 4939 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 0 contigs removed for matching multiple UCE loci
2025-08-12 16:45:16,601 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_fryi: 57 (2.88%) uniques of 1978 contigs, 0 dupe probe matches, 1 UCE loci removed for matching multiple contigs, 6 contigs removed for matching multiple UCE loci
2025-08-12 16:45:17,242 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_germaini: 6 (0.39%) uniques of 1534 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 1 contigs removed for matching multiple UCE loci
2025-08-12 16:45:18,351 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_guri: 6 (0.10%) uniques of 5991 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 0 contigs removed for matching multiple UCE loci
2025-08-12 16:45:19,178 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_kanonganus: 5 (0.09%) uniques of 5415 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 0 contigs removed for matching multiple UCE loci
2025-08-12 16:45:20,228 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_petiti: 31 (0.94%) uniques of 3303 contigs, 0 dupe probe matches, 1 UCE loci removed for matching multiple contigs, 5 contigs removed for matching multiple UCE loci
2025-08-12 16:45:21,477 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_pirassununguensis: 52 (2.20%) uniques of 2361 contigs, 0 dupe probe matches, 1 UCE loci removed for matching multiple contigs, 0 contigs removed for matching multiple UCE loci
2025-08-12 16:45:22,510 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_pretoriae: 21 (0.69%) uniques of 3039 contigs, 0 dupe probe matches, 1 UCE loci removed for matching multiple contigs, 3 contigs removed for matching multiple UCE loci
2025-08-12 16:45:23,093 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_rastratus: 1 (0.06%) uniques of 1807 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 0 contigs removed for matching multiple UCE loci
2025-08-12 16:45:24,344 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_rohdei: 48 (1.19%) uniques of 4050 contigs, 0 dupe probe matches, 1 UCE loci removed for matching multiple contigs, 6 contigs removed for matching multiple UCE loci
2025-08-12 16:45:25,442 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_sp2_RF2025: 47 (2.01%) uniques of 2334 contigs, 0 dupe probe matches, 2 UCE loci removed for matching multiple contigs, 3 contigs removed for matching multiple UCE loci
2025-08-12 16:45:26,431 - phyluce_assembly_match_contigs_to_probes - INFO - Idiops_sp3_RF2025: 8 (0.08%) uniques of 9556 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 1 contigs removed for matching multiple UCE loci
2025-08-12 16:45:27,453 - phyluce_assembly_match_contigs_to_probes - INFO - Moggridgea_crudeni: 4 (0.11%) uniques of 3555 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 0 contigs removed for matching multiple UCE loci
2025-08-12 16:45:28,899 - phyluce_assembly_match_contigs_to_probes - INFO - Neocteniza_toba: 38 (0.49%) uniques of 7777 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 2 contigs removed for matching multiple UCE loci
2025-08-12 16:45:30,227 - phyluce_assembly_match_contigs_to_probes - INFO - Segregara_transvaalensis: 56 (3.78%) uniques of 1480 contigs, 0 dupe probe matches, 1 UCE loci removed for matching multiple contigs, 11 contigs removed for matching multiple UCE loci
2025-08-12 16:45:32,202 - phyluce_assembly_match_contigs_to_probes - INFO - Titanidiops_sp: 91 (3.52%) uniques of 2587 contigs, 0 dupe probe matches, 1 UCE loci removed for matching multiple contigs, 16 contigs removed for matching multiple UCE loci
2025-08-12 16:45:32,202 - phyluce_assembly_match_contigs_to_probes - INFO - -----------------------------------------------------------------
2025-08-12 16:45:32,202 - phyluce_assembly_match_contigs_to_probes - INFO - The LASTZ alignments are in /home/tiagobelintani/uce-treinamento/uce-resultados-busca
2025-08-12 16:45:32,202 - phyluce_assembly_match_contigs_to_probes - INFO - The UCE match database is in /home/tiagobelintani/uce-treinamento/uce-resultados-busca/probe.matches.sqlite
2025-08-12 16:45:32,202 - phyluce_assembly_match_contigs_to_probes - INFO - ======= Completed phyluce_assembly_match_contigs_to_probes ======
```
 #Breve intepretação

<div align="justify"> 
Com os parâmetros definidos para --min_coverage 80 e --min_identity 80 (default), o comando phyluce_assembly_match_contigs_to_probes realizou uma busca relativamente restritiva, exigindo que pelo menos 80% da sonda fosse coberta e que a similaridade mínima entre contig e probe também fosse de 80%. Esses valores favorecem a recuperação de loci de maior qualidade, mas podem reduzir bastante o número de UCEs identificados, especialmente em amostras mais divergentes ou com montagens fragmentadas. Isso explica a grande variação de resultados observada no log, com algumas espécies recuperando acima de 100 loci (ex.: Cteniza_sp com 120) e outras com praticamente nenhum locus (ex.: Idiops_camelus com 0). Além disso, o relatório indica remoções de loci por mapearem em múltiplos contigs ou contigs associados a múltiplos loci, o que evidencia a necessidade de avaliar a integridade e o nível de duplicação das montagens para otimizar parâmetros em futuras análises.
</div>

<p align="center">
  <img src="Imagens/uce_matches_por_especie.png" alt="UCEs únicos por espécie" width="800">
</p>

```bash

O comando `phyluce_assembly_match_contigs_to_probes` identifica quais *contigs* das montagens contêm loci UCE (*Ultra-Conserved Elements*) ao compará-los com um conjunto de sondas (*probes*) de referência.  

Este utilitário possui funções refinadas que podem otimizar a análise, permitindo ajustes finos de parâmetros.  
A escolha de cada configuração depende de fatores como:

- Qualidade das montagens.
- Qualidade das leituras iniciais.
- Similaridade com as sondas utilizadas.
- Cobertura e identidade mínimas desejadas.
```
---
**Explicação dos parâmetros**


```bash
--contigs
Caminho para o diretório contendo as montagens (assemblies), onde cada pasta de amostra deve conter o arquivo contigs.fasta.
Exemplo: /home/user/projeto/assemblies/.

--probes
Arquivo .fasta contendo as sequências das sondas UCE usadas na captura. Deve corresponder ao seu grupo taxonômico.
Exemplo: uce-probes-insetos.fasta.

--output
Diretório onde os resultados serão salvos. O programa criará subpastas e arquivos com informações sobre loci encontrados.

--verbosity (opcional)
Define o nível de mensagens exibidas durante a execução:

INFO → Mostra todas as mensagens (padrão).

WARN → Mostra apenas avisos e erros.

CRITICAL → Mostra apenas mensagens críticas.

--log-path (opcional)
Caminho para salvar o arquivo de log. Útil para revisitar mensagens e parâmetros usados.

--min-coverage (opcional)
Cobertura mínima (em %) que um contig deve ter em relação à sonda para ser considerado um match válido.
Valores altos reduzem falsos positivos, mas podem descartar loci parciais.

--min-identity (opcional)
Percentual mínimo de identidade de bases entre contig e sonda.
Ajustar de acordo com a distância evolutiva do grupo estudado:

Grupos próximos → valores altos (90–95%).

Grupos distantes → valores mais baixos (70–80%).

--dupefile (opcional)
Arquivo onde serão registradas ocorrências de loci duplicados. Útil para inspeção posterior.

--regex (opcional)
Expressão regular para filtrar nomes de amostras ou contigs. Pode ser usada para processar subconjuntos específicos.

--keep-duplicates (opcional)
Se especificado, mantém loci duplicados nos resultados (por padrão, duplicatas são removidas).

--csv (opcional)
Gera saída no formato .csv, facilitando a importação para planilhas ou scripts de análise.
```

**Possíveis problemas e soluções**

| Problema                               | Possível causa                                                 | Como resolver                                                         |
| -------------------------------------- | -------------------------------------------------------------- | --------------------------------------------------------------------- |
| Nenhum locus encontrado                | Sondas incompatíveis com o grupo analisado                     | Use o conjunto de sondas correto para o táxon.                        |
| Poucos loci recuperados                | Parâmetros de `--min-coverage` ou `--min-identity` muito altos | Ajustar valores para permitir correspondências mais relaxadas.        |
| Execução lenta                         | Muitas amostras ou sondas muito grandes                        | Usar mais *cores* e otimizar recursos do cluster.                     |
| Arquivo `contigs.fasta` não encontrado | Estrutura de diretórios incorreta                              | Verificar caminhos em `--contigs` e garantir que os arquivos existam. |

---

### Extraindo Loci UCE (Extracting UCE loci)

Após a identificação dos loci UCE, o passo seguinte consiste em definir quais táxons serão incluídos na análise. Para isso, é necessário criar um arquivo contendo a lista desses táxons e, a partir dele, gerar um arquivo de configuração da matriz de dados. Esse arquivo de configuração especificará, para cada táxon, quais loci foram efetivamente enriquecidos, servindo como referência para a extração das respectivas sequências no formato FASTA.

---

### 1. Criar o conjunto de táxons (*taxon set*)

Decida quais táxons quer incluir. O nome de cada táxon deve ser exatamente igual ao usado nas montagens.  

Exemplo de arquivo `taxon-set.conf`:

```bash
[all]
species_1
specie_2
specie_3
```

Criar e editar config #isso pode variar a partir das estrutura inicial dos arquivos e pastas.

```bash
echo "[all]" > taxa.txt
ls uce-resultados-busca >> taxa.txt
```

```bash
sed 's/\.lastz$//' taxa.txt > taxa.conf
```

Agora é precisamos retirar a linha *probe.matches.sqlite*

```bash
[all]
Cteniza_sp
Ctenolophus_sp
Gorgyrella_namaquensis
Heligmomerus_sp
Idiops_camelus
Idiops_carajas
Idiops_clarus
Idiops_fryi
Idiops_germaini
Idiops_guri
Idiops_kanonganus
Idiops_petiti
Idiops_pirassununguensis
Idiops_pretoriae
Idiops_rastratus
Idiops_rohdei
Idiops_sp2_RF2025
Idiops_sp3_RF2025
Moggridgea_crudeni
Neocteniza_toba
probe.matches.sqlite
Segregara_transvaalensis
Titanidiops_sp
```

*excluir execendentes* 

```bash
rm -r  taxa.txt
```
Verificar taxa.conf

```bash
more taxa.conf
```

```bash
[all]
Cteniza_sp
Ctenolophus_sp
Gorgyrella_namaquensis
Heligmomerus_sp
Idiops_camelus
Idiops_carajas
Idiops_clarus
Idiops_fryi
Idiops_germaini
Idiops_guri
Idiops_kanonganus
Idiops_petiti
Idiops_pirassununguensis
Idiops_pretoriae
Idiops_rastratus
Idiops_rohdei
Idiops_sp2_RF2025
Idiops_sp3_RF2025
Moggridgea_crudeni
Neocteniza_toba

Segregara_transvaalensis
Titanidiops_sp
```


#Agora devemos criar um novo diretorio e um subdiretório

```bash
mkdir -p taxon-set/all
```

#Agora você deve executar 

```bash
 phyluce_assembly_get_match_counts \
        --locus-db ~/uce-treinamento/uce-resultados-busca/probe.matches.sqlite \
        --taxon-list-config  ~/uce-treinamento/taxa.conf \
         --taxon-group 'all' \
        --incomplete-matrix \
        --output  ~/uce-treinamento/taxon-set/all/all-taxa-incomplete.conf
```

#agora nos devemos ter algo semelhante: [phyluce_assembly_get_match_counts.log](https://github.com/TiagoBelintani/Treinamento_Processamento_UCE_UNESP_2025/blob/main/LOGS/phyluce_assembly_get_match_counts.log)

```bash
2025-08-12 17:55:21,146 - phyluce_assembly_get_match_counts - INFO - =========== Starting phyluce_assembly_get_match_counts ==========
2025-08-12 17:55:21,146 - phyluce_assembly_get_match_counts - INFO - Version: 1.7.3
2025-08-12 17:55:21,146 - phyluce_assembly_get_match_counts - INFO - Commit: None
2025-08-12 17:55:21,146 - phyluce_assembly_get_match_counts - INFO - Argument --extend_locus_db: None
2025-08-12 17:55:21,146 - phyluce_assembly_get_match_counts - INFO - Argument --incomplete_matrix: True
2025-08-12 17:55:21,146 - phyluce_assembly_get_match_counts - INFO - Argument --keep_counts: False
2025-08-12 17:55:21,147 - phyluce_assembly_get_match_counts - INFO - Argument --locus_db: /home/tiagobelintani/uce-treinamento/uce-resultados-busca/probe.matches.sqlite
2025-08-12 17:55:21,147 - phyluce_assembly_get_match_counts - INFO - Argument --log_path: None
2025-08-12 17:55:21,147 - phyluce_assembly_get_match_counts - INFO - Argument --optimize: False
2025-08-12 17:55:21,147 - phyluce_assembly_get_match_counts - INFO - Argument --output: /home/tiagobelintani/uce-treinamento/taxon-set/all/all-taxa-incomplete.conf
2025-08-12 17:55:21,147 - phyluce_assembly_get_match_counts - INFO - Argument --random: False
2025-08-12 17:55:21,147 - phyluce_assembly_get_match_counts - INFO - Argument --sample_size: 10
2025-08-12 17:55:21,147 - phyluce_assembly_get_match_counts - INFO - Argument --samples: 10
2025-08-12 17:55:21,147 - phyluce_assembly_get_match_counts - INFO - Argument --silent: False
2025-08-12 17:55:21,147 - phyluce_assembly_get_match_counts - INFO - Argument --taxon_group: all
2025-08-12 17:55:21,147 - phyluce_assembly_get_match_counts - INFO - Argument --taxon_list_config: /home/tiagobelintani/uce-treinamento/taxa.conf
2025-08-12 17:55:21,147 - phyluce_assembly_get_match_counts - INFO - Argument --verbosity: INFO
2025-08-12 17:55:21,289 - phyluce_assembly_get_match_counts - INFO - There are 22 taxa in the taxon-group '[all]' in the config file taxa.conf
2025-08-12 17:55:21,289 - phyluce_assembly_get_match_counts - INFO - Getting UCE names from database
2025-08-12 17:55:21,887 - phyluce_assembly_get_match_counts - INFO - There are 1120 total UCE loci in the database
2025-08-12 17:55:21,988 - phyluce_assembly_get_match_counts - INFO - Getting UCE matches by organism to generate a INCOMPLETE matrix
2025-08-12 17:55:21,989 - phyluce_assembly_get_match_counts - INFO - There are 171 UCE loci in an INCOMPLETE matrix
2025-08-12 17:55:21,989 - phyluce_assembly_get_match_counts - INFO - Writing the taxa and loci in the data matrix to /home/tiagobelintani/uce-treinamento/taxon-set/all/all-taxa-incomplete.conf
2025-08-12 17:56:01,103 - phyluce_assembly_get_match_counts - INFO - =========== Starting phyluce_assembly_get_match_counts ==========
2025-08-12 17:56:01,104 - phyluce_assembly_get_match_counts - INFO - Version: 1.7.3
2025-08-12 17:56:01,104 - phyluce_assembly_get_match_counts - INFO - Commit: None
2025-08-12 17:56:01,104 - phyluce_assembly_get_match_counts - INFO - Argument --extend_locus_db: None
2025-08-12 17:56:01,104 - phyluce_assembly_get_match_counts - INFO - Argument --incomplete_matrix: True
2025-08-12 17:56:01,104 - phyluce_assembly_get_match_counts - INFO - Argument --keep_counts: False
2025-08-12 17:56:01,104 - phyluce_assembly_get_match_counts - INFO - Argument --locus_db: /home/tiagobelintani/uce-treinamento/uce-resultados-busca/probe.matches.sqlite
2025-08-12 17:56:01,104 - phyluce_assembly_get_match_counts - INFO - Argument --log_path: None
2025-08-12 17:56:01,104 - phyluce_assembly_get_match_counts - INFO - Argument --optimize: False
2025-08-12 17:56:01,104 - phyluce_assembly_get_match_counts - INFO - Argument --output: /home/tiagobelintani/uce-treinamento/taxon-set/all/all-taxa-incomplete.conf
2025-08-12 17:56:01,104 - phyluce_assembly_get_match_counts - INFO - Argument --random: False
2025-08-12 17:56:01,104 - phyluce_assembly_get_match_counts - INFO - Argument --sample_size: 10
2025-08-12 17:56:01,105 - phyluce_assembly_get_match_counts - INFO - Argument --samples: 10
2025-08-12 17:56:01,105 - phyluce_assembly_get_match_counts - INFO - Argument --silent: False
2025-08-12 17:56:01,105 - phyluce_assembly_get_match_counts - INFO - Argument --taxon_group: all
2025-08-12 17:56:01,105 - phyluce_assembly_get_match_counts - INFO - Argument --taxon_list_config: /home/tiagobelintani/uce-treinamento/taxa.conf
2025-08-12 17:56:01,105 - phyluce_assembly_get_match_counts - INFO - Argument --verbosity: INFO
2025-08-12 17:56:01,107 - phyluce_assembly_get_match_counts - INFO - There are 22 taxa in the taxon-group '[all]' in the config file taxa.conf
2025-08-12 17:56:01,107 - phyluce_assembly_get_match_counts - INFO - Getting UCE names from database
2025-08-12 17:56:01,113 - phyluce_assembly_get_match_counts - INFO - There are 1120 total UCE loci in the database
2025-08-12 17:56:01,166 - phyluce_assembly_get_match_counts - INFO - Getting UCE matches by organism to generate a INCOMPLETE matrix
2025-08-12 17:56:01,166 - phyluce_assembly_get_match_counts - INFO - There are 171 UCE loci in an INCOMPLETE matrix
2025-08-12 17:56:01,166 - phyluce_assembly_get_match_counts - INFO - Writing the taxa and loci in the data matrix to /home/tiagobelintani/uce-treinamento/taxon-set/all/all-taxa-incomplete.conf
2025-08-12 17:56:01,168 - phyluce_assembly_get_match_counts - INFO - ========== Completed phyluce_assembly_get_match_counts ==========
```

#Breve descrição
<div align="justify"> 
O comando executado buscou contar e organizar os loci UCE encontrados para cada táxon listado no arquivo taxa.conf, utilizando o banco de dados probe.matches.sqlite como fonte de correspondências.

O parâmetro --incomplete_matrix foi ativado, o que significa que a matriz gerada inclui loci presentes em apenas parte dos táxons (e não apenas loci compartilhados por todos). Isso é útil quando queremos aproveitar ao máximo a informação disponível, mesmo que não haja cobertura completa para todas as espécies. Uma matrix completa, seria uma matrix recuparada considerado apenas dados completos. 

O grupo de táxons definido como 'all' contem 22 espécies. No banco de dados,1.120 loci UCE no total foram recuperados, mas apenas 171 loci compuseram a matriz incompleta resultante.

O arquivo final foi salvo em:

```bash
/home/tiagobelintani/uce-treinamento/taxon-set/all/all-taxa-incomplete.conf.
```

**Este número relativamente baixo de loci no conjunto final pode indicar:**

```bash
*Montagens de baixa qualidade ou incompletas para algumas espécies.*

*Reads iniciais com baixa cobertura ou muitos gaps.*

*Diferença filogenética grande em relação às sondas utilizadas (probes), levando a menos correspondências.*

*Parâmetros de identidade ou cobertura usados anteriormente (na etapa match_contigs_to_probes) que podem ter sido restritivos demais.*

</div>
```
#Extrair as sequências FASTA para os loci do conjunto

**Entre no diretório do taxon set e crie uma pasta para logs:**

```bash
cd taxon-sets/all
```
```bash
mkdir -p log
```
*confira o diretório*

```bash
pwd
```

Extraia as sequências FASTA:

```bash
mkdir phyluce_assembly_get_fastas_from_match_counts_job.sh
```

[phyluce_assembly_get_fastas_from_match_counts_job.sh](https://github.com/TiagoBelintani/Treinamento_Processamento_UCE_UNESP_2025/blob/main/Jobs_Conf/phyluce_assembly_get_fastas_from_match_counts_job.sh)


*O job deve ser editado*

```bash
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
```






















`




