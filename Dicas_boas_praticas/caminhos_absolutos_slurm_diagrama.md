# A Importância de Usar Caminhos Absolutos em Execuções no Terminal, SLURM e Scripts

> “O caminho até o sucesso começa… pelo caminho certo.”

## 1. O que é um Caminho Absoluto?

- **Caminho absoluto**: localização completa de um arquivo ou diretório, partindo da raiz `/` do sistema de arquivos.  
  Exemplo:  
  ```bash
  /home/tiago/projetos/analise_dados/script.sh
  ```
- **Caminho relativo**: localização em relação ao diretório atual (`pwd`).  
  Exemplo (partindo de `/home/tiago`):  
  ```bash
  projetos/analise_dados/script.sh
  ```

---

## 2. Por Que Usar Caminhos Absolutos?

1. **Previsibilidade**  
   O terminal ou SLURM não “adivinham” onde você está.  
   Ao rodar um script no SLURM, o `pwd` pode não ser o mesmo de quando você submeteu o job.

2. **Evitar erros silenciosos**  
   Um caminho relativo errado pode não dar erro — pode só gerar saídas no lugar errado ou sobrescrever arquivos.

3. **Reprodutibilidade**  
   Caminhos absolutos garantem que outra pessoa (ou você no futuro) rode o comando com os mesmos resultados.

4. **Independência do ambiente**  
   Em clusters e HPCs, jobs podem rodar em nós diferentes, com diretórios de trabalho variados. Caminhos absolutos reduzem surpresas.

---

## 3. SLURM e Caminhos

No SLURM, seu job pode **iniciar em um diretório inesperado**. Por padrão:
- `sbatch` mantém o `pwd` do momento de submissão — **mas** scripts podem mudar de diretório no meio da execução.
- O nó de execução pode ter diretórios montados de forma diferente (ex.: `/scratch`, `/tmp`).

**Boa prática:**
```bash
#SBATCH --chdir=/home/tiago/projetos/analise_dados
```
ou, dentro do script:
```bash
cd /home/tiago/projetos/analise_dados || exit 1
```

E sempre referencie arquivos com caminhos completos:
```bash
python /home/tiago/projetos/analise_dados/scripts/analise.py     --input /home/tiago/dados/input.csv     --output /home/tiago/resultados/output.csv
```

---

## 4. Cuidados Extras

- **Não misture caminhos relativos e absolutos** no mesmo comando: fica confuso e é difícil depurar.
- **Use variáveis** para caminhos longos, mas sempre defina-as com valores absolutos:
  ```bash
  DATA_DIR="/home/tiago/dados"
  python script.py --input "$DATA_DIR/input.csv"
  ```
- **Ambientes temporários** (`/tmp`, `/scratch`) podem ser limpos automaticamente — use caminhos absolutos e salve resultados em diretórios persistentes.

---

## 5. Checklist Rápido

Antes de executar:
- [ ] Todos os caminhos no script estão absolutos?  
- [ ] O diretório de trabalho está definido (`pwd` esperado)?  
- [ ] Os arquivos de entrada existem (`ls /caminho/arquivo`)?  
- [ ] Os resultados vão para um local seguro e persistente?

---

## 6. Exemplo de Problema Comum

**Caminho relativo (risco):**
```bash
python scripts/analise.py --input dados/input.csv
```
Se o job rodar de `/tmp/` em vez do diretório do projeto → **arquivo não encontrado**.

**Caminho absoluto (seguro):**
```bash
python /home/tiago/projeto/scripts/analise.py     --input /home/tiago/projeto/dados/input.csv
```

---

## 7. Diagrama: Fluxo de Diretórios no SLURM

```text
[Seu Computador Local]
          |
          v
    (ssh / login)
          |
          v
[ Nó de Login do Cluster ]  --(sbatch job.slurm)-->  [ Nó(s) de Execução ]
        PWD: /home/usuario                         PWD: /var/spool/slurm/job12345
                                                      |
                                                      +--> /home  (persistente)
                                                      +--> /scratch (temporário)
                                                      +--> /tmp (volátil, limpo ao final)
```

> Moral: **não confie no `pwd` inicial**. Use caminhos absolutos para arquivos importantes.


---
***Tiago Belintani***
