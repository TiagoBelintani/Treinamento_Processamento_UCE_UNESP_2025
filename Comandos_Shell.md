
#  Comandos Básicos de Shell para o Treinamento

A seguir, uma lista essencial de comandos de terminal para navegar, manipular arquivos e trabalhar eficientemente no ambiente Linux e no cluster do GRID-UNESP.

---

##  Navegação entre diretórios

```bash
pwd          # Mostra o caminho atual (diretório onde você está)
ls           # Lista os arquivos do diretório atual
ls -l        # Lista arquivos em formato detalhado
cd pasta/    # Entra na pasta "pasta"
cd ..        # Volta um nível
cd -         # Volta para o diretório anterior
```

---

## Manipulação de arquivos e diretórios

```bash
mkdir nova_pasta        # Cria uma nova pasta
rm arquivo.txt          # Remove um arquivo
rm -r pasta/            # Remove uma pasta e seu conteúdo
cp arquivo1 arquivo2    # Copia arquivo1 para arquivo2
mv a.txt pasta/         # Move 'a.txt' para dentro de 'pasta/'
```

---

##  Comandos úteis para bioinfo

```bash
head arquivo.txt        # Mostra as 10 primeiras linhas
tail arquivo.txt        # Mostra as 10 últimas linhas
cat arquivo.fasta       # Mostra todo o conteúdo do arquivo
less arquivo.log        # Visualiza longos arquivos com navegação
grep "Rhodnius" dados.txt  # Busca linhas com a palavra "Rhodnius"
```

---

##  Trabalhando com sessões de terminal (screen)

```bash
screen -S nome_sessao   # Cria nova sessão chamada 'nome_sessao'
screen -ls              # Lista sessões abertas
screen -r nome_sessao   # Reconecta a uma sessão existente
Ctrl+A, depois D        # Sai da sessão (sem encerrar)
```

---

## ⚠ Dicas de ouro

- Use **tab** para autocompletar nomes de arquivos e pastas.
- Use **setas ↑ ↓** para navegar por comandos anteriores.
- Use `Ctrl + C` para interromper um processo em execução.
- Combine comandos com `|` e `>` para manipular saídas.

---
Tiago Belintani 2025 
