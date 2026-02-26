# Instalação do Miniconda e PHYLUCE

Este guia descreve o processo para instalar o [Miniconda](https://docs.conda.io/en/latest/miniconda.html) e o [PHYLUCE](https://phyluce.readthedocs.io/en/latest/). Um conjunto de ferramentas para o processamento de dados UCE (Ultra-Conserved Elements).

> Compatível com Linux e macOS. Para Windows, recomenda-se WSL ou uso de máquinas virtuais.


##  Etapa 1 – Instalar o Miniconda

1. Acesse a [página oficial do Miniconda](https://docs.conda.io/en/latest/miniconda.html) e baixe o instalador apropriado para seu sistema.
  
2. Alternativamente, use o terminal:

**Linux:**
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh

```

**macOS (Intel):**
```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

Siga as instruções da instalação (aceitar licença, escolher diretório, permitir modificação no `.bashrc` ou `.zshrc`).

Após instalar, feche e reabra o terminal ou rode:

```bash
source ~/.bashrc  # ou ~/.zshrc, dependendo do shell
```

---

##  Etapa 2 – Criar o ambiente Conda do PHYLUCE

> O PHYLUCE

https://github.com/faircloth-lab/phyluce/releases

```bash

wget https://raw.githubusercontent.com/faircloth-lab/phyluce/v1.7.3/distrib/phyluce-1.7.3-py36-Linux-conda.yml

conda env create -n phyluce_1.7.3 --file phyluce-1.7.3-py36-Linux-conda.yml

```

---

##  Etapa 3 – Testar instalação

Verifique se o `phyluce` está acessível:

```bash
phyluce_probe_get_match_counts_from_probe_run
```
Se aparecer a tela de ajuda, a instalação foi concluída com sucesso.

-
---

## Referências

- [Documentação do PHYLUCE](https://phyluce.readthedocs.io/)
- [Conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

---
