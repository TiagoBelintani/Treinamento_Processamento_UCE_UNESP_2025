# Instalação do Miniconda e PHYLUCE

Este guia descreve o processo para instalar o [Miniconda](https://docs.conda.io/en/latest/miniconda.html) e o [PHYLUCE](https://phyluce.readthedocs.io/en/latest/) — um conjunto de ferramentas para o processamento de dados UCE (Ultra-Conserved Elements).

> Compatível com Linux e macOS. Para Windows, recomenda-se WSL ou uso de máquinas virtuais.

---

## 🧱 Etapa 1 – Instalar o Miniconda

1. Acesse a [página oficial do Miniconda](https://docs.conda.io/en/latest/miniconda.html) e baixe o instalador apropriado para seu sistema.
2. Alternativamente, use o terminal:

**Linux:**
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
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

## 🧪 Etapa 2 – Criar o ambiente Conda do PHYLUCE

> O PHYLUCE depende do Python 2.7 — mantenha esse ambiente isolado.

```bash
conda create -n phyluce python=2.7 -y
conda activate phyluce
```

---

## 🧬 Etapa 3 – Instalar o PHYLUCE

Adicione os canais Conda necessários:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Instale o pacote:

```bash
conda install phyluce -y
```

---

## ✅ Etapa 4 – Testar instalação

Verifique se o `phyluce` está acessível:

```bash
phyluce_probe_get_match_counts_from_probe_run
```

Se aparecer a tela de ajuda, a instalação foi concluída com sucesso.

---

## 📦 Etapa 5 – Exportar/Importar ambiente (opcional)

Para exportar o ambiente:
```bash
conda env export --no-builds > phyluce_env.yml
```

Para que outros usuários recriem o ambiente:
```bash
conda env create -f phyluce_env.yml
```

---

## 📝 Referências

- [Documentação do PHYLUCE](https://phyluce.readthedocs.io/)
- [Conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

---

## ☕ Sugestão

Use ambientes Conda como ambientes de laboratório: bem definidos, limpos e isolados. Nada de misturar reagentes (ou pacotes) sem saber as consequências!
