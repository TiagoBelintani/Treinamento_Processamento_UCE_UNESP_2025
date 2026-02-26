
#  Acesso ao GRID-UNESP com PuTTY e WinSCP (Windows)

Este guia apresenta instruções básicas para configurar e acessar o GRID-UNESP utilizando **PuTTY** (acesso via terminal SSH) e **WinSCP** (transferência de arquivos) em sistemas Windows.

---

##  Requisitos

- Conta de usuário ativa no GRID-UNESP  
- Nome de usuário (login) e senha fornecidos após o cadastro  
- Acesso à internet

---

##  Instalação

### 🔹 PuTTY
1. Acesse: https://www.putty.org  
2. Baixe e instale a versão recomendada para Windows

### 🔹 WinSCP
1. Acesse: https://winscp.net/eng/download.php  
2. Baixe e instale a versão completa

---

##  Conectando com o PuTTY

1. Abra o PuTTY  
2. Em **Host Name (or IP address)**, digite:

```bash
login.grid.unesp.br
```

3. Em **Port**, use:

```bash
22
```

4. Em **Connection type**, selecione:

```bash
SSH
```

5. Clique em **Open**  
6. Na tela preta que aparecer, digite seu login:

```bash
[seu_usuário]
```

7. Digite sua senha (ela não será exibida) e pressione `Enter`.

---

## 🔄 Transferindo arquivos com o WinSCP

1. Abra o WinSCP  
2. Em **Protocolo de arquivo**, selecione:

```plaintext
SCP
```

3. Preencha os campos:
   - **Servidor:**
```bash
login.grid.unesp.br
```

   - **Porta:**
```bash
22
```

   - **Nome de usuário:**
```bash
[seu_usuário]
```

   - **Senha:**
```bash
[sua_senha]
```

4. Clique em **Login**

Você verá duas janelas:
- À esquerda: seus arquivos locais  
- À direita: seus arquivos no cluster  

Arraste arquivos entre os lados para fazer upload/download.

---

##  Dicas

- Use o **PuTTY** para rodar comandos, enviar jobs e navegar via terminal  
- Use o **WinSCP** para transferir scripts, resultados e arquivos de dados  
- Para segurança, evite salvar senhas nos programas  

---

📖 Para mais informações, consulte:  
 [Acessando o cluster – GRID-UNESP](https://www.ncc.unesp.br/gridunesp/docs/v2/manual/01_acessando_o_cluster.html)

---

© 2025, Tiago Belintani – Laboratório de Aracnologia de Rio Claro (LARC), UNESP
