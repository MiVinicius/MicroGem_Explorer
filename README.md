# MicroGen Variant Explorer

**Versão:** 1.0.0
**Domínio:** Bioinformática / Análise de Dados

O MicroGen Variant Explorer é uma ferramenta baseada em Django projetada para automatizar a análise de dados de variantes de Sequenciamento de Nova Geração (NGS) (arquivos VCF) para genomas bacterianos. Ele fornece métricas de Controle de Qualidade (CQ), análise química genômica (razão Ti/Tv), mapeamento de densidade mutacional e anotação funcional.

## Funcionalidades

1.  **Painel de CQ**:
    *   Calcula o total de variantes, contagens de SNP/INDEL e pontuações médias de qualidade.
    *   Identifica variantes de baixa qualidade (QUAL < 20).
    *   Gera um histograma de Distribuição de Pontuação de Qualidade.

2.  **Química Genômica**:
    *   Classifica mutações em Transições (Ti) e Transversões (Tv).
    *   Calcula a razão Ti/Tv para validar a qualidade dos dados biológicos.

3.  **Densidade Mutacional (Hotspots)**:
    *   Realiza análise de janela deslizante (padrão 1kb) para identificar hotspots mutacionais.
    *   Gera um gráfico de densidade visualizando contagens de variantes ao longo do genoma.

4.  **Anotação Funcional**:
    *   Mapeia variantes para genes usando arquivos GFF3.
    *   Relata genes específicos afetados por mutações.

## Instalação

1.  **Clone o repositório** (ou baixe o código-fonte).
2.  **Crie um ambiente virtual**:
    ```bash
    python -m venv venv
    source venv/bin/activate  # No Windows: venv\Scripts\activate
    ```
3.  **Instale as dependências**:
    ```bash
    pip install -r requirements.txt
    ```
    *Dependências incluem: Django, pandas, PyVCF3, matplotlib, seaborn.*

## Uso

A ferramenta é acessada através de uma interface web intuitiva.

1.  **Inicie o Servidor**:
    ```bash
    python manage.py runserver
    ```
2.  **Acesse o Painel**:
    Abra seu navegador e vá para `http://127.0.0.1:8000/`.

3.  **Criar Nova Análise**:
    *   Clique em "Nova Análise" (ou "New Analysis").
    *   Faça upload do seu arquivo VCF (obrigatório).
    *   (Opcional) Faça upload de um arquivo GFF para anotação funcional.
    *   Defina o tamanho da janela para análise de densidade (padrão: 1000bp).
    *   Clique em "Analisar".

4.  **Visualizar Resultados**:
    *   Após o processamento, você será redirecionado para a página de relatório detalhado.
    *   Visualize métricas de CQ, gráficos de qualidade e densidade, e a tabela de anotações.
    *   Baixe o CSV completo das anotações clicando em "Baixar CSV".
    *   Navegue de volta para a lista de análises para ver o histórico.

## Saídas

Os resultados são organizados por análise e acessíveis via interface web:

*   **Relatório Web**: Visualização interativa de todas as métricas.
*   **Gráficos**:
    *   Distribuição de Qualidade (Histograma).
    *   Densidade de Mutação (Gráfico de Linha).
*   **Arquivos para Download**:
    *   `variant_annotations.csv`: Lista completa de variantes e genes afetados.

## Estrutura do Projeto

*   `analysis/`: App Django principal contendo a lógica (`vcf_analyzer.py`, `gff_parser.py`).
*   `data/`: Contém dados de exemplo para teste (`dummy.vcf`, `dummy.gff`).
*   `media/`: Armazena arquivos enviados pelo usuário (`uploads/`) e resultados gerados (`results/`).
*   `manage.py`: Ponto de entrada do Django.
*   `microgen_explorer/`: Contém arquivos de configuração do Django.
*   `requirements.txt`: Lista de dependências do projeto.

