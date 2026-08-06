# FastANI Pipeline (Docker staphb/fastani:1.33)

Pipeline para cálculo de ANI (Average Nucleotide Identity) all-vs-all entre
genomas de E. coli usando a imagem Docker `staphb/fastani:1.33`.

## Pré-requisitos

### Docker

```bash
# Verificar se está instalado e rodando
docker info
# Se não estiver rodando:
sudo systemctl start docker
```

### Pacotes Python (no venv)

```bash
pip install pandas numpy matplotlib seaborn scipy
```

## Uso básico

```bash
git clone https://github.com/p4ul01/fastani_pipeline
cd fastani_pipeline

# Coloque seus genomas (.fna/.fasta/.fa) em ./genomes/
mkdir -p genomes
cp /caminho/para/seus/genomas/*.fna genomes/

# Roda FastANI + análise (default)
bash run_fastani.sh

# Especificar pasta de genomas diferente
bash run_fastani.sh --genomes-dir /outro/caminho/fna

# Apenas análise (se já tem o fastani_output.out)
bash run_fastani.sh --skip-run --analyze
```

O CSV de grupos biológicos é **opcional** (usado só para colorir o heatmap por
probiótico/patogênico/comensal). Um modelo está em `grupos_template.csv.example`
— copie/renomeie para `grupos_template.csv` e edite com seus próprios genomas,
ou aponte para outro arquivo com `--groups`. Sem ele, a análise roda normalmente,
apenas sem as barras de cor por grupo.

## Defaults

| Parâmetro        | Default                                  |
| ----------------- | ------------------------------------------ |
| `--genomes-dir` | `./genomes` (pasta do script)            |
| `--groups`      | `./grupos_template.csv` (opcional)       |
| `--results-dir` | `./results`                              |
| `--threads`     | `8`                                      |
| `--min-frac`    | `0.5` (fração mínima do genoma alinhada) |
| `--image`       | `staphb/fastani:1.33`                    |

## Outputs gerados

Em `results/`:

| Arquivo                    | Descrição                                                            |
| -------------------------- | ---------------------------------------------------------------------- |
| `fastani_output.out`     | Saída bruta do FastANI (tabular: genome1, genome2, ANI, count, total) |
| `fastani_matrix.phylip`  | Matriz no formato PHYLIP (gerada pelo FastANI com`--matrix`)         |
| `genome_list.txt`        | Lista de genomas processados                                           |
| `fastani_matrix.csv`     | Matriz ANI simétrica N×N (gerada pela análise Python)               |
| `fastani_heatmap.png`    | Heatmap + dendrograma (clustermap) com barras de cor por grupo         |
| `fastani_dendrogram.png` | Dendrograma retangular isolado, labels coloridos por grupo             |
| `fastani_summary.csv`    | Estatísticas por genoma (mean/min/max ANI, n_genomas_ani_ge_99)       |
| `possible_outliers.csv`  | Pares com ANI < 95% (possível contaminação ou espécie diferente)   |
| `fastani.log`            | Log completo da execução                                             |

## Interpretação

### ANI values

- **ANI ≥ 95%**: mesma espécie
- **ANI 90-95%**: espécie diferente do mesmo gênero
- **ANI < 90%**: gênero diferente

### Heatmap

- Escala de cor: 95-100% (destaca diferenças intra-E. coli)
- Dendrograma agrupa genomas por similaridade
- Barras laterais coloridas:
  - 🟢 Verde = probiótico
  - 🔴 Vermelho = patogênico
  - 🔵 Azul = comensal
  - ⚪ Cinza = unknown (não mapeado no CSV)

### Outliers

Se `possible_outliers.csv` tiver pares com ANI < 95%, investigar:

- Contaminação da cultura
- Anotação taxonômica errada
- Genoma incompleto (ver `orthologous_fraction` no `.out`)

## Estrutura

```
fastani/
├── run_fastani.sh                 # Pipeline Docker + análise
├── analyze_fastani.py             # Parser + matplotlib/seaborn
├── README.md                      # Este arquivo
├── grupos_template.csv.example    # Modelo do CSV de grupos (opcional)
├── genomes/                       # Coloque aqui seus .fna/.fasta/.fa (não versionado)
└── results/                       # Gerado em runtime
    ├── fastani_output.out
    ├── fastani_matrix.csv
    ├── fastani_heatmap.png
    ├── fastani_dendrogram.png
    ├── fastani_summary.csv
    ├── possible_outliers.csv
    ├── genome_list.txt
    └── fastani.log
```

> Dica: adicione `genomes/` e `results/` ao seu `.gitignore` para não versionar
> dados de genoma e resultados gerados.

## Referências

- **FastANI**: Jain, C., Rodriguez-R, L.M., Phillippy, A.M., Konstantinidis, K.T., Aluru, S. *High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries*. Nature Communications 9, 5114 (2018). https://doi.org/10.1038/s41467-018-07641-9
  Repositório: https://github.com/ParBLiSS/FastANI
- **Imagem Docker**: `staphb/fastani` (State Public Health Bioinformatics consortium)
  Docker Hub: https://hub.docker.com/r/staphb/fastani
  Repositório de build: https://github.com/StaPH-B/docker-builds

## Troubleshooting

### Docker permission denied

```bash
sudo usermod -aG docker $USER
# Fazer logout e login novamente
```

### Imagem não baixa

```bash
docker pull staphb/fastani:1.33
# Se ainda assim falhar, tentar versão mais antiga:
bash run_fastani.sh --image staphb/fastani:1.32
```

### Poucos genomas detectados

- Verificar extensões: `.fna`, `.fasta`, `.fa` (case-sensitive no Linux)
- Verificar se a pasta tem os arquivos:
  ```bash
  ls ./genomes/*.fna | wc -l
  ```

### FastANI muito lento

- Aumentar threads: `--threads 16`
- Reduzir `--min-frac` para `0.2` (menos rigoroso, mais rápido)
- Para 100 genomas: 100×100 = 10.000 comparações (~30 min com 8 threads)

### Análise Python falha

- Verificar dependências: `pip install pandas numpy matplotlib seaborn scipy`
- Verificar se `fastani_output.out` foi gerado: `ls results/fastani_output.out`
