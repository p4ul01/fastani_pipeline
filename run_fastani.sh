#!/bin/bash
# =====================================================================
# run_fastani.sh - Pipeline FastANI all-vs-all via Docker
# =====================================================================
# Usa a imagem staphb/fastani:1.33 para calcular ANI (Average Nucleotide
# Identity) par-a-par entre todos os genomas da pasta informada.
#
# Saida:
#   results/fastani_output.out      - saida bruta do FastANI (tabular)
#   results/fastani_matrix.phylip   - matriz no formato phylip (se --matrix)
#   results/genome_list.txt         - lista de genomas processados
#
# Analise (heatmap + dendrograma): bash run_fastani.sh --analyze
#
# Uso:
#   bash run_fastani.sh                              # usa ./genomes por default
#   bash run_fastani.sh --genomes-dir /caminho/fna    # outra pasta de genomas
#   bash run_fastani.sh --analyze                     # so roda a analise (ja tem o .out)
#   bash run_fastani.sh --skip-run --analyze          # pula FastANI, so analisa
#
# Todos os defaults podem ser sobrescritos via flags (--genomes-dir,
# --groups, --results-dir, --threads, --min-frac, --image).
# =====================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Defaults
GENOMES_DIR="${SCRIPT_DIR}/genomes"
GROUPS_CSV="${SCRIPT_DIR}/grupos_template.csv"
RESULTS_DIR="${SCRIPT_DIR}/results"
DOCKER_IMAGE="staphb/fastani:1.33"
THREADS=8
MIN_FRAC=0.5
SKIP_RUN=0
ANALYZE_ONLY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genomes-dir)   GENOMES_DIR="$2"; shift 2 ;;
        --groups)        GROUPS_CSV="$2"; shift 2 ;;
        --results-dir)   RESULTS_DIR="$2"; shift 2 ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --min-frac)      MIN_FRAC="$2"; shift 2 ;;
        --image)         DOCKER_IMAGE="$2"; shift 2 ;;
        --skip-run)      SKIP_RUN=1; shift ;;
        --analyze)       ANALYZE_ONLY=1; shift ;;
        -h|--help)       sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Argumento desconhecido: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "${RESULTS_DIR}"
LOG_FILE="${RESULTS_DIR}/fastani.log"

say() {
    echo "[$(date '+%H:%M:%S')] $1" | tee -a "${LOG_FILE}"
}
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >> "${LOG_FILE}"
}

# ---------------------------------------------------------------------
# ETAPA 1: Descobrir genomas (auto-detecta .fna/.fasta/.fa)
# ---------------------------------------------------------------------
GENOME_LIST="${RESULTS_DIR}/genome_list.txt"

# Tenta varios patterns de extensao
mapfile -t GENOMES < <(find "${GENOMES_DIR}" -maxdepth 1 -type f \
    \( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \) 2>/dev/null | sort)

N_GENOMES=${#GENOMES[@]}

if [ "${N_GENOMES}" -lt 2 ]; then
    say "ERRO: apenas ${N_GENOMES} genoma(s) encontrado(s) em ${GENOMES_DIR}"
    say "  Extensoes suportadas: .fna, .fasta, .fa"
    exit 1
fi

# Escreve lista de genomas (caminho absoluto dentro do container)
> "${GENOME_LIST}"
for g in "${GENOMES[@]}"; do
    echo "/data/$(basename "${g}")" >> "${GENOME_LIST}"
done

say "Genomas: ${N_GENOMES} arquivos em ${GENOMES_DIR}"
log "Lista de genomas: ${GENOME_LIST}"
log "Genomas:"
for g in "${GENOMES[@]}"; do log "  - $(basename "${g}")"; done

# ---------------------------------------------------------------------
# ETAPA 2: Verificar Docker
# ---------------------------------------------------------------------
if ! command -v docker >/dev/null 2>&1; then
    say "ERRO: docker nao encontrado no PATH"
    exit 1
fi

if ! docker info >/dev/null 2>&1; then
    say "ERRO: docker daemon nao esta rodando (ou usuario sem permissao)"
    say "  Tente: sudo systemctl start docker"
    say "  Ou adicione o usuario ao grupo docker: sudo usermod -aG docker \$USER"
    exit 1
fi

# Verifica se a imagem existe; se nao, faz pull
if ! docker image inspect "${DOCKER_IMAGE}" >/dev/null 2>&1; then
    say "Baixando imagem ${DOCKER_IMAGE}..."
    if ! docker pull "${DOCKER_IMAGE}" >> "${LOG_FILE}" 2>&1; then
        say "ERRO: falha ao baixar ${DOCKER_IMAGE}"
        exit 1
    fi
fi
say "Docker: ${DOCKER_IMAGE}"

# ---------------------------------------------------------------------
# ETAPA 3: Rodar FastANI all-vs-all
# ---------------------------------------------------------------------
if [ "${SKIP_RUN}" -eq 1 ] || [ "${ANALYZE_ONLY}" -eq 1 ]; then
    say "Pulando FastANI (--skip-run ou --analyze)"
else
    OUTPUT_FILE="${RESULTS_DIR}/fastani_output.out"
    MATRIX_FILE="${RESULTS_DIR}/fastani_matrix.phylip"

    # Remove saidas antigas se existirem (FastANI nao sobrescreve)
    rm -f "${OUTPUT_FILE}" "${MATRIX_FILE}"

    say "FastANI all-vs-all (${N_GENOMES}x${N_GENOMES} = $((N_GENOMES*N_GENOMES)) comparacoes)..."

    set +e
    docker run --rm \
        --user "$(id -u):$(id -g)" \
        -v "${GENOMES_DIR}:/data:ro" \
        -v "${RESULTS_DIR}:/out" \
        "${DOCKER_IMAGE}" \
        fastANI \
        --ql /out/genome_list.txt \
        --rl /out/genome_list.txt \
        -o /out/fastani_output.out \
        --matrix \
        -t "${THREADS}" \
        --minFraction "${MIN_FRAC}" \
        >> "${LOG_FILE}" 2>&1
    DOCKER_EXIT=$?
    set -e

    if [ ${DOCKER_EXIT} -ne 0 ]; then
        say "ERRO: FastANI falhou com codigo ${DOCKER_EXIT}"
        say "  Verificando log detalhado..."
        tail -n 30 "${LOG_FILE}" | while read line; do
            say "  LOG: ${line}"
        done
        exit 1
    fi

    # Conta comparacoes bem-sucedidas
    if [ -f "${OUTPUT_FILE}" ]; then
        N_COMPARISONS=$(wc -l < "${OUTPUT_FILE}")
        say "FastANI: OK (${N_COMPARISONS} comparacoes)"
    else
        say "ERRO: arquivo de saida nao foi criado em ${OUTPUT_FILE}"
        exit 1
    fi
fi

# ---------------------------------------------------------------------
# ETAPA 4: Analise (matriz + heatmap + dendrograma)
# ---------------------------------------------------------------------
ANALYSIS_SCRIPT="${SCRIPT_DIR}/analyze_fastani.py"
if [ ! -f "${ANALYSIS_SCRIPT}" ]; then
    say "ERRO: ${ANALYSIS_SCRIPT} nao encontrado"
    exit 1
fi

# Verifica dependencias Python
MISSING=()
python3 -c "import pandas"      2>/dev/null || MISSING+=("pandas")
python3 -c "import numpy"       2>/dev/null || MISSING+=("numpy")
python3 -c "import matplotlib"  2>/dev/null || MISSING+=("matplotlib")
python3 -c "import seaborn"     2>/dev/null || MISSING+=("seaborn")
python3 -c "import scipy"       2>/dev/null || MISSING+=("scipy")

if [ ${#MISSING[@]} -gt 0 ]; then
    say "ERRO: pacotes Python faltando: ${MISSING[*]}"
    say "  Instale: pip install ${MISSING[*]}"
    exit 1
fi

OUTPUT_FILE="${RESULTS_DIR}/fastani_output.out"
if [ ! -f "${OUTPUT_FILE}" ]; then
    say "ERRO: ${OUTPUT_FILE} nao encontrado - rode o FastANI primeiro"
    exit 1
fi

if [ ! -f "${GROUPS_CSV}" ]; then
    say "Aviso: CSV de grupos nao encontrado em ${GROUPS_CSV} (opcional)"
    say "  A analise seguira sem colorir os genomas por grupo biologico."
    say "  Veja grupos_template.csv.example para o formato esperado."
fi

say "Analise: gerando matriz + heatmap..."
if python3 "${ANALYSIS_SCRIPT}" \
    --input "${OUTPUT_FILE}" \
    --groups "${GROUPS_CSV}" \
    --output-dir "${RESULTS_DIR}" \
    >> "${LOG_FILE}" 2>&1; then
    say "Analise: OK -> ${RESULTS_DIR}/"
else
    say "ERRO: analise falhou (ver ${LOG_FILE})"
    exit 1
fi

# Resumo final
say ""
say "Concluido!"
say "  Resultados: ${RESULTS_DIR}/"
say "  Log:        ${LOG_FILE}"
say ""
say "Arquivos gerados:"
for f in fastani_matrix.csv fastani_heatmap.png fastani_dendrogram.png \
         fastani_summary.csv possible_outliers.csv; do
    [ -f "${RESULTS_DIR}/${f}" ] && say "  - ${f}"
done
