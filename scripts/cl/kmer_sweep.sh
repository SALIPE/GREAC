#!/bin/bash
#$ -S /bin/bash
#$ -N greac_kmer_sweep
#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs
#
# Sweep 3-D do GREAC (k x window x threshold) repetido REPEATS vezes, com os
# datasets copiados para o disco local do nó (/tmp2/felipe) antes de correr.
#
# A repetição é o laço EXTERNO: cada uma recopia e refaz o split train/test,
# então as repetições são partições independentes (é isso que o desvio padrão
# mede), e todas as combinações dentro de uma repetição partilham o mesmo split.
# Cada combinação é acrescentada a
#     $PROJECTHOME/output-sweep-<org>/parameter_sweep_<org>.csv
# e a média/desvio por combinação é refeita no fim de cada repetição em
#     $PROJECTHOME/output-sweep-<org>/parameter_summary_<org>.csv
#
# Submissão (todos os organismos num job):
#     qsub kmer_sweep.sh
# Um job por organismo (array job, ordem de ALL_ORGANISMS):
#     qsub -t 1-5 kmer_sweep.sh
# Subconjunto de organismos:
#     qsub kmer_sweep.sh hbv denv
# Parâmetros pelo ambiente:
#     qsub -v REPEATS=20,K_LIST="5 6 7 8",WINDOW_MAX=0.004 kmer_sweep.sh hbv
# Dividir as repetições entre jobs simultâneos (um por nó), com RUN_ID distinto:
#     qsub -v REPEAT_START=1,REPEATS=50,RUN_ID=a kmer_sweep.sh hiv
#     qsub -v REPEAT_START=51,REPEATS=50,RUN_ID=b kmer_sweep.sh hiv
#   Cada job grava seu próprio parameter_sweep_<org>.part-<RUN_ID>.csv (append
#   concorrente no mesmo arquivo sobre NFS não é atômico) e tem seu próprio cache
#   do GREAC. O summary e o notebook juntam todos os shards.
#   Use REPEAT_START explícito nesse caso: o auto lê o que já está gravado e dois
#   jobs simultâneos escolheriam a mesma repetição inicial.
#
# Continuar de onde cada organismo parou (lê a última repetição dos shards):
#     qsub -v REPEAT_START=auto,REPEATS=10 kmer_sweep.sh hbv hiv
# Continuar a partir de uma repetição específica:
#     qsub -v REPEAT_START=42,REPEATS=10 kmer_sweep.sh hiv
# REPEAT_START=1 (padrão) é um sweep novo e arquiva o CSV anterior. Qualquer outro
# valor acrescenta ao CSV, e o job é recusado se uma repetição pedida já estiver
# gravada: repeti-la com outro split substituiria os resultados antigos.
# Uma repetição cortada no meio fica com as combinações que terminaram; o auto
# segue para a próxima em vez de misturar dois splits no mesmo número.
#
# ATENÇÃO: o custo é REPEATS x |K_LIST| x |windows| x |thresholds| por organismo.
# Com os defaults: 10 x 6 x 3 x 7 = 1260 combinações de treino+classificação.

source /home/a61491/.bashrc

set -u

FEHOME=/home/a61491
PROJECTHOME=$FEHOME/GREAC/GREAC
DATASETS=$FEHOME/datasets/original
BALANCEDATASET=$FEHOME/Fasta-splitter/FastaSplitter
TEMPROOT=/tmp2/felipe/greac

# Identifica o shard e o cache deste job. Jobs simultâneos precisam de RUN_ID
# diferentes; o padrão já separa por job do SGE.
RUN_ID=${RUN_ID:-${JOB_ID:-$$}${SGE_TASK_ID:+-$SGE_TASK_ID}}

REPEATS=${REPEATS:-100}
REPEAT_START=${REPEAT_START:-1}   # 1 = sweep novo | N = continua em N | auto

WINDOW=${WINDOW:-0.001}
WINDOW_MAX=${WINDOW_MAX:-0.002}
WINDOW_STEP=${WINDOW_STEP:-0.0005}
THRESHOLD_MIN=${THRESHOLD_MIN:-0.5}
THRESHOLD_MAX=${THRESHOLD_MAX:-0.8}
THRESHOLD_STEP=${THRESHOLD_STEP:-0.05}
K_LIST=${K_LIST:-"5 6 7 8 9 10"}


ALL_ORGANISMS=(denv hbv hiv mkpx sars)

dataset_of() {
    case "$1" in
        denv) echo "$DATASETS/denv/data_clean" ;;
        hbv)  echo "$DATASETS/HBV/data_clean" ;;
        hiv)  echo "$DATASETS/hiv/data_clean" ;;
        mkpx) echo "$DATASETS/mkpx/data_clean" ;;
        sars) echo "$DATASETS/sars_cov2/data_clean" ;;
        *) return 1 ;;
    esac
}


if [ "${1:-}" = "-v" ]; then
    echo "❌ -v é opção do qsub. Execução direta: $0 [organismo...]," >&2
    echo "   com os parâmetros pelo ambiente: REPEATS=2 $0 sars" >&2
    exit 1
fi

# Organismos: linha de comando, depois array job do SGE, depois todos
if [ "$#" -gt 0 ]; then
    ORGANISMS=("$@")
elif [ -n "${SGE_TASK_ID:-}" ] && [ "${SGE_TASK_ID}" != "undefined" ]; then
    ORGANISMS=("${ALL_ORGANISMS[$((SGE_TASK_ID - 1))]}")
else
    ORGANISMS=("${ALL_ORGANISMS[@]}")
fi

# Diretório local exclusivo deste job/tarefa: jobs paralelos no mesmo nó não
# pisam no split uns dos outros
JOBTEMP="$TEMPROOT/${JOB_ID:-$$}"
cleanup() {
    rm -rf "$JOBTEMP"
    echo "[$(date '+%F %T')] 🧹 $JOBTEMP removido"
}
trap cleanup EXIT INT TERM
mkdir -p "$JOBTEMP"

if [ "$REPEAT_START" != "auto" ] && ! [[ "$REPEAT_START" =~ ^[1-9][0-9]*$ ]]; then
    echo "❌ REPEAT_START deve ser um inteiro >= 1 ou 'auto', recebido: $REPEAT_START" >&2
    exit 1
fi

sweep_csv() { echo "$PROJECTHOME/output-sweep-$1/parameter_sweep_$1.part-$RUN_ID.csv"; }

# Shards do sweep: o arquivo base e os .part-<run-id> dos jobs. Os arquivados
# (parameter_sweep_<org>_<timestamp>.csv) não entram.
sweep_shards() {
    local dir=$PROJECTHOME/output-sweep-$1
    ls "$dir/parameter_sweep_$1.csv" "$dir/parameter_sweep_$1".part-*.csv 2>/dev/null
}

# Maior repetição já gravada em qualquer shard do organismo (0 se não houver)
last_rep() {
    local shards
    shards=$(sweep_shards "$1")
    [ -n "$shards" ] || { echo 0; return; }
    # shellcheck disable=SC2086
    awk -F, 'FNR == 1 { c = 0; for (i = 1; i <= NF; i++) if ($i == "rep") c = i; next }
             c && $c + 0 > m { m = $c + 0 } END { print m + 0 }' $shards
}

# Repetição inicial de cada organismo. Tudo é checado antes de começar: um
# conflito aborta o job inteiro em vez de pular um organismo em silêncio.
declare -A START
CONFLICT=0
for ORGANISM in "${ORGANISMS[@]}"; do
    LAST=$(last_rep "$ORGANISM")
    case "$REPEAT_START" in
        1)
            START[$ORGANISM]=1
            STAMP=$(date +%Y%m%d%H%M%S)
            ARCHIVE=$PROJECTHOME/output-sweep-$ORGANISM/arquivados
            for SHARD in $(sweep_shards "$ORGANISM"); do
                mkdir -p "$ARCHIVE"
                mv "$SHARD" "$ARCHIVE/${STAMP}_$(basename "$SHARD")"
                echo "[$(date '+%F %T')] 📦 $ORGANISM: $(basename "$SHARD") → arquivados/ (reps 1..$LAST)"
            done
            ;;
        auto)
            START[$ORGANISM]=$((LAST + 1))
            ;;
        *)
            if (( REPEAT_START <= LAST )); then
                echo "❌ $ORGANISM: REPEAT_START=$REPEAT_START, mas o CSV já tem até a repetição $LAST." >&2
                echo "   Use REPEAT_START=$((LAST + 1)) ou REPEAT_START=auto." >&2
                CONFLICT=1
            fi
            START[$ORGANISM]=$REPEAT_START
            ;;
    esac
done
(( CONFLICT )) && exit 1

echo "[$(date '+%F %T')] $REPEATS repetições | organismos: ${ORGANISMS[*]} | k: $K_LIST"
for ORGANISM in "${ORGANISMS[@]}"; do
    echo "   $ORGANISM: repetições ${START[$ORGANISM]}..$((START[$ORGANISM] + REPEATS - 1))"
done
echo "   window $WINDOW..$WINDOW_MAX (passo $WINDOW_STEP) | threshold $THRESHOLD_MIN..$THRESHOLD_MAX (passo $THRESHOLD_STEP)"
echo "   nó $(hostname) | threads ${JULIA_NUM_THREADS:-padrão} | temp $JOBTEMP"
echo "   shard: parameter_sweep_<org>.part-$RUN_ID.csv "

for ((I = 0; I < REPEATS; I++)); do
    echo "########## [$(date '+%F %T')] rodada $((I + 1)) / $REPEATS ##########"

    for ORGANISM in "${ORGANISMS[@]}"; do
        REP=$((START[$ORGANISM] + I))
        if ! SOURCE=$(dataset_of "$ORGANISM"); then
            echo "[$(date '+%F %T')] ❌ $ORGANISM: organismo desconhecido, ignorado" >&2
            continue
        fi
        if [ ! -d "$SOURCE" ]; then
            echo "[$(date '+%F %T')] ❌ $ORGANISM: $SOURCE não existe, ignorado" >&2
            continue
        fi

        # Cópia nova a cada repetição: garante uma partição realmente nova em
        # vez de o splitter reaproveitar o train/test da repetição anterior
        TEMP_DATA="$JOBTEMP/$ORGANISM"
        rm -rf "$TEMP_DATA"
        echo "[$(date '+%F %T')] 📦 $ORGANISM: copiando $SOURCE -> $TEMP_DATA"
        if ! cp -r "$SOURCE" "$TEMP_DATA"; then
            echo "[$(date '+%F %T')] ❌ $ORGANISM: cópia falhou, ignorado" >&2
            continue
        fi

        if ! "$BALANCEDATASET/testcl.sh" "$TEMP_DATA"; then
            echo "[$(date '+%F %T')] ❌ $ORGANISM: split falhou, ignorado" >&2
            continue
        fi

        echo "[$(date '+%F %T')] 🔄 $ORGANISM: sweep da repetição $REP"
        ( cd "$PROJECTHOME" && julia --project src/GREAC.jl --no-cache \
            --group-name "$ORGANISM" \
            -w "$WINDOW" fit-parameters \
            --train-dir "$TEMP_DATA/train" \
            --test-dir "$TEMP_DATA/test" \
            -k $K_LIST \
            --window-max "$WINDOW_MAX" \
            --window-step "$WINDOW_STEP" \
            --threshold-min "$THRESHOLD_MIN" \
            --threshold-max "$THRESHOLD_MAX" \
            --threshold-step "$THRESHOLD_STEP" \
            --rep "$REP" \
            --run-id "$RUN_ID" \
            --classifier ) \
            || echo "[$(date '+%F %T')] ❌ $ORGANISM: julia terminou com erro na repetição $REP" >&2

        rm -rf "$TEMP_DATA"
        echo "[$(date '+%F %T')] ✅ $ORGANISM: repetição $REP concluída"
    done
done

echo "[$(date '+%F %T')] 🎉 sweep concluído"
for ORGANISM in "${ORGANISMS[@]}"; do
    echo "   $PROJECTHOME/output-sweep-$ORGANISM/parameter_summary_$ORGANISM.csv"
done
