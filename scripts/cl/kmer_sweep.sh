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
# Continuar num segundo job de onde o primeiro parou:
#     qsub -v REPEAT_START=11,REPEATS=10 kmer_sweep.sh hbv
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

REPEATS=${REPEATS:-10}
REPEAT_START=${REPEAT_START:-1}

WINDOW=${WINDOW:-0.002}
WINDOW_MAX=${WINDOW_MAX:-0.003}
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
JOBTEMP="$TEMPROOT/${JOB_ID:-$$}_${SGE_TASK_ID:-0}"
cleanup() {
    rm -rf "$JOBTEMP"
    echo "[$(date '+%F %T')] 🧹 $JOBTEMP removido"
}
trap cleanup EXIT INT TERM
mkdir -p "$JOBTEMP"
# O GREAC grava o cache (kmerset, outmasks, modelo) em $HOME/.project_cache.
# Apontar o HOME do Julia para o disco local do nó mantém esse I/O fora do NFS
# e isola jobs paralelos do mesmo organismo.
JULIA_HOME_LOCAL="$JOBTEMP/home"
mkdir -p "$JULIA_HOME_LOCAL"

REPEAT_END=$((REPEAT_START + REPEATS - 1))
echo "[$(date '+%F %T')] reps $REPEAT_START..$REPEAT_END | organismos: ${ORGANISMS[*]} | k: $K_LIST"
echo "   window $WINDOW..$WINDOW_MAX (passo $WINDOW_STEP) | threshold $THRESHOLD_MIN..$THRESHOLD_MAX (passo $THRESHOLD_STEP)"
echo "   nó $(hostname) | threads $JULIA_NUM_THREADS | temp $JOBTEMP | julia $JULIA_BIN"

# Sweep novo (REPEAT_START=1) arquiva o CSV anterior; uma continuação acrescenta
if [ "$REPEAT_START" = "1" ]; then
    for ORGANISM in "${ORGANISMS[@]}"; do
        SWEEP_CSV=$PROJECTHOME/output-sweep-$ORGANISM/parameter_sweep_$ORGANISM.csv
        if [ -f "$SWEEP_CSV" ]; then
            mv "$SWEEP_CSV" "${SWEEP_CSV%.csv}_$(date +%Y%m%d%H%M%S).csv"
            echo "[$(date '+%F %T')] 📦 $ORGANISM: sweep anterior arquivado"
        fi
    done
fi

for REP in $(seq "$REPEAT_START" "$REPEAT_END"); do
    echo "########## [$(date '+%F %T')] repetição $REP / $REPEAT_END ##########"

    for ORGANISM in "${ORGANISMS[@]}"; do
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
        ( cd "$PROJECTHOME" && HOME="$JULIA_HOME_LOCAL" "$JULIA_BIN" --project src/GREAC.jl --no-cache \
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
