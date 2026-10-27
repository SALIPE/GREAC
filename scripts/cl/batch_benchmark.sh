#!/bin/bash
#$ -S /bin/bash
#$ -N greac_batch_benchmark
#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs
#
# Executa o benchmark REPEATS vezes com os MESMOS parâmetros (k, janela,
# threshold), refazendo o split train/test a cada repetição, com o dataset
# copiado para o disco local do nó (/tmp2/felipe).
#
# A repetição é o laço externo: cada uma é uma partição train/test independente,
# e é daí que vem o desvio padrão. Cada execução acrescenta uma linha a
#     $PROJECTHOME/output-benchmark-<org>/benchmark_results_<org>.part-<RUN_ID>.csv
# e o Julia refaz, ao fim de cada execução,
#     $PROJECTHOME/output-benchmark-<org>/benchmark_summary_<org>.csv
# com média, desvio e as matrizes de confusão somadas de TODOS os shards.
#
# UM organismo por job: k, janela e threshold são específicos do organismo.
#
# Submissão (os parâmetros são obrigatórios):
#     qsub -v WINDOW=0.002,KMER=8,THRESHOLD=0.6,REPEATS=100 batch_benchmark.sh denv
# Dividir as repetições entre jobs simultâneos (um por nó), com RUN_ID distinto:
#     qsub -v WINDOW=0.002,KMER=8,THRESHOLD=0.6,REPEAT_START=1,REPEATS=50,RUN_ID=a batch_benchmark.sh denv
#     qsub -v WINDOW=0.002,KMER=8,THRESHOLD=0.6,REPEAT_START=51,REPEATS=50,RUN_ID=b batch_benchmark.sh denv
#   Cada job grava seu próprio shard e tem seu próprio cache do GREAC. Use
#   REPEAT_START explícito nesse caso: o auto lê o que já está gravado e dois
#   jobs simultâneos escolheriam a mesma repetição inicial.
# Continuar de onde parou:
#     qsub -v WINDOW=0.002,KMER=8,THRESHOLD=0.6,REPEAT_START=auto,REPEATS=50 batch_benchmark.sh denv
#
# REPEAT_START=1 (padrão) é um lote novo e arquiva os shards anteriores. Qualquer
# outro valor acrescenta, e o job é recusado se a repetição pedida já estiver
# gravada: repeti-la com outro split faria a linha antiga contar duas vezes.
#
# MEMBERSHIPS=1 grava também o classifications_<org>.part-<RUN_ID>.csv (uma linha
# por sequência por execução; com 100 execuções vira centenas de MB no NFS).

source /home/a61491/.bashrc

set -u

# Cronometragem com o `time` do bash, só da execução do GREAC: é esse o tempo do
# método, sem a cópia do dataset nem o split. A linha vai para o log de saída do
# cluster (-o) prefixada com TIME (filtre com `grep '^TIME' <log>`) e os mesmos
# números são gravados no shard de tempos, para a média e o desvio.
# fd 3 = saída original, fd 4 = erro original: a saída e os erros do comando
# cronometrado continuam indo para os logs de sempre, e só o que o `time`
# escreve é capturado aqui.
exec 3>&1 4>&2
cronometrar_benchmark() {
    local REP=$1
    shift
    local TIMEFORMAT='%3R %3U %3S'
    local TEMPOS STATUS REAL USER SYS

    TEMPOS=$( { time "$@" 1>&3 2>&4 ; } 2>&1 )
    STATUS=$?
    read -r REAL USER SYS <<< "$TEMPOS"

    local RESULTADO=ok
    [ "$STATUS" -eq 0 ] || RESULTADO=falha

    echo "TIME script=batch_benchmark org=$ORGANISM rep=$REP k=$KMER w=$WINDOW t=$THRESHOLD status=$RESULTADO real=$REAL user=$USER sys=$SYS"
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$(date '+%F %T')" "$ORGANISM" "$REP" "$RUN_ID" "$RESULTADO" \
        "$REAL" "$USER" "$SYS" "$KMER" "$WINDOW" "$THRESHOLD" "$METRIC" \
        "$JULIA_NUM_THREADS" "$(hostname)" >> "$TIMES_SHARD"

    return "$STATUS"
}

FEHOME=/home/a61491
PROJECTHOME=$FEHOME/GREAC/GREAC
DATASETS=$FEHOME/datasets/original
BALANCEDATASET=$FEHOME/Fasta-splitter/FastaSplitter
TEMPROOT=/tmp2/felipe/greac

# Identifica o shard e o cache deste job. Jobs simultâneos precisam de RUN_ID
# diferentes; o padrão já separa por job do SGE.
RUN_ID=${RUN_ID:-${JOB_ID:-$$}${SGE_TASK_ID:+-$SGE_TASK_ID}}

REPEATS=${REPEATS:-100}
REPEAT_START=${REPEAT_START:-1}   # 1 = lote novo | N = continua em N | auto
METRIC=${METRIC:-manhattan}
MEMBERSHIPS=${MEMBERSHIPS:-1}

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
    echo "❌ -v é opção do qsub. Execução direta: $0 <organismo>," >&2
    echo "   com os parâmetros pelo ambiente: WINDOW=0.002 KMER=8 THRESHOLD=0.6 $0 denv" >&2
    exit 1
fi

# Um organismo por job: os parâmetros do benchmark são específicos dele
if [ "$#" -ne 1 ]; then
    echo "❌ Informe exatamente um organismo (k/janela/threshold são por organismo)." >&2
    echo "   Ex: qsub -v WINDOW=0.002,KMER=8,THRESHOLD=0.6 $0 denv" >&2
    exit 1
fi
ORGANISM=$1

# Sem defaults: rodar 100 execuções com o parâmetro errado passa despercebido
for VAR in WINDOW KMER THRESHOLD; do
    if [ -z "${!VAR:-}" ]; then
        echo "❌ $VAR não definido. Ex: qsub -v WINDOW=0.002,KMER=8,THRESHOLD=0.6 $0 $ORGANISM" >&2
        exit 1
    fi
done

if [ "$REPEAT_START" != "auto" ] && ! [[ "$REPEAT_START" =~ ^[1-9][0-9]*$ ]]; then
    echo "❌ REPEAT_START deve ser um inteiro >= 1 ou 'auto', recebido: $REPEAT_START" >&2
    exit 1
fi

if ! SOURCE=$(dataset_of "$ORGANISM"); then
    echo "❌ Organismo desconhecido: $ORGANISM" >&2
    exit 1
fi
if [ ! -d "$SOURCE" ]; then
    echo "❌ $ORGANISM: $SOURCE não existe" >&2
    exit 1
fi

# Diretório local exclusivo deste job: jobs paralelos no mesmo nó não pisam no
# split uns dos outros
JOBTEMP="$TEMPROOT/bench-${JOB_ID:-$$}${SGE_TASK_ID:+-$SGE_TASK_ID}"
cleanup() {
    rm -rf "$JOBTEMP"
    echo "[$(date '+%F %T')] 🧹 $JOBTEMP removido"
}
trap cleanup EXIT INT TERM
mkdir -p "$JOBTEMP"

export JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-${NSLOTS:-1}}


OUTDIR_REL=./output-benchmark-$ORGANISM
OUTDIR=$PROJECTHOME/output-benchmark-$ORGANISM
mkdir -p "$OUTDIR"

# Tempos de execução do GREAC: um shard por job, mais o resumo com média e
# desvio padrão de TODOS os shards do lote.
TIMES_SHARD=$OUTDIR/benchmark_times_$ORGANISM.part-$RUN_ID.csv
TIMES_SUMMARY=$OUTDIR/benchmark_times_summary_$ORGANISM.csv

# Shards de tempo do lote (os arquivados não entram)
times_shards() {
    ls "$OUTDIR/benchmark_times_$ORGANISM".part-*.csv 2>/dev/null
}

# Shards do lote: o arquivo base e os .part-<run-id> dos jobs. Os arquivados
# (em arquivados/) não entram.
bench_shards() {
    ls "$OUTDIR/benchmark_results_$ORGANISM.csv" \
       "$OUTDIR/benchmark_results_$ORGANISM".part-*.csv 2>/dev/null
}

# Maior repetição já gravada em qualquer shard (0 se não houver)
LAST=0
SHARDS=$(bench_shards)
if [ -n "$SHARDS" ]; then
    # shellcheck disable=SC2086
    LAST=$(awk -F, 'FNR == 1 { c = 0; for (i = 1; i <= NF; i++) if ($i == "rep") c = i; next }
                    c && $c + 0 > m { m = $c + 0 } END { print m + 0 }' $SHARDS)
fi

case "$REPEAT_START" in
    1)
        START=1
        STAMP=$(date +%Y%m%d%H%M%S)
        ARCHIVE=$OUTDIR/arquivados
        for OLD in $(bench_shards) $(times_shards) \
                   "$OUTDIR/benchmark_summary_$ORGANISM.csv" "$TIMES_SUMMARY"; do
            [ -f "$OLD" ] || continue
            mkdir -p "$ARCHIVE"
            mv "$OLD" "$ARCHIVE/${STAMP}_$(basename "$OLD")"
            echo "[$(date '+%F %T')] 📦 $(basename "$OLD") → arquivados/ (reps 1..$LAST)"
        done
        ;;
    auto)
        START=$((LAST + 1))
        ;;
    *)
        if (( REPEAT_START <= LAST )); then
            echo "❌ REPEAT_START=$REPEAT_START, mas os shards já têm até a repetição $LAST." >&2
            echo "   Use REPEAT_START=$((LAST + 1)) ou REPEAT_START=auto." >&2
            exit 1
        fi
        START=$REPEAT_START
        ;;
esac
END=$((START + REPEATS - 1))

MEMBERSHIPS_FLAG=--no-memberships
[ "$MEMBERSHIPS" = "1" ] && MEMBERSHIPS_FLAG=""

echo "[$(date '+%F %T')] benchmark $ORGANISM | repetições $START..$END"
echo "   k=$KMER janela=$WINDOW threshold=$THRESHOLD métrica=$METRIC"
echo "   nó $(hostname) | threads $JULIA_NUM_THREADS | temp $JOBTEMP"
echo "   shard: benchmark_results_$ORGANISM.part-$RUN_ID.csv"

# Só a execução do GREAC: é isto que o `time` mede, sem a cópia nem o split.
rodar_greac() {
    local REP=$1
    ( cd "$PROJECTHOME" && julia --project src/GREAC.jl --no-cache \
        --group-name "$ORGANISM" \
        -w "$WINDOW" benchmark \
        --train-dir "$TEMP_DATA/train" \
        --test-dir "$TEMP_DATA/test" \
        -m "$METRIC" \
        -k "$KMER" \
        --threshold "$THRESHOLD" \
        -o "$OUTDIR_REL" \
        --rep "$REP" \
        --run-id "$RUN_ID" \
        $MEMBERSHIPS_FLAG \
        --classifier )
}

[ -s "$TIMES_SHARD" ] || echo "timestamp,organism,rep,run_id,status,real_s,user_s,sys_s,k,window,threshold,metric,threads,host" > "$TIMES_SHARD"

for REP in $(seq "$START" "$END"); do
    echo "########## [$(date '+%F %T')] repetição $REP / $END ##########"

    # Cópia nova a cada repetição: garante uma partição realmente nova em vez de
    # o splitter reaproveitar o train/test da repetição anterior. Fora do `time`:
    # não faz parte do método.
    TEMP_DATA="$JOBTEMP/$ORGANISM"
    rm -rf "$TEMP_DATA"
    if ! cp -r "$SOURCE" "$TEMP_DATA"; then
        echo "[$(date '+%F %T')] ❌ cópia falhou, repetição $REP ignorada" >&2
        continue
    fi

    if ! "$BALANCEDATASET/testcl.sh" "$TEMP_DATA"; then
        echo "[$(date '+%F %T')] ❌ split falhou, repetição $REP ignorada" >&2
        continue
    fi

    cronometrar_benchmark "$REP" rodar_greac "$REP" \
        || echo "[$(date '+%F %T')] ❌ repetição $REP falhou, seguindo" >&2

    rm -rf "$TEMP_DATA"
done

# Resumo dos tempos: média e desvio padrão amostral (ddof=1) sobre TODOS os
# shards do lote, contando só as execuções que terminaram bem. Média e desvio
# são calculados em uma passada (Welford), para não perder precisão somando
# quadrados de valores grandes.
resumir_tempos() {
    local SHARDS
    SHARDS=$(times_shards)
    [ -n "$SHARDS" ] || { echo "[$(date '+%F %T')] sem shards de tempo para resumir" >&2; return 1; }

    # shellcheck disable=SC2086
    awk -F, -v org="$ORGANISM" -v k="$KMER" -v w="$WINDOW" -v t="$THRESHOLD" \
            -v metric="$METRIC" -v agora="$(date '+%F %T')" '
        FNR == 1 { next }                      # cabeçalho de cada shard
        $5 != "ok" { falhas++; next }
        {
            n++
            for (c = 6; c <= 8; c++) {
                x = $c + 0
                d = x - media[c]
                media[c] += d / n
                sq[c] += d * (x - media[c])
            }
            total += $6 + 0
            if (n == 1 || $6 + 0 < minimo) minimo = $6 + 0
            if (n == 1 || $6 + 0 > maximo) maximo = $6 + 0
        }
        END {
            if (n == 0) { print "sem execuções bem sucedidas" > "/dev/stderr"; exit 1 }
            for (c = 6; c <= 8; c++) dp[c] = (n > 1) ? sprintf("%.3f", sqrt(sq[c] / (n - 1))) : ""
            print "organism,n,failures,real_mean_s,real_std_s,real_min_s,real_max_s,real_total_s," \
                  "user_mean_s,user_std_s,sys_mean_s,sys_std_s,k,window,threshold,metric,generated_at"
            printf "%s,%d,%d,%.3f,%s,%.3f,%.3f,%.3f,%.3f,%s,%.3f,%s,%s,%s,%s,%s,%s\n",
                org, n, falhas + 0, media[6], dp[6], minimo, maximo, total,
                media[7], dp[7], media[8], dp[8], k, w, t, metric, agora
        }
    ' $SHARDS > "$TIMES_SUMMARY.tmp" && mv "$TIMES_SUMMARY.tmp" "$TIMES_SUMMARY"
}

if resumir_tempos; then
    echo "[$(date '+%F %T')] ⏱  tempos do GREAC (média ± dp das execuções bem sucedidas):"
    column -s, -t "$TIMES_SUMMARY" 2>/dev/null || cat "$TIMES_SUMMARY"
else
    rm -f "$TIMES_SUMMARY.tmp"
fi

echo "[$(date '+%F %T')] 🎉 lote concluído"
echo "   $OUTDIR/benchmark_results_$ORGANISM.part-$RUN_ID.csv"
echo "   $OUTDIR/benchmark_summary_$ORGANISM.csv"
echo "   $TIMES_SHARD"
echo "   $TIMES_SUMMARY"
