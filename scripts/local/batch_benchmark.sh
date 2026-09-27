#!/bin/bash
# Executa o benchmark REPS vezes com os MESMOS parâmetros (k, janela, threshold),
# refazendo o split train/test a cada repetição, e agrega tudo em
#   GREAC/output-benchmark-<GROUPNAME>/benchmark_results_<GROUPNAME>.csv  (uma linha por execução)
#   GREAC/output-benchmark-<GROUPNAME>/benchmark_summary_<GROUPNAME>.csv  (média, desvio e matrizes somadas)
#
# O summary é refeito pelo Julia ao fim de cada execução, então resultados
# parciais continuam válidos se isto for interrompido.
#
# Repetições pelo ambiente:
#   REPS=100 $0 hbv 0.002 7 0.7                  # novo (arquiva o CSV anterior)
#   REPEAT_START=101 REPS=50 $0 hbv 0.002 7 0.7  # continua na repetição 101
#   REPEAT_START=auto REPS=50 $0 hbv 0.002 7 0.7 # continua de onde o CSV parou
# Repetir um número de repetição já gravado é recusado: o split é outro, e a
# linha antiga contaria duas vezes na média.
#
# MEMBERSHIPS=1 grava também o classifications_<grupo>.csv (uma linha por
# sequência por execução; com 100 execuções isso vira centenas de MB).
set -u

PROJECTHOME=~/Desktop/GREAC/GREAC

DATASETS=~/Desktop/datasets/original
DATASETS_ncbi=~/Desktop/datasets/new_sars_denv/denv_sars_cov_2

BALANCEDATASET=~/Desktop/Fasta-splitter/FastaSplitter

if [ $# -lt 4 ]; then
    echo "❌ Erro: Argumentos insuficientes"
    echo "Uso: $0 <GROUPNAME> <WINDOW> <KMER> <THRESHOLD> [METRIC]"
    echo "Ex:  REPS=100 $0 hbv 0.002 7 0.7"
    exit 1
fi

GROUPNAME=$1
WINDOW=$2
KMER=$3
THRESHOLD=$4
METRIC=${5:-manhattan}
REPS=${REPS:-100}
REPEAT_START=${REPEAT_START:-1}   # 1 = novo | N = continua em N | auto
MEMBERSHIPS=${MEMBERSHIPS:-1}

if [ "$REPEAT_START" != "auto" ] && ! [[ "$REPEAT_START" =~ ^[1-9][0-9]*$ ]]; then
    echo "❌ Erro: REPEAT_START deve ser um inteiro >= 1 ou 'auto', recebido: $REPEAT_START"
    exit 1
fi

case $GROUPNAME in
    denv) 
        SOURCE=$DATASETS/denv/data_clean
        SOURCE_TESTE=$DATASETS_ncbi/denv_ncbi_clean ;;

    hbv)  SOURCE=$DATASETS/HBV/data_clean ;;
    hiv)  SOURCE=$DATASETS/hiv/data_clean ;;
    sars) 
        SOURCE=$DATASETS/sars_cov2/data_clean 
        SOURCE_TESTE=$DATASETS_ncbi/sars_ncbi_clean;;

    mkpx) SOURCE=$DATASETS/mkpx/data_clean ;;
    bees[0-9]*)
        chr=${GROUPNAME#bees}
        if (( chr >= 1 && chr <= 16 )); then
            SOURCE=$DATASETS/bees/data_$chr
        else
            echo "❌ Erro: Número fora do intervalo permitido (1–16): $chr"
            exit 1
        fi
        ;;
    *)
        echo "❌ Erro: GROUPNAME desconhecido: $GROUPNAME"
        exit 1
        ;;
esac

if [ ! -d "$SOURCE" ]; then
    echo "❌ Erro: Diretório SOURCE não existe: $SOURCE"
    exit 1
fi

TRAIN=$SOURCE/train
TESTDIR=$SOURCE/test
#TESTDIR=$SOURCE_TESTE/test
OUTDIR=./output-benchmark-$GROUPNAME
RESULTS_CSV=$PROJECTHOME/output-benchmark-$GROUPNAME/benchmark_results_$GROUPNAME.csv

# Maior repetição já gravada (0 se não houver CSV)
LAST=0
if [ -s "$RESULTS_CSV" ]; then
    LAST=$(awk -F, 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "rep") c = i; next }
                    c && $c + 0 > m { m = $c + 0 } END { print m + 0 }' "$RESULTS_CSV")
fi

case "$REPEAT_START" in
    1)    START=1 ;;
    auto) START=$((LAST + 1)) ;;
    *)
        if (( REPEAT_START <= LAST )); then
            echo "❌ Erro: REPEAT_START=$REPEAT_START, mas o CSV já tem até a repetição $LAST."
            echo "   Use REPEAT_START=$((LAST + 1)) ou REPEAT_START=auto."
            exit 1
        fi
        START=$REPEAT_START
        ;;
esac
END=$((START + REPS - 1))

echo "📋 Benchmark repetido:"
echo "   - GROUPNAME: $GROUPNAME"
echo "   - SOURCE:    $SOURCE"
echo "   - PARÂMETROS: k=$KMER janela=$WINDOW threshold=$THRESHOLD métrica=$METRIC"
echo "   - REPS:      $START..$END"
echo "   - SAÍDA:     GREAC/output-benchmark-$GROUPNAME/"

# Só uma execução nova arquiva; uma continuação acrescenta ao mesmo CSV.
if [ "$START" = "1" ]; then
    for OLD in "$RESULTS_CSV" "${RESULTS_CSV%results_$GROUPNAME.csv}summary_$GROUPNAME.csv"; do
        [ -f "$OLD" ] || continue
        mv "$OLD" "${OLD%.csv}_$(date +%Y%m%d%H%M%S).csv"
        echo "📦 $(basename "$OLD") anterior (reps 1..$LAST) arquivado"
    done
fi

MEMBERSHIPS_FLAG=--no-memberships
[ "$MEMBERSHIPS" = "1" ] && MEMBERSHIPS_FLAG=""

for REP in $(seq "$START" "$END"); do
    echo "=== [$(date '+%F %T')] Repetição $REP ($START..$END)"

    if [ -x "$BALANCEDATASET/run_maxtrain.sh" ]; then
       # $BALANCEDATASET/run_maxtrain.sh $SOURCE
        $BALANCEDATASET/test.sh $SOURCE
    else
        echo "⚠️  $BALANCEDATASET/run_maxtrain.sh não encontrado, usando o split atual"
    fi

    ( cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache \
        --group-name "$GROUPNAME" \
        -w "$WINDOW" benchmark \
        --train-dir "$TRAIN" \
        --test-dir "$TESTDIR" \
        -m "$METRIC" \
        -k "$KMER" \
        --threshold "$THRESHOLD" \
        -o "$OUTDIR" \
        --rep "$REP" \
        $MEMBERSHIPS_FLAG \
        --classifier ) \
        || echo "❌ Repetição $REP falhou, seguindo"
done

echo "Resultados: GREAC/output-benchmark-$GROUPNAME/benchmark_results_$GROUPNAME.csv"
echo "Resumo:     GREAC/output-benchmark-$GROUPNAME/benchmark_summary_$GROUPNAME.csv"
