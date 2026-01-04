#!/bin/sh

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

SALIPE_HOME=/home/a61491
TEMPDIR=/tmp2/felipe

# Datasets directory
DENGUE=$SALIPE_HOME/datasets/dengue/data
HBV=$SALIPE_HOME/datasets/HBV/data
SARS=$SALIPE_HOME/datasets/sars_cov2/data
MONKEYPOX=$SALIPE_HOME/datasets/mkpx/data
HIV=$SALIPE_HOME/datasets/hiv/data

# Reference paths
REF_HIV=$SALIPE_HOME/datasets/hiv/hiv1_refseq.fasta
REF_HBV=$SALIPE_HOME/datasets/HBV/refseq.fasta
REF_DENV=$SALIPE_HOME/datasets/dengue/refseq.fasta
REF_SARS=$SALIPE_HOME/datasets/sars_cov2/SARS-CoV2_wuhan_refseq.fasta
REF_MONKEYPOX=$SALIPE_HOME/datasets/mkpx/GCF_000857045.1_ViralProj15142_genomic.fna

# Execution parameters
GREAC=$SALIPE_HOME/GREAC/scripts/cl/benchmark.sh
GREAC_FIT=$SALIPE_HOME/GREAC/scripts/cl/fit_parameters.sh
BALANCEDATASET=$SALIPE_HOME/Fasta-splitter/FastaSplitter

cleanup() {
    echo "🧹 Cleaning TEMP..."
    if [ -d "$TEMP_DATASET" ]; then
        rm -rf "$TEMP_DATASET"
        echo "✅ TEMP removed: $TEMP_DATASET"
    fi
}

# Registrar função de cleanup para execução ao sair
trap cleanup EXIT INT TERM

if [ $# -lt 4 ]; then
    echo "❌ Erro: Argumentos insuficientes"
    echo "Uso: $0 <GROUPNAME> <WINDOW> <KMERSIZE> <THRESHOLD>"
    exit 1
fi

GROUPNAME=$1
WINDOW=$2
KMERSIZE=$3 
THRESHOLD=$4

echo "📋 Parâmetros recebidos:"
echo "   - GROUPNAME: $GROUPNAME"
echo "   - WINDOW: $WINDOW"
echo "   - KMERSIZE: $KMERSIZE"
echo "   - THRESHOLD: $THRESHOLD"

case $GROUPNAME in
    denv)
        SOURCE=$DENGUE
        REF_TOTAL=$REF_DENV
        echo "✅ Dataset DENGUE selected: $SOURCE"
        ;;
    hbv)
        SOURCE=$HBV
        REF_TOTAL=$REF_HBV
        echo "✅ Dataset HBV selected: $SOURCE"
        ;;
    bees[0-9]*)
        if [[ $GROUPNAME =~ ^bees([0-9]+)$ ]]; then
            chr="${BASH_REMATCH[1]}"
            if (( chr >= 1 && chr <= 16 )); then
                REF_BEES="REF_BEES_${chr}"
                REF_TOTAL=$(eval echo \$$REF_BEES)
                SOURCE_VAR="BEES_${chr}"
                SOURCE=$(eval echo \$$SOURCE_VAR)
                echo "✅ Dataset BEES selected: $SOURCE"
            else
                echo "❌ Erro: Número fora do intervalo permitido (1–16): $chr"
                exit 1
            fi
        else
            echo "❌ Erro: Formato inválido para bees: $GROUPNAME"
            exit 1
        fi
        ;;
    hiv)
        SOURCE=$HIV
        REF_TOTAL=$REF_HIV
        echo "✅ Dataset HIV selected: $SOURCE"
        ;;
    hiv2)
        SOURCE=$HIV2
        REF_TOTAL=$REF_HIV
        echo "✅ Dataset HIV2 selected: $SOURCE"
        ;;
    sars)
        SOURCE=$SARS
        REF_TOTAL=$REF_SARS
        echo "✅ Dataset SARS selected: $SOURCE"
        ;;
    monkeypox)
        SOURCE=$MONKEYPOX
        REF_TOTAL=$REF_MONKEYPOX
        echo "✅ Dataset MONKEYPOX selected: $SOURCE"
        ;;
    *)
        echo "❌ Error: GROUPNAME invaild: $GROUPNAME"
        exit 1
        ;;
esac

if [ ! -d "$SOURCE" ]; then
    echo "❌ Error: Directory SOURCE don't exists: $SOURCE"
    exit 1
fi

METRIC=manhattan

# $$ Job number
TEMP_DATASET="$TEMPDIR/${GROUPNAME}_$$"
TEMP_DATA="$TEMP_DATASET/data"

echo "📁 Preparando diretório temporário..."
echo "   - TEMP_DATASET: $TEMP_DATASET"


if [ -d "$TEMP_DATASET" ]; then
    echo "⚠️  Removendo diretório temporário existente..."
    rm -rf "$TEMP_DATASET"
fi

mkdir -p $TEMP_DATA

echo "📦 Copying data to TEMP..."
if [ -d "$SOURCE" ]; then
    cp -r "$SOURCE/"* "$TEMP_DATA/" 2>/dev/null || echo "⚠️  Nenhum arquivo encontrado"
else
    echo "❌ Erro: Diretório de dados não existe: $SOURCE"
    exit 1
fi


echo "✅ Data copied TEMPDIR"

for i in {1..1}; do
    echo "=========================================="
    echo "🔄 Iteração $i de 1"
    echo "=========================================="
    
    $BALANCEDATASET/testcl.sh $TEMP_DATA
    
    $GREAC_FIT $TEMP_DATA/train $TEMP_DATA/test $GROUPNAME $WINDOW $KMERSIZE $TEMP_REF
    
    if [ $? -eq 0 ]; then
        echo "✅ Iteração $i concluída com sucesso"
    else
        echo "❌ Erro na iteração $i"
    fi
    echo ""
done

echo "=========================================="
echo "🎉 Processamento concluído!"
echo "=========================================="

# A função cleanup() será chamada automaticamente ao sair
exit 0