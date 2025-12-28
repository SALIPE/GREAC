#!/bin/sh

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

# Função para limpar TEMPDIR ao sair
cleanup() {
    echo "🧹 Limpando diretório temporário..."
    if [ -d "$TEMP_DATASET" ]; then
        rm -rf "$TEMP_DATASET"
        echo "✅ Diretório temporário removido: $TEMP_DATASET"
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
        echo "✅ Dataset DENGUE selecionado: $SOURCE"
        ;;
    hbv)
        SOURCE=$HBV
        REF_TOTAL=$REF_HBV
        echo "✅ Dataset HBV selecionado: $SOURCE"
        ;;
    bees[0-9]*)
        if [[ $GROUPNAME =~ ^bees([0-9]+)$ ]]; then
            chr="${BASH_REMATCH[1]}"
            if (( chr >= 1 && chr <= 16 )); then
                REF_BEES="REF_BEES_${chr}"
                REF_TOTAL=$(eval echo \$$REF_BEES)
                SOURCE_VAR="BEES_${chr}"
                SOURCE=$(eval echo \$$SOURCE_VAR)
                echo "✅ Dataset BEES selecionado: $SOURCE"
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
        echo "✅ Dataset HIV selecionado: $SOURCE"
        ;;
    hiv2)
        SOURCE=$HIV2
        REF_TOTAL=$REF_HIV
        echo "✅ Dataset HIV2 selecionado: $SOURCE"
        ;;
    sars)
        SOURCE=$SARS
        REF_TOTAL=$REF_SARS
        echo "✅ Dataset SARS selecionado: $SOURCE"
        ;;
    monkeypox)
        SOURCE=$MONKEYPOX
        REF_TOTAL=$REF_MONKEYPOX
        echo "✅ Dataset MONKEYPOX selecionado: $SOURCE"
        ;;
    *)
        echo "❌ Erro: GROUPNAME inválido: $GROUPNAME"
        exit 1
        ;;
esac

# Verificar se SOURCE existe
if [ ! -d "$SOURCE" ]; then
    echo "❌ Erro: Diretório SOURCE não existe: $SOURCE"
    exit 1
fi

METRIC=manhattan

# Configurar diretório temporário para este job
TEMP_DATASET="$TEMPDIR/${GROUPNAME}_$$"
TEMP_TRAIN="$TEMP_DATASET/train"
TEMP_TEST="$TEMP_DATASET/test"

echo "📁 Preparando diretório temporário..."
echo "   - TEMP_DATASET: $TEMP_DATASET"

# Limpar qualquer resíduo anterior
if [ -d "$TEMP_DATASET" ]; then
    echo "⚠️  Removendo diretório temporário existente..."
    rm -rf "$TEMP_DATASET"
fi

# Criar estrutura de diretórios temporários
mkdir -p "$TEMP_TRAIN"
mkdir -p "$TEMP_TEST"

# Copiar dados para TEMPDIR
echo "📦 Copiando dados para diretório temporário..."
echo "   - Copiando TRAIN..."
if [ -d "$SOURCE/train" ]; then
    cp -r "$SOURCE/train/"* "$TEMP_TRAIN/" 2>/dev/null || echo "⚠️  Nenhum arquivo de treino encontrado"
else
    echo "❌ Erro: Diretório de treino não existe: $SOURCE/train"
    exit 1
fi

echo "   - Copiando TEST..."
if [ -d "$SOURCE/test" ]; then
    cp -r "$SOURCE/test/"* "$TEMP_TEST/" 2>/dev/null || echo "⚠️  Nenhum arquivo de teste encontrado"
else
    echo "❌ Erro: Diretório de teste não existe: $SOURCE/test"
    exit 1
fi

# Copiar referência para TEMPDIR (se necessário)
TEMP_REF="$TEMP_DATASET/$(basename $REF_TOTAL)"
if [ -f "$REF_TOTAL" ]; then
    echo "   - Copiando referência..."
    cp "$REF_TOTAL" "$TEMP_REF"
else
    echo "❌ Erro: Arquivo de referência não existe: $REF_TOTAL"
    exit 1
fi

echo "✅ Dados copiados com sucesso para TEMPDIR"
echo ""
echo "📊 Estatísticas dos dados temporários:"
echo "   - Train: $(find "$TEMP_TRAIN" -type f | wc -l) arquivos"
echo "   - Test: $(find "$TEMP_TEST" -type f | wc -l) arquivos"
echo ""

# Loop de execução
for i in {1..1}; do
    echo "=========================================="
    echo "🔄 Iteração $i de 1"
    echo "=========================================="
    
    # Balancear dataset (se necessário, usar dados temporários)
    echo "⚖️  Balanceando dataset..."
    $BALANCEDATASET/test.sh $TEMP_DATASET
    
    # Executar GREAC com dados temporários
    echo "🚀 Executando GREAC_FIT..."
    $GREAC_FIT "$TEMP_TRAIN" "$TEMP_TEST" "$GROUPNAME" "$WINDOW" "$KMERSIZE" "$TEMP_REF"
    
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