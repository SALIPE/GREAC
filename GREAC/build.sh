#!/bin/bash

set -e  # Parar em caso de erro

echo "🚀 Iniciando compilação com juliac.jl..."

# Cores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Verificar se Julia está instalado
if ! command -v julia &> /dev/null; then
    echo -e "${RED}❌ Julia não está instalado!${NC}"
    exit 1
fi

# Verificar versão do Julia
JULIA_VERSION=$(julia -e "println(VERSION)")
echo -e "${BLUE}📍 Julia versão: $JULIA_VERSION${NC}"

# Instalar dependências se necessário
echo -e "${YELLOW}📦 Instalando/atualizando dependências...${NC}"
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"


# Compilar usando juliac.jl
echo -e "${YELLOW}⚙️ Compilando executável nativo...${NC}"

# Detectar sistema operacional e arquitetura
OS=$(uname -s)
ARCH=$(uname -m)

case $OS in
    Linux*)
        PLATFORM="linux"
        EXE_SUFFIX=""
        ;;
    Darwin*)
        PLATFORM="macos"
        EXE_SUFFIX=""
        ;;
    MINGW*|CYGWIN*|MSYS*)
        PLATFORM="windows"
        EXE_SUFFIX=".exe"
        ;;
    *)
        PLATFORM="unknown"
        EXE_SUFFIX=""
        ;;
esac

case $ARCH in
    x86_64|amd64)
        ARCH_NAME="x64"
        ;;
    aarch64|arm64)
        ARCH_NAME="arm64"
        ;;
    *)
        ARCH_NAME="$ARCH"
        ;;
esac

EXE_NAME="greac${EXE_SUFFIX}"
DIST_NAME="greac-${PLATFORM}-${ARCH_NAME}${EXE_SUFFIX}"

echo -e "${BLUE}🎯 Compilando para: ${PLATFORM}-${ARCH_NAME}${NC}"

# Compilar com juliac.jl
julia --project=./ build_greac.jl
