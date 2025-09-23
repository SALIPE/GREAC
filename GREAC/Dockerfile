FROM ubuntu:24.10

# Define variáveis de ambiente para evitar perguntas interativas durante a instalação
ENV DEBIAN_FRONTEND=noninteractive

# Atualiza a lista de pacotes e instala as dependências necessárias
# para baixar e instalar o Julia via juliaup.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    gnupg \
    gcc \
    g++ \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN wget -qO- https://install.julialang.org | bash -s -- -y

ENV PATH="/root/.juliaup/bin:${PATH}"

WORKDIR /app

COPY . /app

CMD ["/bin/bash", "build.sh"]
