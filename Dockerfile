# stage 1
FROM continuumio/miniconda3:latest AS builder

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    libdbus-1-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set always_yes yes

COPY base_environment.yml env2_environment.yml /tmp/

RUN conda env update -n base -f /tmp/base_environment.yml && \
    conda env create -n env2 -f /tmp/env2_environment.yml && \
    conda clean --all -y && \
    find /opt/conda/ -type f -name '*.a' -delete && \
    find /opt/conda/ -type f -name '__pycache__' -exec rm -rf {} +

# stage 2
FROM continuumio/miniconda3:latest

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    bowtie \
    bedtools \
    postgresql-client \
    libdbus-1-3 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/conda /opt/conda

WORKDIR /app

COPY . /app
RUN chmod +x /app/scripts/entrypoint.sh

ENV PATH="/opt/conda/bin:$PATH"

ENTRYPOINT ["/app/scripts/entrypoint.sh"]