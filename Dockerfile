# FROM continuumio/miniconda3

# WORKDIR /app

# # fix corrupted default channels
# RUN conda config --remove-key channels || true && \
#     conda config --add channels conda-forge && \
#     conda config --add channels defaults && \
#     conda config --show channels

# RUN apt-get update && apt-get install -y \
#     build-essential pkg-config \
#     libdbus-1-dev bowtie bedtools \
#     && rm -rf /var/lib/apt/lists/*

# COPY base_environment.yml /tmp/base_environment.yml
# COPY env2_environment.yml /tmp/env2_environment.yml
# COPY agat_environment.yml /tmp/agat_environment.yml

# RUN conda env update -n base -f /tmp/base_environment.yml && conda clean --all -y

# RUN conda config --add channels https://repo.anaconda.com/pkgs/main && \
#     conda config --add channels https://repo.anaconda.com/pkgs/r && \
#     conda config --set always_yes yes

# RUN conda env create -f /tmp/env2_environment.yml && conda clean --all -y
# RUN conda env create -f /tmp/agat_environment.yml && conda clean --all -y


# #chay db
# RUN apt-get update && apt-get install -y postgresql-client && rm -rf /var/lib/apt/lists/*


# COPY . /app

# RUN chmod +x /app/scripts/entrypoint.sh && \
#     chmod +x /app/scripts/init_db.sh

# ENTRYPOINT ["/app/scripts/entrypoint.sh"]





# stage 1
FROM mambaorg/micromamba:latest AS builder

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    libdbus-1-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN micromamba config append channels conda-forge && \
    micromamba config set always_yes yes

COPY base_environment.yml env2_environment.yml /tmp/

RUN micromamba install -y -n base -f /tmp/base_environment.yml && \
    micromamba create -y -n env2 -f /tmp/env2_environment.yml && \
    micromamba clean --all -y && \
    find /opt/conda/ -type f -name '*.a' -delete && \
    find /opt/conda/ -type f -name '__pycache__' -exec rm -rf {} +

# stage 2
FROM mambaorg/micromamba:latest

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

ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH="/opt/conda/bin:$PATH"

ENTRYPOINT ["/app/scripts/entrypoint.sh"]