FROM continuumio/miniconda3

WORKDIR /app

# fix corrupted default channels
RUN conda config --remove-key channels || true && \
    conda config --add channels conda-forge && \
    conda config --add channels defaults && \
    conda config --show channels

RUN apt-get update && apt-get install -y \
    build-essential pkg-config \
    libdbus-1-dev bowtie bedtools \
    && rm -rf /var/lib/apt/lists/*

COPY base_environment.yml /tmp/base_environment.yml
COPY env2_environment.yml /tmp/env2_environment.yml
COPY agat_environment.yml /tmp/agat_environment.yml

RUN conda env update -n base -f /tmp/base_environment.yml && conda clean --all -y

RUN conda config --add channels https://repo.anaconda.com/pkgs/main && \
    conda config --add channels https://repo.anaconda.com/pkgs/r && \
    conda config --set always_yes yes

RUN conda env create -f /tmp/env2_environment.yml && conda clean --all -y
RUN conda env create -f /tmp/agat_environment.yml && conda clean --all -y


#chay db
RUN apt-get update && apt-get install -y postgresql-client && rm -rf /var/lib/apt/lists/*


COPY . /app

RUN chmod +x /app/scripts/entrypoint.sh && \
    chmod +x /app/scripts/init_db.sh

ENTRYPOINT ["/app/scripts/entrypoint.sh"]