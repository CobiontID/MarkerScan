FROM condaforge/miniforge3:24.9.2-0

ENV SHELL=/bin/bash

RUN mamba create -q -y -c conda-forge -c bioconda -n snakemake snakemake=6.10.0 python=3.9.7 singularity tabulate=0.8.10 \
  && conda clean --all -y

RUN  echo "source activate snakemake" > ~/.bashrc
ENV PATH=/opt/conda/envs/snakemake/bin:$PATH

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y \
         locales \
         libtbb2 \
         unzip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /Marker_pipeline

COPY . .

RUN snakemake -p --cores 1 --use-conda --conda-prefix /opt/conda/envs --conda-create-envs-only -s create-envs.smk
