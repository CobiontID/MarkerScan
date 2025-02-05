FROM condaforge/miniforge3:24.9.2-0

ENV SHELL=/bin/bash

RUN mamba create -q -y -c conda-forge -c bioconda -n snakemake snakemake=6.10.0 python=3.9.7 singularity \
  && conda clean --all -y

RUN  echo "source activate snakemake" > ~/.bashrc
ENV PATH=/opt/conda/envs/snakemake/bin:$PATH

WORKDIR /Marker_pipeline

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y \
         locales \
         libtbb2 \
         unzip \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/CobiontID/Marker-pipeline.git . && \
    sed -i 's/dependencies:/dependencies:\n  - python=3.9/' envs/hifiasm.yaml && \
    sed -i '/- defaults/d'  envs/hifiasm.yaml

RUN snakemake -p --cores 1 --use-conda --conda-prefix /opt/conda/envs --conda-create-envs-only -s  /MarkerScanPipeline/create-envs.smk
