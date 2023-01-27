FROM ubuntu:22.04
WORKDIR /MarkerScanPipeline
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y \
         locales \
         libtbb2 \
         unzip \
         git \
         snakemake \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/CobiontID/MarkerScan.git /MarkerScanPipeline
RUN snakemake -p --cores 1 --use-conda --conda-prefix /opt/conda/ --conda-frontend conda --conda-create-envs-only -s /MarkerScanPipeline/create-envs.smk
