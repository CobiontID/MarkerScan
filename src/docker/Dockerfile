FROM ubuntu:22.04
WORKDIR /MarkerScanPipeline
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y \
         locales \
         libtbb2 \
         unzip \
         update \
         wget \
         git \
         snakemake \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# Install mamba
RUN conda install -n base -c conda-forge mamba

RUN git clone https://github.com/CobiontID/MarkerScan.git /MarkerScanPipeline
RUN snakemake -p --cores 1 --use-conda --conda-prefix /opt/conda/ --conda-create-envs-only -s /MarkerScanPipeline/create-envs.smk
