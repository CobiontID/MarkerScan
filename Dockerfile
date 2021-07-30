WORKDIR /Marker_pipeline
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y \
         locales \
         libtbb2 \
         unzip \
    && rm -rf /var/lib/apt/lists/*
RUN git clone https://github.com/CobiontID/Marker-pipeline.git .
RUN snakemake -p --cores 1 --use-conda --conda-prefix /opt/conda/ --conda-create-envs-only -s create-envs.smk