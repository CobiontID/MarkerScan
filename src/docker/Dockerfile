FROM ubuntu:22.04
WORKDIR /Marker_pipeline
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y \
         locales \
         libtbb2 \
         unzip \
         git \
    && rm -rf /var/lib/apt/lists/*

#remove this code once it is a public repo
RUN mkdir /root/.ssh/

# Copy over private key, and set permissions
# Warning! Anyone who gets their hands on this image will be able
# to retrieve this private key file from the corresponding image layer
ADD id_rsa /root/.ssh/id_rsa

# Create known_hosts
RUN touch /root/.ssh/known_hosts
# Add bitbuckets key
RUN ssh-keyscan bitbucket.org >> /root/.ssh/known_hosts

RUN git clone git@github.com:CobiontID/Marker-pipeline.git .
RUN snakemake -p --cores 1 --use-conda --conda-prefix /opt/conda/ --conda-create-envs-only -s create-envs.smk
