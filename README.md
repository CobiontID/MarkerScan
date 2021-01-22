# Marker_pipeline
This repository contains a Snakemake pipeline which determines the species composition of the sample by SSU presence and seperates them accordingly. 

## Required input
Please copy:
1. all scripts for the directory /scripts
2. all yaml containing information regarding external programs required to be downloaded by conda in /envs 
3. the hmmerprofile SSU_Prok_Euk_Microsporidia.hmm 

to your location of choice.

To run the pipeline a yaml file containing all external parameters is needed, an example is shown below.

```
reads: /lustre/scratch116/tol/projects/darwin/sub_projects/cobiont_detection/pipeline/hmm_pipeline/readfiles/ilBlaLact1fasta.gz
genome: /lustre/scratch116/vr/projects/vgp/build/insects/ilBlaLact1/assemblies/hicanu.20200327/ilBlaLact1.unitigs.fasta
shortname: ilBlaLact1
sci_name: Blastobasis lacticolella
workingdirectory: $WORKINGDIR
scriptdir: $SCRIPTDIR
datadir: $DATADIR
envsdir: $ENVSDIR
SSUHMMfile: $SSU/SSU_Prok_Euk_Microsporidia.hmm
```

## Script to launch the pipeline

```
#run using: bsub -o snakemake.output.%J -e snakemake.error.%J -n 10 -R"select[mem>25000] rusage[mem=25000]" -M25000 ./run_snakemake.sh $configfile

# make sure conda activate is available in script
eval "$(conda shell.bash hook)"

# activate the Conda environment
conda activate snakemake
snakemake --configfile $1 --cores 10 --use-conda --conda-prefix $condaprefix
```
