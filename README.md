# Marker_pipeline
This repository contains a Snakemake pipeline which determines the species composition of the sample by SSU presence and seperates them accordingly. 

## Required input
Please copy:
1. snakefile
2. all scripts for the directory /scripts
3. all yaml files containing information regarding external programs required to be downloaded by conda in /envs 
4. the hmmerprofile SSU_Prok_Euk_Microsporidia.hmm 

to your location of choice.

To run the pipeline a yaml file containing all external parameters is needed, an example is shown below.

```
reads: /lustre/scratch116/tol/projects/darwin/sub_projects/cobiont_detection/pipeline/hmm_pipeline/readfiles/ilBlaLact1fasta.gz
genome: /lustre/scratch116/vr/projects/vgp/build/insects/ilBlaLact1/assemblies/hicanu.20200327/ilBlaLact1.unitigs.fasta
shortname: ilBlaLact1 (dTOL acronym)
sci_name: Blastobasis lacticolella 
workingdirectory: $WORKINGDIR
scriptdir: $SCRIPTDIR
datadir: $DATADIR
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

## Overview of the pipeline

For now very little intermediate files are removed. Once the pipeline is integrated, file removal will need to be optimized.

1. Run nhmmer with SSU_Prok_Euk_Microsporidia.hmm across the draft assembly
2. Extract coordinates of matches: "{workingdirectory}/{shortname}.ProkSSU.readsinfo"
3. Get SSU locus sequence: "{workingdirectory}/{shortname}.ProkSSU.fa" and collapse by 99% ID using cd-hit: "{workingdirectory}/{shortname}.ProkSSU.reduced.fa"
4. Download SILVA DB if new version available on https://ftp.arb-silva.de/current/ARB_files/ into {datadir}/silva
5. Classify SSU regions using SILVA. Taxonomy per sequence is found in "{workingdirectory}/{shortname}.ProkSSU.reduced.SILVA.tax"
6. Determine the species composition of sample and for which families the procedure continues

In the meantime some files are downloaded if necessary:
1. Download all refseq organellar sequences from https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/ and https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/ and store in {datadir}/organelles if files not exist or older than 30 days
2. Download all genbank organellar sequences for apicomplexans (common contaminant, but sequence information is rare) via e-utils and store in {datadir}/apicomplexa if files not exist or older than 30 days
3. Download NCBI taxonomy (both names.dmp/nodes.dmp and nucl_wgs.accession2taxid/nucl_gb.accession2taxid) if not latest version. Curl statement will download checksum file from ftp and compare download dates. Database is updated regularly, so downloading is often required.
4. Download genomes for the closest relatives of the target species available. This step is not stored in {datadir} but in {workingdir}, so downloading will be performed for every sample. This step could be optimized if pipeline is going to be run for several closely related species (e.g. moths from same family). Next, this fasta file is split and masked using duskmasker. Outputfile: {workingdirectory}/kraken.relatives.masked.ffn.

The following part of the pipeline will be done very every detected family in the 'metagenomic' composition of the sample.
1. Download all available genomes (refseq if bacterial, all if eukaryotic) for the detected families and store in {datadir}/genera if files not exist or older than 30 days. Copy over to {workingdirectory}/genera. 

As these steps all write their output to a shared directory, running several samples simultaneously could cause problems. Downloading genomes using the NCBI datasets tool could with several queries at once has also proven to be unstable.

This is followed by these steps:
1. All fasta files of the contaminant families are combined, split and masked using duskmasker. Now masking is done every time pipeline is run. This could be optimized if masking was already done when downloading. Result is in {workingdirectory}/kraken.tax.masked.ffn.
2. A kraken database is created using the masked genomes of closely related species target and combined masked fasta file of family members contaminant. {workingdirectory}/krakendb
3. Kraken2 is run. Outputfiles are {workingdirectory}/kraken.output, {workingdirectory}/kraken.report
4. All reads are mapped to the draft assembly. {workingdirectory}/AllReadsGenome.paf

The following part of the pipeline will be done very every detected family in the 'metagenomic' composition of the sample.
1. Reads are extracted per bin. {workingdirectory}/{genus}/kraken.fa
2. Busco is run on the reads {workingdirectory}/{genus}/buscoReads/busco
3. Kraken reads are mapped to draft assembly. Fully aligned contigs {workingdirectory}/{genus}/{genus}.ctgs and corresponding reads {workingdirectory}/{genus}/{genus}.reads
4. Run Busco on these contigs. {workingdirectory}/{genus}/busco/busco
5. Run Nucmer {workingdirectory}/genera/{genus}.kraken.tax.ffn on these contigs. {workingdirectory}/{genus}/{genus}_vs_contigs.overview.txt
7. Combine these results and define certain set of reads to remove. {workingdirectory}/{genus}/{genus}.final_reads.fa --> concatenated across families in **{workingdirectory}/final_reads_removal.fa** and the remaining reads are in **{workingdirectory}/final_reads_target.fa.gz**.

Moreover, also a re-assembly is done.
1. Reads of draft contigs which are not fully aligned are added to the kraken reads
2. Assembly is done using hifiasm: {workingdirectory}/{genus}/hifiasm/
3. Busco on re-assembled contigs: {workingdirectory}/{genus}/buscoAssembly
4. Nucmer {workingdirectory}/genera/{genus}.kraken.tax.ffn against re-assembled contigs
5. Map reads to re-assembled contigs: {workingdirectory}/{genus}/{genus}.putative_reads.fa --> concatenated across families in **{workingdirectory}/putative_reads_removal.fa**
6. Draw circos plot

Combine all results and generate report file **{workingdirectory}/{shortname}.report.pdf**
