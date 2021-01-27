#run using: bsub -o snakemake.output.%J -e snakemake.error.%J -n 10 -R"select[mem>25000] rusage[mem=25000]" -M25000 ./run_snakemake.sh

# make sure conda activate is available in script
eval "$(conda shell.bash hook)"

# activate the Conda environment
conda activate /nfs/users/nfs_e/ev3/tools/miniconda3/envs/snakemake

snakemake --configfile $1 --cores 10 --use-conda --conda-prefix /lustre/scratch116/tol/projects/darwin/sub_projects/cobiont_detection/pipeline/miniconda/
