
for f in config_*.yaml
do
	echo $f 
	foo=${f#"config_"}
	foo=${foo%".yaml"}
	echo $foo
	mkdir $foo
	#./run_snakemake.sh $f
	#rm -r $foo/genera/
	bsub -o snakemake.output.%J -e snakemake.error.%J -n 10 -R"select[mem>25000] rusage[mem=25000]" -M25000 ./run_snakemake.sh $f
done
