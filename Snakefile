"""
Pipeline to detect contaminants in Hifi reads
----------------------------------------------------
Requirements 
 - Conda (https://conda.io/docs/commands/conda-install.html)
 - SnakeMake (http://snakemake.readthedocs.io/en/stable/)
Basic usage:
  snakemake -p --use-conda --conda-prefix condadir --configfile config.yaml
"""
scriptdir = config["scriptdir"]
reads = config["reads"]
datadir = config["datadir"]
sciname_goi = config["sci_name"]
SSUHMMfile = config["SSUHMMfile"]
genome = config["genome"]
microsporidiadb = config["microsporidiadb"]
acaridb = config["acaridb"]
full=config["full"]
pwd=config["workingdirectory"]

rule all:
	input:
		expand("{pwd}/{name}.ProkSSU.reduced.fa",pwd=config["workingdirectory"], name=config["shortname"]),
		expand("{pwd}/{name}.ProkSSU.reduced.SILVA.genus.txt",pwd=config["workingdirectory"], name=config["shortname"]),		
		expand("{pwd}/kraken.tax.masked.ffn",pwd=config["workingdirectory"]),
		expand("{pwd}/kraken.report",pwd=config["workingdirectory"]),
		expand("{pwd}/final_assembly.fa",pwd=config["workingdirectory"]),
		expand("{pwd}/final_reads_removal.fa",pwd=config["workingdirectory"]),
		expand("{pwd}/putative_reads_removal.fa",pwd=config["workingdirectory"]),
		expand("{pwd}/{name}.report.pdf",pwd=config["workingdirectory"], name=config["shortname"]),

rule HMMscan_SSU:
	"""
	Run HMMscan with prokaryotic+viral HMM (RF00177+RF01959)
	"""
	output:
		dom = temporary("{workingdirectory}/{shortname}.ProkSSU.domout"), 
		log = temporary("{workingdirectory}/{shortname}.HMMscan.log")
	threads: 10
	conda: "envs/hmmer.yaml"
	shell:
		"""
		nhmmscan --cpu {threads} --noali --tblout {output.dom} -o {output.log} {SSUHMMfile} {genome}
		"""

rule FetchHMMReads:
	"""
	Fetch detected reads with prokaryotic 16S signature
	"""
	input:
		dom = "{workingdirectory}/{shortname}.ProkSSU.domout"
	output:
		readsinfo = "{workingdirectory}/{shortname}.ProkSSU.readsinfo",
		readsinfomicro = "{workingdirectory}/{shortname}.ProkSSU.microsporidia.readsinfo",
		readslist = temporary("{workingdirectory}/{shortname}.ProkSSU.readslist"),
		readslistmicro = temporary("{workingdirectory}/{shortname}.ProkSSU.microsporidia.readslist")
	shell:
		"""
		python {scriptdir}/GetReadsSSU_nhmmscan.py -i {input.dom} | grep -v 'RF02542.afa' > {output.readsinfo} || true
		python {scriptdir}/GetReadsSSU_nhmmscan.py -i {input.dom} | grep 'RF02542.afa' > {output.readsinfomicro} || true
		cut -f1 {output.readsinfo} > {output.readslist}
		cut -f1 {output.readsinfomicro} > {output.readslistmicro}
		"""

rule Fetch16SLoci:
	"""
	Get fasta sequences for detected reads with prokaryotic 16S signature and extract 16S locus
	"""
	input:
		readslist = "{workingdirectory}/{shortname}.ProkSSU.readslist",
		readsinfo = "{workingdirectory}/{shortname}.ProkSSU.readsinfo",
		readsinfomicro = "{workingdirectory}/{shortname}.ProkSSU.microsporidia.readsinfo",
		readslistmicro = "{workingdirectory}/{shortname}.ProkSSU.microsporidia.readslist"
	output:
		fasta16S = temporary("{workingdirectory}/{shortname}.ProkSSU.reads.fa"),
		fasta16SLoci = temporary("{workingdirectory}/{shortname}.ProkSSU.fa"),
		fasta16SLociReduced = "{workingdirectory}/{shortname}.ProkSSU.reduced.fa",
		fasta16Smicro = temporary("{workingdirectory}/{shortname}.ProkSSU.reads.microsporidia.fa"),
		fasta16SLocimicro = temporary("{workingdirectory}/{shortname}.ProkSSU.microsporidia.fa"),
		fasta16SLociReducedmicro = "{workingdirectory}/{shortname}.ProkSSU.microsporidia.reduced.fa"
	log: "{workingdirectory}/{shortname}.cdhit.log"
	conda:	"envs/cdhit.yaml"
	shell:
		"""
		seqtk subseq {genome} {input.readslist} > {output.fasta16S}
		python {scriptdir}/FetchSSUReads.py -i {input.readsinfo} -f {output.fasta16S} -o {output.fasta16SLoci}
		cd-hit-est -i {output.fasta16SLoci} -o {output.fasta16SLociReduced} -c 0.99 -T 1 -G 0 -aS 1 2> {log}
		if [ -s {input.readslistmicro} ]; then
			seqtk subseq {genome} {input.readslistmicro} > {output.fasta16Smicro}
			python {scriptdir}/FetchSSUReads.py -i {input.readsinfomicro} -f {output.fasta16Smicro} -o {output.fasta16SLocimicro}
			cd-hit-est -i {output.fasta16SLocimicro} -o {output.fasta16SLociReducedmicro} -c 0.99 -T 1 -G 0 -aS 1 2> {log}
		else
			touch {output.fasta16Smicro}
			touch {output.fasta16SLocimicro}
			touch {output.fasta16SLociReducedmicro}
		fi
		"""

rule DownloadSILVA:
	"""
	Download latest release SILVA DB
	"""
	output:
		donesilva = temporary("{workingdirectory}/silva_download.done.txt")
	shell:
		"""
		var=$(curl -L https://ftp.arb-silva.de/current/ARB_files/ | grep 'SSURef_opt.arb.gz.md5' | cut -f2 -d '\"')
		curl -R https://ftp.arb-silva.de/current/ARB_files/$var --output {datadir}/silva/$var
		filename=$(basename $var .md5)
		filenameshort=$(basename $filename .gz)
		if [ ! -d {datadir}/silva ]; then
  			mkdir {datadir}/silva
		fi
		if [ -f {datadir}/silva/SILVA_SSURef.arb ]; then
			if [ {datadir}/silva/$var -nt {datadir}/silva/SILVA_SSURef.arb ]; then
				curl -R https://ftp.arb-silva.de/current/ARB_files/$filename --output {datadir}/silva/$filename
				gunzip {datadir}/silva/$filename
				mv {datadir}/silva/$filenameshort {datadir}/silva/SILVA_SSURef.arb
			fi
		else
			curl -R https://ftp.arb-silva.de/current/ARB_files/$filename --output {datadir}/silva/$filename
			gunzip {datadir}/silva/$filename
			mv {datadir}/silva/$filenameshort {datadir}/silva/SILVA_SSURef.arb
		fi
		touch {output.donesilva}
		"""

rule DownloadOrganelles:
	"""
	Download gff flatfiles and fna of plastid and mitochondria from ftp release NCBI
	"""
	input:
		taxnames = expand("{datadir}/taxonomy/names.dmp",datadir=config["datadir"])
	output:
		doneorganelles = temporary("{workingdirectory}/organelles_download.done.txt")
	conda:	"envs/cdhit.yaml"
	shell:
		"""
		if [ ! -d {datadir}/organelles ]; then
  			mkdir {datadir}/organelles
		fi
		if [ -s {datadir}/organelles/organelles.lineage.txt ]; then
			before=$(date -d 'today - 180 days' +%s)
			timestamp=$(stat -c %y {datadir}/organelles/organelles.lineage.txt | cut -f1 -d ' ')
			timestampdate=$(date -d $timestamp +%s)
        	if [ $before -ge $timestampdate ]; then
                rm {datadir}/organelles/*
				mt=$(curl -L https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/ | grep -E 'genomic.gbff|genomic.fna' | cut -f2 -d '\"')
				pt=$(curl -L https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/ | grep -E 'genomic.gbff|genomic.fna' | cut -f2 -d '\"')
				for file in $mt;
				do
					curl -R https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/$file --output {datadir}/organelles/$file
				done
				for file in $pt;
				do
					curl -R https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/$file  --output {datadir}/organelles/$file
				done
				python {scriptdir}/OrganelleLineage.py -d {datadir}/organelles/ -na {input.taxnames} -o {datadir}/organelles/organelles.lineage.txt
				cat {datadir}/organelles/*genomic.fna.gz | gunzip > {datadir}/organelles/organelles.fna
				rm {datadir}/organelles/*genomic.fna.gz
				rm {datadir}/organelles/*gbff.gz
			fi
		else
			mt=$(curl -L https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/ | grep -E 'genomic.gbff|genomic.fna' | cut -f2 -d '\"')
			pt=$(curl -L https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/ | grep -E 'genomic.gbff|genomic.fna' | cut -f2 -d '\"')
			for file in $mt;
			do
				curl -R https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/$file --output {datadir}/organelles/$file
			done
			for file in $pt;
			do
				curl -R https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/$file  --output {datadir}/organelles/$file
			done
			python {scriptdir}/OrganelleLineage.py -d {datadir}/organelles/ -na {input.taxnames} -o {datadir}/organelles/organelles.lineage.txt
			cat {datadir}/organelles/*genomic.fna.gz | gunzip > {datadir}/organelles/organelles.fna
			rm {datadir}/organelles/*genomic.fna.gz
			rm {datadir}/organelles/*gbff.gz
		fi	
		touch {output.doneorganelles}
		"""

rule DownloadApicomplexa:
	"""
	Download gff flatfiles and fna of plastid and mitochondria from ftp release NCBI
	"""
	input:
		taxnames = expand("{datadir}/taxonomy/names.dmp",datadir=config["datadir"])
	output:
		done_api = temporary("{workingdirectory}/apicomplexa_download.done.txt")
	conda:	"envs/eutils.yaml"
	shell:
		"""
		if [ ! -d {datadir}/apicomplexa ]; then
  			mkdir {datadir}/apicomplexa
		fi
		if [ -s {datadir}/apicomplexa/apicomplexa.lineage.ffn ]; then
        	before=$(date -d 'today - 180 days' +%s)
        	timestamp=$(stat -c %y {datadir}/apicomplexa/apicomplexa.lineage.ffn | cut -f1 -d ' ')
        	timestampdate=$(date -d $timestamp +%s)
        	if [ $before -ge $timestampdate ]; then
                rm {datadir}/apicomplexa/*
				esearch -db nucleotide -query "apicoplast[Title] complete genome[Title] txid5794 [Organism]" | efilter -source insd | efetch -format fasta > {datadir}/apicomplexa/apicoplast.fasta
				esearch -db nucleotide -query "mitochondrion[Title] complete genome[Title] txid5794 [Organism]" | efilter -source insd | efetch -format fasta > {datadir}/apicomplexa/mito.fasta
				python {scriptdir}/ApicomplexaLineage.py -d {datadir}/apicomplexa/ -na {input.taxnames} -o {datadir}/apicomplexa/apicomplexa.lineage.ffn
			fi
		else
			esearch -db nucleotide -query "apicoplast[Title] complete genome[Title] txid5794 [Organism]" | efilter -source insd | efetch -format fasta > {datadir}/apicomplexa/apicoplast.fasta
			esearch -db nucleotide -query "mitochondrion[Title] complete genome[Title] txid5794 [Organism]" | efilter -source insd | efetch -format fasta > {datadir}/apicomplexa/mito.fasta
			python {scriptdir}/ApicomplexaLineage.py -d {datadir}/apicomplexa/ -na {input.taxnames} -o {datadir}/apicomplexa/apicomplexa.lineage.ffn
		fi	
		touch {output.done_api}
		"""

rule DownloadNCBITaxonomy:
	"""
	Download current version of NCBI taxonomy
	"""
	input:
		taxdir = expand("{datadir}/taxonomy/",datadir=config["datadir"]),
	output:
		donefile = temporary("{workingdirectory}/taxdownload.done.txt")
	shell:
		"""
		if [ ! -d {datadir}/taxonomy ]; then
			mkdir {datadir}/taxonomy
		fi
		if [ -s {datadir}/taxonomy/names.dmp ]; then
			before=$(date -d 'today - 180 days' +%s)
			timestamp=$(stat -c %y {datadir}/taxonomy/names.dmp | cut -f1 -d ' ')
			timestampdate=$(date -d $timestamp +%s)
		fi
        if [ ! -s {datadir}/taxonomy/names.dmp ] || [ $before -ge $timestampdate ]; then
			curl -R https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz.md5 --output {input.taxdir}/taxdump.tar.gz.md5
			curl -R https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz --output taxdump.tar.gz
			tar -C {input.taxdir} -xzf taxdump.tar.gz names.dmp nodes.dmp
			curl -R https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz.md5 --output {input.taxdir}/nucl_gb.accession2taxid.gz.md5
			curl -R https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz --output nucl_gb.accession2taxid.gz
			curl -R https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz --output nucl_wgs.accession2taxid.gz
			gzip -dc nucl_gb.accession2taxid.gz > {input.taxdir}/nucl_gb.accession2taxid
			gzip -dc nucl_wgs.accession2taxid.gz > {input.taxdir}/nucl_wgs.accession2taxid
			rm nucl_gb.accession2taxid.gz nucl_wgs.accession2taxid.gz taxdump.tar.gz
		fi
		touch {output.donefile}
		"""

rule ClassifySSU:
	"""
	Classify all extracted (and reduced) 16S loci using SILVA DB to determine genera present
	"""
	input:
		fasta16SLociReduced = "{workingdirectory}/{shortname}.ProkSSU.reduced.fa",
		fasta16SLociReducedmicro = "{workingdirectory}/{shortname}.ProkSSU.microsporidia.reduced.fa",
		donetaxon = "{workingdirectory}/taxdownload.done.txt",
		donesilva = "{workingdirectory}/silva_download.done.txt"
	output:
		SILVA_output_embl = temporary("{workingdirectory}/{shortname}.ProkSSU.reduced.SILVA.embl.csv"),
		SILVA_output_silva = temporary("{workingdirectory}/{shortname}.ProkSSU.reduced.SILVA.silva.csv"),
		SILVA_output_ltp = temporary("{workingdirectory}/{shortname}.ProkSSU.reduced.SILVA.ltp.csv"),
		SILVA_output = "{workingdirectory}/{shortname}.ProkSSU.reduced.SILVA.csv",
		SILVA_tax = "{workingdirectory}/{shortname}.ProkSSU.reduced.SILVA.tax",
		blastout = "{workingdirectory}/{shortname}.ProkSSU.reduced.microsporidia.blast.txt",
		blastgenus = "{workingdirectory}/{shortname}.ProkSSU.reduced.microsporidia.genus.txt",
		aclist = temporary("{workingdirectory}/{shortname}.ProkSSU.acari.list.txt"),
		fastaAcari = "{workingdirectory}/{shortname}.ProkSSU.acari.fa",
		blastacari = "{workingdirectory}/{shortname}.ProkSSU.reduced.acari.blast.txt",
		blastgenusAcari = "{workingdirectory}/{shortname}.ProkSSU.reduced.acari.genus.txt",
		SILVA16Sgenus = "{workingdirectory}/{shortname}.ProkSSU.reduced.SILVA.genus.txt"
	params:
		taxnames = expand("{datadir}/taxonomy/names.dmp",datadir=config["datadir"]),
		taxnodes = expand("{datadir}/taxonomy/nodes.dmp",datadir=config["datadir"])
	conda: "envs/sina.yaml"
	threads: 10
	shell:
		"""
		sina -i {input.fasta16SLociReduced} -o {output.SILVA_output_embl} --db {datadir}/silva/SILVA_SSURef.arb --search --search-min-sim 0.9 -p {threads} --lca-fields tax_embl_ebi_ena --outtype csv --lca-quorum 0.8 --search-max-result 20
		sina -i {input.fasta16SLociReduced} -o {output.SILVA_output_silva} --db {datadir}/silva/SILVA_SSURef.arb --search --search-min-sim 0.9 -p {threads} --lca-fields tax_slv --outtype csv --lca-quorum 0.8 --search-max-result 20
		sina -i {input.fasta16SLociReduced} -o {output.SILVA_output_ltp} --db {datadir}/silva/SILVA_SSURef.arb --search --search-min-sim 0.9 -p {threads} --lca-fields tax_ltp --outtype csv --lca-quorum 0.8 --search-max-result 20
		cat {output.SILVA_output_embl} {output.SILVA_output_silva} {output.SILVA_output_ltp} > {output.SILVA_output}
		cut -f1,8 -d',' {output.SILVA_output} | tr ',' '\t' > {output.SILVA_tax}
		cut -f2 {output.SILVA_tax} | grep -v 'lca_tax_embl_ebi_ena' | grep -v 'lca_tax_slv' | grep -v 'lca_tax_ltp' | sort | uniq > {output.SILVA16Sgenus} && [[ -s {output.SILVA16Sgenus} ]]
		if [ -s {input.fasta16SLociReducedmicro} ]; then
			blastn -db {microsporidiadb} -query {input.fasta16SLociReducedmicro} -out {output.blastout} -outfmt 6
			python {scriptdir}/ParseBlastLineage.py -b {output.blastout} -na {params.taxnames} -no {params.taxnodes} > {output.blastgenus}
			cat {output.blastgenus} >> {output.SILVA16Sgenus}
		else
			touch {output.blastout}
			touch {output.blastgenus}
		fi
		if grep -Fq 'Acari' in {output.SILVA_tax}; then
			cat {output.SILVA_tax} | grep 'Acari' | cut -f1 | sort | uniq > {output.aclist}
			python {scriptdir}/FetchSSUFasta.py -f {input.fasta16SLociReduced} -i {output.aclist} -o {output.fastaAcari}
			blastn -db {acaridb} -query {output.fastaAcari} -out {output.blastacari} -outfmt 6
			python {scriptdir}/ParseBlastLineage.py -b {output.blastacari} -na {params.taxnames} -no {params.taxnodes} > {output.blastgenusAcari}
			cat {output.blastgenusAcari} >> {output.SILVA16Sgenus}
		else
			touch {output.aclist} {output.blastacari} {output.fastaAcari} {output.blastgenusAcari}
		fi
		"""

rule MapAllReads2Assembly:
	input:
		krakenffnall = "{workingdirectory}/kraken.tax.ffn"
	output:
		paffile = temporary("{workingdirectory}/AllReadsGenome.paf"),
		mapping = temporary("{workingdirectory}/AllReadsGenome.ctgs"),
		reads = temporary("{workingdirectory}/AllReadsGenome.reads")
	threads: 10
	conda: "envs/minimap.yaml"
	shell:
		"""
		if [ -s {input.krakenffnall} ]; then
			minimap2 -x map-pb -t {threads} {genome} {reads}  > {output.paffile}
			python {scriptdir}/PafAlignment.py -p {output.paffile} -o {output.mapping} -r {output.reads}
		else
			touch {output.paffile} {output.paffile} {output.reads}
		fi
		"""

checkpoint GetGenera:
	"""
	Get genera which were detected in SILVA DB 16 screen
	"""
	input:
		SILVA16Sgenus = expand("{pwd}/{name}.ProkSSU.reduced.SILVA.genus.txt",pwd=config["workingdirectory"], name=config["shortname"]),
		donetaxon = "{workingdirectory}/taxdownload.done.txt",
	output:
		generadir = directory("{workingdirectory}/genera")
	params:
		taxnames = expand("{datadir}/taxonomy/names.dmp",datadir=config["datadir"]),
		taxnodes = expand("{datadir}/taxonomy/nodes.dmp",datadir=config["datadir"])
	conda: "envs/dataset.yaml"
	shell:
		"""
		mkdir {output.generadir}
		if [ {full} ]; then
			python {scriptdir}/DetermineGenera.py -i {input.SILVA16Sgenus} -t family -na {params.taxnames} -no {params.taxnodes} -od {output} -suf SSU.genera_taxonomy.txt -g '{sciname_goi}'
			while read p
			do
				shortname=`echo $p | cut -d, -f1`	
				echo $p > {output.generadir}/genus.$shortname.txt
			done < {output.generadir}/euk.SSU.genera_taxonomy.txt
			while read p
			do
				pattern=" |'"
				if ! [[ $p =~ $pattern ]]
				then
					echo "Bacteria/Archaea" > {output.generadir}/genus.$p.txt
					#touch  {output.generadir}/genus.$p.txt
				fi
			done < {output.generadir}/prok.SSU.genera_taxonomy.txt
			rm {output.generadir}/euk.SSU.genera_taxonomy.txt
			rm {output.generadir}/prok.SSU.genera_taxonomy.txt
		fi
		"""

rule DownloadRefSeqGenus:
	"""
	Download RefSeq genomes (per species) of selected genera from 16S screen
	"""
	input:
		generafiles = "{workingdirectory}/genera/genus.{genus}.txt",
		doneorganelles = "{workingdirectory}/organelles_download.done.txt",
		done_api = "{workingdirectory}/apicomplexa_download.done.txt"
	params:
		taxname = "{genus}"
	conda: "envs/dataset.yaml"
	output:
		krakenffnall = temporary("{workingdirectory}/genera/{genus}.kraken.tax.ffn"),
		orglist = temporary("{workingdirectory}/genera/{genus}.organelles.list"),
		orgfasta = temporary("{workingdirectory}/genera/{genus}.organelles.ffn"),
		apifile = temporary("{workingdirectory}/genera/{genus}.additional.ffn"),
		donefile = temporary("{workingdirectory}/{genus}.refseqdownload.done.txt")
	shell:
		"""
		if [ ! -d {datadir}/genera ]; then
  			mkdir {datadir}/genera
		fi
		if [ -s {datadir}/genera/{params.taxname}.kraken.tax.ffn ]; then
			before=$(date -d 'today - 30 days' +%s)
			timestamp=$(stat -c %y {datadir}/genera/{params.taxname}.kraken.tax.ffn | cut -f1 -d ' ')
			timestampdate=$(date -d $timestamp +%s)
			if [ $before -ge $timestampdate ]; then
				if [ -d {datadir}/genera/{params.taxname} ]; then
					rm -r {datadir}/genera/{params.taxname}
				fi
				mkdir {datadir}/genera/{params.taxname}
				if grep -q Eukaryota {input.generafiles}; then
					python {scriptdir}/FetchGenomesRefSeq.py --refseq no --taxname {input.generafiles} --dir {datadir}/genera/{params.taxname} > {datadir}/genera/{params.taxname}/{params.taxname}.refseq.log
				else
					python {scriptdir}/FetchGenomesRefSeq.py --refseq yes --taxname {input.generafiles} --dir {datadir}/genera/{params.taxname} > {datadir}/genera/{params.taxname}/{params.taxname}.refseq.log
				fi
				if [ -s {datadir}/genera/{params.taxname}/{params.taxname}.refseq.log ]; then
					#unzip -d {datadir}/genera/{params.taxname}/{params.taxname}.Refseq {datadir}/genera/{params.taxname}/RefSeq.{params.taxname}.zip
					python {scriptdir}/AddTaxIDKraken.py -d {datadir}/genera/{params.taxname}/{params.taxname}.Refseq -o {datadir}/genera/{params.taxname}.kraken.tax.ffn
				fi
			fi
		else
			if [ ! -d {datadir}/genera/{params.taxname} ]; then
				mkdir {datadir}/genera/{params.taxname}
			fi
			if grep -q Eukaryota {input.generafiles}; then
				python {scriptdir}/FetchGenomesRefSeq.py --refseq no --taxname {input.generafiles} --dir {datadir}/genera/{params.taxname} > {datadir}/genera/{params.taxname}/{params.taxname}.refseq.log
			else
				python {scriptdir}/FetchGenomesRefSeq.py --refseq yes --taxname {input.generafiles} --dir {datadir}/genera/{params.taxname} > {datadir}/genera/{params.taxname}/{params.taxname}.refseq.log
			fi
			if [ -s {datadir}/genera/{params.taxname}/{params.taxname}.refseq.log ]; then
				#unzip -d {datadir}/genera/{params.taxname}/{params.taxname}.Refseq {datadir}/genera/{params.taxname}/RefSeq.{params.taxname}.zip
				python {scriptdir}/AddTaxIDKraken.py -d {datadir}/genera/{params.taxname}/{params.taxname}.Refseq -o {datadir}/genera/{params.taxname}.kraken.tax.ffn
			else
				touch {datadir}/genera/{params.taxname}.kraken.tax.ffn
			fi
		fi
		touch {output.apifile}
		if grep -q Eukaryota {input.generafiles}; then
			grep {params.taxname} {datadir}/organelles/organelles.lineage.txt > {output.orglist} || true
			python {scriptdir}/FastaSelect.py -f {datadir}/organelles/organelles.fna -l {output.orglist} -o {output.orgfasta}
			if grep -q Apicomplexa {input.generafiles}; then
				cp {datadir}/apicomplexa/apicomplexa.lineage.ffn {output.apifile}
			fi
		else
			touch {output.orglist}
			touch {output.orgfasta}
		fi
		cat {datadir}/genera/{params.taxname}.kraken.tax.ffn {output.orgfasta} {output.apifile} > {output.krakenffnall}
		touch {output.donefile}
		"""

def aggregate_kraken(wildcards):
	checkpoint_output=checkpoints.GetGenera.get(**wildcards).output[0]
	return expand ("{workingdirectory}/genera/{genus}.kraken.tax.ffn", workingdirectory=config["workingdirectory"], genus=glob_wildcards(os.path.join(checkpoint_output, 'genus.{genus}.txt')).genus)

rule concatenate_kraken_input:
	input:
		aggregate_kraken
	output:
		temporary("{workingdirectory}/kraken.tax.ffn")
	shell:
		"""
		if [ -n "{input}" ]
		then
			cat {input} > {output}
		else
			touch {output}
		fi
		"""

rule DownloadGenusRel:
	"""
	Download assemblies of closely related species to species of interest
	"""
	input:
		donetaxon = "{workingdirectory}/taxdownload.done.txt",
		krakenffnall = "{workingdirectory}/kraken.tax.ffn"
	output:
		novel_pwd = directory("{workingdirectory}/relatives/"),
		refseqlog = "{workingdirectory}/relatives/relatives.refseq.log",
		refseqdir = directory("{workingdirectory}/relatives/relatives.Refseq"),
		krakenffnrel = "{workingdirectory}/relatives/relatives.kraken.tax.ffn"
	params:
		taxnames = expand("{datadir}/taxonomy/names.dmp",datadir=config["datadir"]),
		taxnodes = expand("{datadir}/taxonomy/nodes.dmp",datadir=config["datadir"])
	conda: "envs/dataset.yaml"
	shell:
		"""
		if [ ! -d {datadir}/relatives ]; then
  			mkdir {datadir}/relatives
		fi
		if [ -s {input.krakenffnall} ]; then
			python {scriptdir}/FetchGenomesRefSeqRelatives.py --taxname '{sciname_goi}' --dir {output.novel_pwd} -na {params.taxnames} -no {params.taxnodes} > {output.refseqlog}
			python {scriptdir}/AddTaxIDKraken.py -d {output.refseqdir} -o {output.krakenffnrel}
		else
			mkdir {output.refseqdir}
			touch {output.refseqlog} {output.krakenffnrel}
		fi
		"""

checkpoint SplitFasta:
	"""
	Split downloaded assemblies fasta file depending on the number of cores
	"""
	input:
		krakenffnall = "{workingdirectory}/kraken.tax.ffn"
	output:
		splitdir = temporary(directory("{workingdirectory}/split_fasta/")),
		splitdone = temporary("{workingdirectory}/split_fasta.done.txt")
	shell:
		"""
		mkdir -p {output.splitdir}
		python {scriptdir}/FastaSplit.py -f {input.krakenffnall} -s 30000 -o {output.splitdir}
		touch {output.splitdone}
		"""

rule doMasking:
	"""
	Rule to mask repetitive regions in fasta file
	"""
	input:
		fastafile = "{workingdirectory}/split_fasta/kraken.tax.{num}.fa"
	output:
		maskedfile = temporary("{workingdirectory}/split_fasta/kraken.tax.{num}.masked.fa")
	conda: "envs/kraken.yaml"
	threads: 1
	shell:
		"""
		dustmasker -in {input.fastafile} -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' > {output.maskedfile}
		"""

def aggregate_masking(wildcards):
	checkpoint_output=checkpoints.SplitFasta.get(**wildcards).output[0]
	return expand ("{workingdirectory}/split_fasta/kraken.tax.{num}.masked.fa", workingdirectory=config["workingdirectory"], num=glob_wildcards(os.path.join(checkpoint_output, 'kraken.tax.{num}.fa')).num)

rule concatenate_masking:
	input:
		aggregate_masking
	output:
		"{workingdirectory}/kraken.tax.masked.ffn"
	shell:
		"""
		if [ -n "{input}" ]
		then
			cat {input} > {output}
		else
			touch {output}
		fi
		"""

rule CreateKrakenDB:
	"""
	Create Kraken DB for all downloaded refseq genomes
	"""
	input:
		donefile = "{workingdirectory}/taxdownload.done.txt",
		krakenffnall = "{workingdirectory}/kraken.tax.masked.ffn",
		krakenffnrel = "{workingdirectory}/relatives/relatives.kraken.tax.ffn",
		splitdone = "{workingdirectory}/split_fasta.done.txt",
		splitdir = "{workingdirectory}/split_fasta/",
		krakenfasta = "{workingdirectory}/kraken.tax.ffn"
	output:
		krakendb = directory("{workingdirectory}/krakendb")
	threads: 10
	conda: "envs/kraken.yaml"
	shell:
		"""
		if [ -s {input.krakenffnall} ]
		then
			mkdir {output.krakendb}
			mkdir {output.krakendb}/taxonomy
			cp {datadir}/taxonomy/names.dmp {datadir}/taxonomy/nodes.dmp {datadir}/taxonomy/nucl_gb.accession2taxid {datadir}/taxonomy/nucl_wgs.accession2taxid  {output.krakendb}/taxonomy
			kraken2-build --threads {threads} --add-to-library {input.krakenffnall} --db {output.krakendb} --no-masking
			kraken2-build --threads {threads} --add-to-library {input.krakenffnrel} --db {output.krakendb} --no-masking
			kraken2-build --threads {threads} --build --kmer-len 50 --db {output.krakendb}
		else
			mkdir {output.krakendb}
		fi
		rm -r {input.splitdir}
		rm {input.krakenfasta}
		"""

rule RunKraken:
	"""
	Run Kraken on Hifi reads
	"""
	input:
		krakenffnall = "{workingdirectory}/kraken.tax.masked.ffn",
		krakendb = "{workingdirectory}/krakendb"
	output:
		krakenout = "{workingdirectory}/kraken.output",
		krakenreport = "{workingdirectory}/kraken.report"
	threads: 10
	conda: "envs/kraken.yaml"
	shell:
		"""
		if [ -s {input.krakenffnall} ]
		then
			if [[ {reads} == *gz ]] 
			then
				kraken2 --gzip-compressed --threads {threads} --report {output.krakenreport} --db {input.krakendb} {reads} > {output.krakenout}
			else
				kraken2 --threads {threads} --report {output.krakenreport} --db {input.krakendb} {reads} > {output.krakenout}
			fi
			rm -r {input.krakendb}/taxonomy/*
			rm -r {input.krakendb}/library/added/*
		else
			touch {output.krakenout}
			touch {output.krakenreport}
		fi
		"""

rule ExtractReadsKraken:
	"""
	For each genus extract the classified reads and get into fasta format
	"""
	input:
		krakenout = "{workingdirectory}/kraken.output",
		krakenreport = "{workingdirectory}/kraken.report",
		generafiles = "{workingdirectory}/genera/genus.{genus}.txt"
	output:
		krakenreads = "{workingdirectory}/{genus}/kraken.reads",
		krakenfa = "{workingdirectory}/{genus}/kraken.fa"
	conda: "envs/seqtk.yaml"
	shell:
		"""
		python {scriptdir}/KrakenReadsPerGenus.py -i {input.krakenout} -rep {input.krakenreport} -g {input.generafiles} -r {output.krakenreads}
		seqtk subseq {reads} {output.krakenreads} > {output.krakenfa}
		"""

rule Map2Assembly:
	input:
		krakenffnall = "{workingdirectory}/kraken.tax.masked.ffn",
		krakenfa = "{workingdirectory}/{genus}/kraken.fa"
	output:
		paffile = temporary("{workingdirectory}/{genus}/{genus}.paf"),
		mapping = "{workingdirectory}/{genus}/{genus}.ctgs",
		contiglist = temporary("{workingdirectory}/{genus}/{genus}.ctgs.list"),
		reads = temporary("{workingdirectory}/{genus}/{genus}.reads"),
		fasta = temporary("{workingdirectory}/{genus}/{genus}.ctgs.fa")
	threads: 10
	conda: "envs/minimap.yaml"
	shell:
		"""
		if [ -s {input.krakenffnall} ]
		then
			minimap2 -x map-pb -t {threads} {genome} {input.krakenfa}  > {output.paffile}
			python {scriptdir}/PafAlignment.py -p {output.paffile} -o {output.mapping} -r {output.reads}
			grep -v 'NOT COMPLETE' {output.mapping} | cut -f1 | sort | uniq > {output.contiglist} || true
			seqtk subseq {genome} {output.contiglist} > {output.fasta}
		else
			touch {output.paffile} {output.mapping} {output.contiglist} {output.reads} {output.fasta}
		fi
		"""

rule RunBusco:
	"""
	Detect number of BUSCO genes per contig
	"""
	input:
		circgenome = "{workingdirectory}/{genus}/{genus}.ctgs.fa",
		donetaxon = "{workingdirectory}/taxdownload.done.txt"
	params:
		buscodir = directory("{workingdirectory}/{genus}/busco"),
		taxnames = expand("{datadir}/taxonomy/names.dmp",datadir=config["datadir"]),
		taxnodes = expand("{datadir}/taxonomy/nodes.dmp",datadir=config["datadir"])
	output:
		buscodbs = temporary("{workingdirectory}/{genus}/info_dbs.txt"),
		buscoini = temporary("{workingdirectory}/{genus}/config_busco.ini"),
		#proteins = "{workingdirectory}/{genus}/busco/busco/prodigal_output/predicted_genes/predicted.faa",
		completed = temporary("{workingdirectory}/{genus}/busco/done.txt")
	conda: "envs/busco.yaml"
	threads:
		10
	shell:
		"""
		if [ -s {input.circgenome} ]; then
			busco --list-datasets > {output.buscodbs}
			python {scriptdir}/BuscoConfig.py -na {params.taxnames} -no {params.taxnodes} -f {input.circgenome} -d {params.buscodir} -dl {datadir}/busco_data/ -c {threads} -db {output.buscodbs} -o {output.buscoini}
			busco --config {output.buscoini} -f || true
		else
			touch {output.buscodbs}
			touch {output.buscoini}
		fi
		touch {output.completed}
		"""

rule NucmerRefSeqContigs:
	"""
	Alignment all contigs against reference genomes
	"""
	input:
		circgenome = "{workingdirectory}/{genus}/{genus}.ctgs.fa",
		buscotable = "{workingdirectory}/{genus}/busco/done.txt",
		refseqmasked = "{workingdirectory}/genera/{genus}.kraken.tax.ffn"
	output:
		completed = temporary("{workingdirectory}/{genus}/nucmer_contigs.done.txt"),
		nucmerdelta = temporary("{workingdirectory}/{genus}/{genus}_vs_contigs.delta"),
        nucmercoords = temporary("{workingdirectory}/{genus}/{genus}_vs_contigs.coords.txt"),
		nucmercontigs = "{workingdirectory}/{genus}/{genus}_vs_contigs.overview.txt"
	conda: "envs/nucmer.yaml"
	shell:
		"""
		if [ -s {input.circgenome} ]; then
			nucmer --maxmatch --delta {output.nucmerdelta} {input.circgenome} {input.refseqmasked}
			show-coords -c -l -L 100 -r -T {output.nucmerdelta} > {output.nucmercoords}
			python {scriptdir}/ParseNucmer.py -n {output.nucmercoords} -o {output.nucmercontigs}
		else
			touch {output.nucmerdelta}
			touch {output.nucmercoords}
			touch {output.nucmercontigs}
		fi
		touch {output.completed}
		"""

rule ClusterBusco:
	"""
	Detect number of genomes in assembly based on busco genes and coverage
	"""
	input:
		assemblyinfo = "{workingdirectory}/{genus}/{genus}.ctgs",
		completed = "{workingdirectory}/{genus}/busco/done.txt",
		nucmercontigs = "{workingdirectory}/{genus}/{genus}_vs_contigs.overview.txt",
		circgenome = "{workingdirectory}/{genus}/{genus}.ctgs.fa",
		krakenfa = "{workingdirectory}/{genus}/kraken.fa",
		krakenreads = "{workingdirectory}/{genus}/kraken.reads",
		reads = "{workingdirectory}/{genus}/{genus}.reads"
	output:
		summary = "{workingdirectory}/{genus}/busco/completeness_per_contig.txt",
		finalassembly = "{workingdirectory}/{genus}/{genus}.finalassembly.fa",
		contigsid = temporary("{workingdirectory}/{genus}/{genus}.ids.txt"),
		readids = temporary("{workingdirectory}/{genus}/{genus}.readsids.txt"),
		finalreads = "{workingdirectory}/{genus}/{genus}.final_reads.fa",
		nucmercontiglist = temporary("{workingdirectory}/{genus}/{genus}.nucmer.contigs.txt"),
		buscocontiglist = temporary("{workingdirectory}/{genus}/{genus}.busco.contigs.txt"),
		unmapped = temporary("{workingdirectory}/{genus}/{genus}.unmapped.reads"),
		unmappedfa = temporary("{workingdirectory}/{genus}/{genus}.unmapped.fa"),
	conda: "envs/seqtk.yaml"
	shell:
		"""
		if [ -s {input.circgenome} ]; then
			python {scriptdir}/ParseBuscoTableMapping.py -d {input.completed} -i {input.assemblyinfo} -o {output.summary} 
			grep -v 'NOT COMPLETE' {input.nucmercontigs} | cut -f1 | sort | uniq > {output.nucmercontiglist} || true
			cut -f1 {output.summary} | sort | uniq | grep -v '^#' > {output.buscocontiglist} || true
			cat  {output.buscocontiglist} {output.nucmercontiglist} | sort | uniq > {output.contigsid} || true
			if [ -s {output.contigsid}  ]; then
				seqtk subseq {input.circgenome} {output.contigsid} > {output.finalassembly}
				python {scriptdir}/SelectReads.py -r {input.reads} -o {output.readids} -c {output.contigsid}
				seqtk subseq {input.krakenfa} {output.readids} > {output.finalreads}
				comm -23 <(sort {input.krakenreads}) <(sort {output.readids}) > {output.unmapped}
				seqtk subseq {input.krakenfa} {output.unmapped} > {output.unmappedfa}
			else
				touch {output.finalassembly}
				touch {output.readids}
				touch {output.finalreads}
				touch {output.unmapped}
				cp {input.krakenfa} {output.unmappedfa}
			fi
		else
			touch {output.summary}
			touch {output.finalassembly}
			touch {output.contigsid}
			touch {output.readids}
			touch {output.finalreads}
			touch {output.nucmercontiglist}
			touch {output.buscocontiglist}
			touch {output.unmapped}
			touch {output.unmappedfa}
		fi
		"""

rule AddMappingReads:
	"""
	Add all reads mapping to contigs detected in Map2Assembly
	"""
	input:
		readsmap = "{workingdirectory}/AllReadsGenome.reads",
		mapping = "{workingdirectory}/{genus}/{genus}.ctgs",
		krakenfa = "{workingdirectory}/{genus}/kraken.reads"
	output:
		readslist = temporary("{workingdirectory}/{genus}/{genus}.allreads"),
		finalreads = temporary("{workingdirectory}/{genus}/{genus}.finalreads"),
		finalreadfasta = "{workingdirectory}/{genus}/{genus}.reads2assemble.fa"
	conda: "envs/seqtk.yaml"
	shell:
		"""
		python {scriptdir}/MappedContigs.py -m {input.mapping} -r {input.readsmap} > {output.readslist}
		cat {output.readslist} {input.krakenfa} | sort | uniq > {output.finalreads}
		seqtk subseq {reads} {output.finalreads} > {output.finalreadfasta}
 		"""

rule Hifiasm:
	"""
	Run hifiasm assembly on kraken classfied reads
	"""
	input:
		finalreadfasta = "{workingdirectory}/{genus}/{genus}.reads2assemble.fa"
	params:
		assemblyprefix = "{workingdirectory}/{genus}/hifiasm/hifiasm"
	output:
		completed = temporary("{workingdirectory}/{genus}/assembly.done.txt"),
		dirname = directory("{workingdirectory}/{genus}/hifiasm"),
		gfa = "{workingdirectory}/{genus}/hifiasm/hifiasm.p_ctg.gfa",
		fasta = "{workingdirectory}/{genus}/hifiasm/hifiasm.p_ctg.fasta",
		fai = "{workingdirectory}/{genus}/hifiasm/hifiasm.p_ctg.fasta.fai"
	threads: 10
	conda: "envs/hifiasm.yaml"
	shell:
		"""
		if [ ! -d {output.dirname} ]; then
  			mkdir {output.dirname}
		fi
		if [ -s {input.finalreadfasta} ]; then
			linecount=$(grep -c '>' < {input.finalreadfasta})
			if [ $linecount -le 100000 ]; then
				hifiasm -o {params.assemblyprefix} -t {threads} {input.finalreadfasta} -D 10 -l 1 -s 0.999 || true
				if [ -s {output.gfa} ]; then
					awk '/^S/{{print ">"$2"\\n"$3}}' {output.gfa} | fold > {output.fasta} || true
					faidx {output.fasta}
				else
					touch {output.fasta}
					touch {output.fai}
					touch {output.gfa} 
				fi
			else
				touch {output.fasta}
				touch {output.fai}
				touch {output.gfa} 
			fi 
		else
			touch {output.gfa} 
			touch {output.fasta}
			touch {output.fai} 
		fi
		touch {output.completed}
		"""

rule RunBuscoAssembly:
	"""
	Detect number of BUSCO genes per contig
	"""
	input:
		circgenome = "{workingdirectory}/{genus}/hifiasm/hifiasm.p_ctg.fasta",
		donetaxon = "{workingdirectory}/taxdownload.done.txt"
	params:
		buscodir = directory("{workingdirectory}/{genus}/buscoAssembly"),
		taxnames = expand("{datadir}/taxonomy/names.dmp",datadir=config["datadir"]),
		taxnodes = expand("{datadir}/taxonomy/nodes.dmp",datadir=config["datadir"])
	output:
		buscodbs = temporary("{workingdirectory}/{genus}/info_dbs_assembly.txt"),
		buscoini = temporary("{workingdirectory}/{genus}/config_busco_assembly.ini"),
		completed = temporary("{workingdirectory}/{genus}/buscoAssembly/done.txt"),
	conda: "envs/busco.yaml"
	threads:
		10
	shell:
		"""
		if [ -s {input.circgenome} ]; then
			busco --list-datasets > {output.buscodbs}
			python {scriptdir}/BuscoConfig.py -na {params.taxnames} -no {params.taxnodes} -f {input.circgenome} -d {params.buscodir} -dl {datadir}/busco_data/ -c {threads} -db {output.buscodbs} -o {output.buscoini}
			busco --config {output.buscoini} -f || true
		else
			touch {output.buscodbs}
			touch {output.buscoini}
		fi
		touch {output.completed}
		"""

rule NucmerRefSeqHifiasm:
	"""
	Alignment all contigs against reference genomes
	"""
	input:
		circgenome = "{workingdirectory}/{genus}/hifiasm/hifiasm.p_ctg.fasta",
		buscotable = "{workingdirectory}/{genus}/buscoAssembly/done.txt",
		refseqmasked = "{workingdirectory}/genera/{genus}.kraken.tax.ffn"
	output:
		completed = temporary("{workingdirectory}/{genus}/nucmer_hifiasm.done.txt"),
		nucmerdelta = temporary("{workingdirectory}/{genus}/{genus}_vs_hifiasm.delta"),
        nucmercoords = temporary("{workingdirectory}/{genus}/{genus}_vs_hifiasm.coords.txt"),
		nucmercontigs = "{workingdirectory}/{genus}/{genus}_vs_hifiasm.overview.txt"
	conda: "envs/nucmer.yaml"
	shell:
		"""
		if [ -s {input.circgenome} ]; then
			nucmer --maxmatch --delta {output.nucmerdelta} {input.circgenome} {input.refseqmasked}
			show-coords -c -l -L 100 -r -T {output.nucmerdelta} > {output.nucmercoords}
			python {scriptdir}/ParseNucmer.py -n {output.nucmercoords} -o {output.nucmercontigs}
		else
			touch {output.nucmerdelta}
			touch {output.nucmercoords}
			touch {output.nucmercontigs}
		fi
		touch {output.completed}
		"""

rule Map2AssemblyHifiasm:
	input:
		krakenfa = "{workingdirectory}/{genus}/{genus}.reads2assemble.fa",
		assemblyfasta = "{workingdirectory}/{genus}/hifiasm/hifiasm.p_ctg.fasta",
		completed = "{workingdirectory}/{genus}/buscoAssembly/done.txt",
		unmapped = "{workingdirectory}/{genus}/{genus}.unmapped.reads",
		nucmercontigs = "{workingdirectory}/{genus}/{genus}_vs_hifiasm.overview.txt",
		readfile = "{workingdirectory}/{genus}/buscoReads.txt"
	output:
		summary = "{workingdirectory}/{genus}/buscoAssembly/completeness_per_contig.txt",
		buscocontiglist = temporary("{workingdirectory}/{genus}/{genus}.buscoAssembly.contigs.txt"),
		nucmercontiglist = temporary("{workingdirectory}/{genus}/{genus}.NucmerAssembly.contigs.txt"),
		contiglist = temporary("{workingdirectory}/{genus}/{genus}.Assembly.contigs.txt"),
		paffile = temporary("{workingdirectory}/{genus}/{genus}.assembly.paf"),
		fasta = "{workingdirectory}/{genus}/{genus}.putative_assembly.fa",
		mapping = temporary("{workingdirectory}/{genus}/{genus}.assembly.ctgs"),
		reads = temporary("{workingdirectory}/{genus}/{genus}.assembly.reads"),
		reads_unmapped = temporary("{workingdirectory}/{genus}/{genus}.assembly.unmapped.reads"),
		readsfasta = "{workingdirectory}/{genus}/{genus}.putative_reads.fa",
		busco_assembly = temporary("{workingdirectory}/{genus}/buscoReadsAssembly.txt"),
		busco_assembly_hifi = temporary("{workingdirectory}/{genus}/buscoReadsAssemblyHifi.txt")
	threads: 10
	conda: "envs/minimap.yaml"
	shell:
		"""
		if [ -s {input.assemblyfasta} ]; then
			python {scriptdir}/ParseBuscoTableMapping.py -d {input.completed} -i {input.assemblyfasta} -o {output.summary} 
			grep -v 'NOT COMPLETE' {input.nucmercontigs} | cut -f1 | sort | uniq > {output.nucmercontiglist} || true
		else
			touch {output.nucmercontiglist} 
			touch {output.summary}
		fi
		cut -f1 {output.summary} | sort | uniq | grep -v '^#' > {output.buscocontiglist} || true
		cat {output.buscocontiglist} {output.nucmercontiglist} | sort | uniq > {output.contiglist}
		seqtk subseq {input.assemblyfasta} {output.contiglist} > {output.fasta}
		minimap2 -x map-pb -t {threads} {output.fasta} {input.krakenfa}  > {output.paffile}
		python {scriptdir}/PafAlignment.py -p {output.paffile} -o {output.mapping} -r {output.reads}
		comm -12 <(sort {input.unmapped}) <(cut -f2 {output.reads} | tr ',' '\n' | sort | uniq) > {output.reads_unmapped}
		comm -12 <(sort {input.unmapped}) <(sort {input.readfile}) > {output.busco_assembly}
		cat {output.reads_unmapped} {output.busco_assembly} | sort | uniq > {output.busco_assembly_hifi}
		seqtk subseq {input.krakenfa} {output.busco_assembly_hifi} > {output.readsfasta}
		"""

rule RunBuscoReads:
	"""
	Detect number of BUSCO genes per contig
	"""
	input:
		circgenome = "{workingdirectory}/{genus}/{genus}.reads2assemble.fa",
		donetaxon = "{workingdirectory}/taxdownload.done.txt"
	params:
		buscodir = directory("{workingdirectory}/{genus}/buscoReads"),
		genus = "{genus}",
		workingdirectory = "{workingdirectory}",
		taxnames = expand("{datadir}/taxonomy/names.dmp",datadir=config["datadir"]),
		taxnodes = expand("{datadir}/taxonomy/nodes.dmp",datadir=config["datadir"])
	output:
		renamedfa = temporary("{workingdirectory}/{genus}/kraken.renamed.fa"),
		convtable = temporary("{workingdirectory}/{genus}/kraken.convtable.txt"),
		buscodbs = temporary("{workingdirectory}/{genus}/info_dbs_reads.txt"),
		buscoini = temporary("{workingdirectory}/{genus}/config_busco_reads.ini"),
		completed = temporary("{workingdirectory}/{genus}/buscoReads/done.txt"),
		readfile = temporary("{workingdirectory}/{genus}/buscoReads.txt")
	conda: "envs/busco.yaml"
	threads:
		10
	shell:
		"""
		if [ -s {input.circgenome} ]; then
			linecount=$(grep -c '>' < {input.circgenome})
			if [ $linecount -le 100000 ]; then
				python {scriptdir}/RenameFastaHeader.py -i {input.circgenome} -o {output.convtable} > {output.renamedfa}
				busco --list-datasets > {output.buscodbs}
				python {scriptdir}/BuscoConfig.py -na {params.taxnames} -no {params.taxnodes} -f {output.renamedfa} -d {params.buscodir} -dl {datadir}/busco_data/ -c {threads} -db {output.buscodbs} -o {output.buscoini}
				busco --config {output.buscoini} -f || true
				touch {output.completed}
				python {scriptdir}/ParseBuscoTableMappingRead.py -d {output.completed} -c {output.convtable} -o {output.readfile}
			else
				touch {output.renamedfa} {output.convtable} {output.buscodbs} {output.buscoini} {output.readfile}
			fi
		else 
			touch {output.renamedfa}
			touch {output.convtable}
			touch {output.buscodbs}
			touch {output.buscoini}
			touch {output.readfile}
		fi
		touch {output.completed}
		"""

rule DrawCircos:
	"""
	Draw circos plot of re-assembly
	"""
	input:
		assemblyfasta = "{workingdirectory}/{genus}/hifiasm/hifiasm.p_ctg.fasta",
		contiglist = "{workingdirectory}/{genus}/{genus}.Assembly.contigs.txt",
		completed = "{workingdirectory}/{genus}/buscoAssembly/done.txt"
	params:
		dirname = "{workingdirectory}/{genus}/"
	output:
		karyo = temporary("{workingdirectory}/{genus}/circos.karyo"),
		cdsfile = temporary("{workingdirectory}/{genus}/busco.cds.dat"),
		linkfile = temporary("{workingdirectory}/{genus}/links.dat"),
		conffile = "{workingdirectory}/{genus}/circos.conf",
		figure = "{workingdirectory}/{genus}/circos.png"
	conda: "envs/circos.yaml"
	shell:
		"""
		if [ -s {input.contiglist} ]; then
			linecount=$(wc -l < {input.contiglist})
			if [ $linecount -le 200 ]; then
				python {scriptdir}/input_circos.py -f {input.assemblyfasta} -c {input.contiglist} -b {input.completed} -k {output.karyo} -d {output.cdsfile} -l {output.linkfile}
				python {scriptdir}/config_circos.py -k {output.karyo} -d {output.cdsfile} -l {output.linkfile} > {output.conffile}
				circos -conf {output.conffile} -outputdir {params.dirname}
			else
				touch {output.karyo}
				touch {output.cdsfile}
				touch {output.linkfile}
				touch {output.conffile}
				touch {output.figure}
			fi
		else
			touch {output.karyo}
			touch {output.cdsfile}
			touch {output.linkfile}
			touch {output.conffile}
			touch {output.figure}
		fi
		"""
		
def aggregate_assemblies(wildcards):
	checkpoint_output=checkpoints.GetGenera.get(**wildcards).output[0]
	return expand ("{workingdirectory}/{genus}/{genus}.finalassembly.fa", workingdirectory=config["workingdirectory"], genus=glob_wildcards(os.path.join(checkpoint_output, 'genus.{genus}.txt')).genus)

rule concatenate_asm:
	input:
		aggregate_assemblies
	output:
		"{workingdirectory}/final_assembly.fa"
	shell:
		"""
		if [ -n "{input}" ]
		then
			cat {input} > {output}
		else
			touch {output}
		fi
		"""

def aggregate_readsets(wildcards):
	checkpoint_output=checkpoints.GetGenera.get(**wildcards).output[0]
	return expand ("{workingdirectory}/{genus}/{genus}.final_reads.fa", workingdirectory=config["workingdirectory"], genus=glob_wildcards(os.path.join(checkpoint_output, 'genus.{genus}.txt')).genus)

rule concatenate_reads:
	input:
		aggregate_readsets
	output:
		"{workingdirectory}/final_reads_removal.fa"
	shell:
		"""
		if [ -n "{input}" ]
		then
			cat {input} > {output}
		else
			touch {output}
		fi
		"""

'''
def aggregate_readsetslist(wildcards):
	checkpoint_output=checkpoints.GetGenera.get(**wildcards).output[0]
	return expand ("{workingdirectory}/{genus}/{genus}.readsids.txt", workingdirectory=config["workingdirectory"], genus=glob_wildcards(os.path.join(checkpoint_output, 'genus.{genus}.txt')).genus)

rule concatenate_readlist:
	input:
		aggregate_readsetslist
	output:
		cont = "{workingdirectory}/final_reads_removal.txt",
		target = "{workingdirectory}/final_reads_target.fa.gz"
	shell:
		"""
		if [ -n "{input}" ]
		then
			cat {input} > {output.cont}
			if [[ {reads} == *gz ]] 
			then
				zcat {reads} | paste - - - - | grep -v -F -f {output.cont} | tr "\t" "\n" | gzip > {output.target}
			else
				cat {reads} | grep -v -F -f {output.cont} | gzip > {output.target}
			fi
		else
			touch {output.cont}
			cp {reads} {output.target}
		fi
		"""
'''

def aggregate_readsets_putative(wildcards):
	checkpoint_output=checkpoints.GetGenera.get(**wildcards).output[0]
	return expand ("{workingdirectory}/{genus}/{genus}.putative_reads.fa", workingdirectory=config["workingdirectory"], genus=glob_wildcards(os.path.join(checkpoint_output, 'genus.{genus}.txt')).genus)

rule concatenate_reads_putative:
	input:
		aggregate_readsets_putative
	output:
		"{workingdirectory}/putative_reads_removal.fa"
	shell:
		"""
		if [ -n "{input}" ]
		then
			cat {input} > {output}
		else
			touch {output}
		fi
		"""

def aggregate_figures(wildcards):
	checkpoint_output=checkpoints.GetGenera.get(**wildcards).output[0]
	return expand ("{workingdirectory}/{genus}/circos.png", workingdirectory=config["workingdirectory"], genus=glob_wildcards(os.path.join(checkpoint_output, 'genus.{genus}.txt')).genus)

rule concatenate_figures:
	input:
		aggregate_figures
	output:
		temporary("{workingdirectory}/figures_done.txt")
	shell:
		"touch {output}"

rule create_report:
	input:
		finalrem = "{workingdirectory}/final_reads_removal.fa",
		krakenout = "{workingdirectory}/kraken.output",
		putrem = "{workingdirectory}/putative_reads_removal.fa",
		figs = "{workingdirectory}/figures_done.txt"
	params:
		datadir = expand("{datadir}/genera/",datadir=config["datadir"])
	output:
		rep = "{workingdirectory}/{shortname}.report.pdf"
	conda: "envs/fpdf.yaml"
	shell:
		"""
		python {scriptdir}/ReportFile.py -o {output.rep} -r {input.finalrem} -d {params.datadir}
		gzip {input.krakenout}
		"""
