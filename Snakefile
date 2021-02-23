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
#envsdir = config["envsdir"]
sciname_goi = config["sci_name"]
SSUHMMfile = config["SSUHMMfile"]
genome = config["genome"]

rule all:
	input:
		expand("{pwd}/{name}.ProkSSU.reduced.fa",pwd=config["workingdirectory"], name=config["shortname"]),
		expand("{pwd}/{name}.ProkSSU.reduced.SILVA.genus.txt",pwd=config["workingdirectory"], name=config["shortname"]),		
		expand("{pwd}/kraken.tax.masked.ffn",pwd=config["workingdirectory"]),
		expand("{pwd}/kraken.output",pwd=config["workingdirectory"]),
		expand("{pwd}/final_assembly.fa",pwd=config["workingdirectory"]),
		expand("{pwd}/final_reads_removal.fa",pwd=config["workingdirectory"]),
		expand("{pwd}/putative_reads_removal.fa",pwd=config["workingdirectory"]),
		expand("{pwd}/{name}.report.pdf",pwd=config["workingdirectory"], name=config["shortname"])

rule HMMscan_SSU:
	"""
	Run HMMscan with prokaryotic+viral HMM (RF00177+RF01959)
	"""
	output:
		dom = "{workingdirectory}/{shortname}.ProkSSU.domout", 
		log = "{workingdirectory}/{shortname}.HMMscan.log"
	threads: 10
	conda: "envs/hmmer.yaml"
	shell:
		"""
		nhmmer --cpu {threads} --noali --tblout {output.dom} -o {output.log} {SSUHMMfile} {genome}
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
		readslist = "{workingdirectory}/{shortname}.ProkSSU.readslist",
		readslistmicro = "{workingdirectory}/{shortname}.ProkSSU.microsporidia.readslist"
	shell:
		"""
		python {scriptdir}/GetReadsSSU_nhmmer.py -i {input.dom} | grep -v 'RF02542.afa' > {output.readsinfo} || true
		python {scriptdir}/GetReadsSSU_nhmmer.py -i {input.dom} | grep 'RF02542.afa' > {output.readsinfomicro} || true
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
		fasta16S = "{workingdirectory}/{shortname}.ProkSSU.reads.fa",
		fasta16SLoci = "{workingdirectory}/{shortname}.ProkSSU.fa",
		fasta16SLociReduced = "{workingdirectory}/{shortname}.ProkSSU.reduced.fa",
		fasta16Smicro = "{workingdirectory}/{shortname}.ProkSSU.reads.microsporidia.fa",
		fasta16SLocimicro = "{workingdirectory}/{shortname}.ProkSSU.microsporidia.fa",
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
		donesilva = "{workingdirectory}/silva_download.done.txt"
	shell:
		"""
		var=$(curl -L https://ftp.arb-silva.de/current/ARB_files/ | grep 'SSURef_opt.arb.gz.md5' | cut -f2 -d '\"')
		curl -R https://ftp.arb-silva.de/current/ARB_files/$var --output {datadir}/silva/$var
		filename=$(basename $var .md5)
		filenameshort=$(basename $filename .gz)
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
		doneorganelles = "{workingdirectory}/organelles_download.done.txt"
	conda:	"envs/cdhit.yaml"
	shell:
		"""
		if [ ! -d {datadir}/organelles ]; then
  			mkdir {datadir}/organelles
		fi
		if [ -s {datadir}/organelles/organelles.lineage.txt ]; then
        	before=$(date -d 'today - 30 days' +%s)
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

rule DownloadNCBITaxonomy:
	"""
	Download current version of NCBI taxonomy
	"""
	input:
		taxdir = directory(expand("{datadir}/taxonomy/",datadir=config["datadir"])),
	output:
		#taxnames = "{datadir}/taxonomy/names.dmp",
		#taxnodes = "{datadir}/taxonomy/nodes.dmp",
		#accessionfile = "{datadir}/taxonomy/nucl_gb.accession2taxid",
		donefile = "{workingdirectory}/taxdownload.done.txt"
	shell:
		"""
		curl -R ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz.md5 --output {input.taxdir}/taxdump.tar.gz.md5
		if [ -f {input.taxdir}/names.dmp ]; then
			if [ {input.taxdir}/taxdump.tar.gz.md5 -nt {input.taxdir}/names.dmp ]; then
				curl -R ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz --output taxdump.tar.gz
				tar -xvf taxdump.tar.gz
				mv names.dmp {input.taxdir}/names.dmp
				mv nodes.dmp {input.taxdir}/nodes.dmp
				rm division.dmp gencode.dmp citations.dmp delnodes.dmp merged.dmp readme.txt gc.prt taxdump.*
			fi
		else
			curl -R ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz --output taxdump.tar.gz
			tar -xvf taxdump.tar.gz
			mv names.dmp {input.taxdir}/names.dmp
			mv nodes.dmp {input.taxdir}/nodes.dmp
			rm division.dmp gencode.dmp citations.dmp delnodes.dmp merged.dmp readme.txt gc.prt taxdump.*
		fi
		curl -R ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz.md5 --output nucl_gb.accession2taxid.gz.md5
		if [ -f {input.taxdir}/nucl_gb.accession2taxid ]; then
			if [ nucl_gb.accession2taxid.gz.md5 -nt {input.taxdir}/nucl_gb.accession2taxid.gz.md5 ]; then
				mv nucl_gb.accession2taxid.gz.md5 {input.taxdir}/nucl_gb.accession2taxid.gz.md5
				curl -R ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz --output nucl_gb.accession2taxid.gz
				curl -R ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz --output nucl_wgs.accession2taxid.gz
				gunzip nucl_gb.accession2taxid.gz nucl_wgs.accession2taxid.gz
				mv nucl_wgs.accession2taxid {input.taxdir}/nucl_wgs.accession2taxid
				mv nucl_gb.accession2taxid {input.taxdir}/nucl_gb.accession2taxid
			else
				rm nucl_gb.accession2taxid.gz.md5
			fi
		else
			curl -R ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz --output nucl_gb.accession2taxid.gz
			curl -R ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz --output nucl_wgs.accession2taxid.gz
			gunzip nucl_gb.accession2taxid.gz nucl_wgs.accession2taxid.gz
			mv nucl_wgs.accession2taxid {input.taxdir}/nucl_wgs.accession2taxid
			mv nucl_gb.accession2taxid {input.taxdir}/nucl_gb.accession2taxid			
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
		taxnames = expand("{datadir}/taxonomy/names.dmp",datadir=config["datadir"]),
		taxnodes = expand("{datadir}/taxonomy/nodes.dmp",datadir=config["datadir"]),
		donesilva = "{workingdirectory}/silva_download.done.txt"
	output:
		SILVA_output_embl = "{workingdirectory}/{shortname}.ProkSSU.reduced.SILVA.embl.csv",
		SILVA_output_silva = "{workingdirectory}/{shortname}.ProkSSU.reduced.SILVA.silva.csv",
		SILVA_output = "{workingdirectory}/{shortname}.ProkSSU.reduced.SILVA.csv",
		SILVA_tax = "{workingdirectory}/{shortname}.ProkSSU.reduced.SILVA.tax",
		blastout = "{workingdirectory}/{shortname}.ProkSSU.reduced.microsporidia.blast.txt",
		blastgenus = "{workingdirectory}/{shortname}.ProkSSU.reduced.microsporidia.genus.txt",
		SILVA16Sgenus = "{workingdirectory}/{shortname}.ProkSSU.reduced.SILVA.genus.txt"
	conda: "envs/sina.yaml"
	threads: 10
	shell:
		"""
		sina -i {input.fasta16SLociReduced} -o {output.SILVA_output_embl} --db {datadir}/silva/SILVA_SSURef.arb --search --search-min-sim 0.9 -p {threads} --lca-fields tax_embl_ebi_ena --outtype csv
		sina -i {input.fasta16SLociReduced} -o {output.SILVA_output_silva} --db {datadir}/silva/SILVA_SSURef.arb --search --search-min-sim 0.9 -p {threads} --lca-fields tax_slv --outtype csv
		cat {output.SILVA_output_embl} {output.SILVA_output_silva} > {output.SILVA_output}
		cut -f1,8 -d',' {output.SILVA_output} | tr ',' '\t' > {output.SILVA_tax}
		cut -f2 {output.SILVA_tax} | grep -v 'lca_tax_embl_ebi_ena' | grep -v 'lca_tax_slv' | sort | uniq > {output.SILVA16Sgenus} && [[ -s {output.SILVA16Sgenus} ]]
		if [ -s {input.fasta16SLociReducedmicro} ]; then
			blastn -db {datadir}/silva/MicrosporidiaSSU_NCBI -query {input.fasta16SLociReducedmicro} -out {output.blastout} -outfmt 6
			python {scriptdir}/ParseBlastLineage.py -b {output.blastout} -na {input.taxnames} -no {input.taxnodes} > {output.blastgenus}
			cat {output.blastgenus} >> {output.SILVA16Sgenus}
		else
			touch {output.blastout}
			touch {output.blastgenus}
		fi
		"""

checkpoint GetGenera:
	"""
	Get genera which were detected in SILVA DB 16 screen
	"""
	input:
		#SILVA16Sgenus = "{workingdirectory}/{shortname}.ProkSSU.reduced.SILVA.genus.txt"
		SILVA16Sgenus = expand("{pwd}/{name}.ProkSSU.reduced.SILVA.genus.txt",pwd=config["workingdirectory"], name=config["shortname"]),
		#screenfile = "{workingdirectory}/refseq202.screen",
		#asminfo = expand("{datadir}/mash/assembly_summary_refseq.txt",datadir=config["datadir"]),
		taxnames = expand("{datadir}/taxonomy/names.dmp",datadir=config["datadir"]),
		taxnodes = expand("{datadir}/taxonomy/nodes.dmp",datadir=config["datadir"])
	output:
		generadir = directory("{workingdirectory}/genera")
	shell:
		"""
		mkdir {output.generadir}
		python {scriptdir}/DetermineGenera.py -i {input.SILVA16Sgenus} -t family -na {input.taxnames} -no {input.taxnodes} -od {output} -suf SSU.genera_taxonomy.txt -g '{sciname_goi}'
		while read p
		do
			echo "Eukaryota" > {output.generadir}/genus.$p.txt
		done < {output.generadir}/euk.SSU.genera_taxonomy.txt
		while read p
		do
			echo "Bacteria/Archaea" > {output.generadir}/genus.$p.txt
			#touch  {output.generadir}/genus.$p.txt
		done < {output.generadir}/prok.SSU.genera_taxonomy.txt
		rm {output.generadir}/euk.SSU.genera_taxonomy.txt
		rm {output.generadir}/prok.SSU.genera_taxonomy.txt
		"""

rule DownloadRefSeqGenus:
	"""
	Download RefSeq genomes (per species) of selected genera from 16S screen
	"""
	input:
		generafiles = "{workingdirectory}/genera/genus.{genus}.txt",
		doneorganelles = "{workingdirectory}/organelles_download.done.txt"
	params:
		taxname = "{genus}"
	output:
		#novel_pwd = directory("{datadir}/genera/{genus}"),
		#refseqlog = "{datadir}/genera/{genus}/{genus}.refseq.log",
		#downloadlog = "{datadir}/genera/{genus}/{genus}.download.log",
		#refseqdir = directory("{datadir}/genera/{genus}/{genus}.Refseq"),
		#refseqdir_orig = temp("{datadir}/genera/{genus}/RefSeq.{genus}.zip"),
		krakenffnall = "{workingdirectory}/genera/{genus}.kraken.tax.ffn",
		orglist = "{workingdirectory}/genera/{genus}.organelles.list",
		orgfasta = "{workingdirectory}/genera/{genus}.organelles.ffn",
		donefile = "{workingdirectory}/{genus}.refseqdownload.done.txt"
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
					unzip -d {datadir}/genera/{params.taxname}/{params.taxname}.Refseq {datadir}/genera/{params.taxname}/RefSeq.{params.taxname}.zip
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
				unzip -d {datadir}/genera/{params.taxname}/{params.taxname}.Refseq {datadir}/genera/{params.taxname}/RefSeq.{params.taxname}.zip
				python {scriptdir}/AddTaxIDKraken.py -d {datadir}/genera/{params.taxname}/{params.taxname}.Refseq -o {datadir}/genera/{params.taxname}.kraken.tax.ffn
			else
				touch {datadir}/genera/{params.taxname}.kraken.tax.ffn
			fi
		fi
		if grep -q Eukaryota {input.generafiles}; then
			grep {params.taxname} {datadir}/organelles/organelles.lineage.txt > {output.orglist} || true
			python {scriptdir}/FastaSelect.py -f {datadir}/organelles/organelles.fna -l {output.orglist} -o {output.orgfasta}
		else
			touch {output.orglist}
			touch {output.orgfasta}
		fi
		cat {datadir}/genera/{params.taxname}.kraken.tax.ffn {output.orgfasta} > {output.krakenffnall}
		touch {output.donefile}
		"""

def aggregate_kraken(wildcards):
	checkpoint_output=checkpoints.GetGenera.get(**wildcards).output[0]
	return expand ("{workingdirectory}/genera/{genus}.kraken.tax.ffn", workingdirectory=config["workingdirectory"], genus=glob_wildcards(os.path.join(checkpoint_output, 'genus.{genus}.txt')).genus)

rule concatenate_kraken_input:
	input:
		aggregate_kraken
	output:
		"{workingdirectory}/kraken.tax.ffn"
	shell:
		"cat {input} > {output}"


rule DownloadGenusRel:
	"""
	Download assemblies of closely related species to species of interest
	"""
	input:
		taxnames = expand("{datadir}/taxonomy/names.dmp",datadir=config["datadir"]),
		taxnodes = expand("{datadir}/taxonomy/nodes.dmp",datadir=config["datadir"])
	output:
		novel_pwd = directory("{workingdirectory}/relatives/"),
		refseqlog = "{workingdirectory}/relatives/relatives.refseq.log",
		refseqdir = directory("{workingdirectory}/relatives/relatives.Refseq"),
		refseqdir_orig = "{workingdirectory}/relatives/RefSeq.relatives.zip",
		krakenffnrel = "{workingdirectory}/relatives/relatives.kraken.tax.ffn"
	shell:
		"""
		python {scriptdir}/FetchGenomesRefSeqRelatives.py --taxname '{sciname_goi}' --dir {output.novel_pwd} -na {input.taxnames} -no {input.taxnodes} > {output.refseqlog}
		unzip -d {output.refseqdir} {output.refseqdir_orig}
		python {scriptdir}/AddTaxIDKraken.py -d {output.refseqdir} -o {output.krakenffnrel}
		"""

checkpoint SplitFasta:
	"""
	Split downloaded assemblies fasta file depending on the number of cores
	"""
	input:
		krakenffnall = "{workingdirectory}/kraken.tax.ffn"
	output:
		splitdir = directory("{workingdirectory}/split_fasta/")
	shell:
		"""
		python {scriptdir}/FastaSplit.py -f {input.krakenffnall} -s 5000 -o {output.splitdir}
		"""

rule doMasking:
	"""
	Rule to mask repetitive regions in fasta file
	"""
	input:
		fastafile = "{workingdirectory}/split_fasta/kraken.tax.{num}.fa"
	output:
		maskedfile = "{workingdirectory}/split_fasta/kraken.tax.{num}.masked.fa"
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
		cat {input} > {output}
		"""

checkpoint SplitFastaRel:
	"""
	Split downloaded assemblies fasta file depending on the number of cores
	"""
	input:
		krakenffnall = "{workingdirectory}/relatives/relatives.kraken.tax.ffn"
	output:
		splitdir = directory("{workingdirectory}/split_fasta_rel/")
	shell:
		"""
		python {scriptdir}/FastaSplit.py -f {input.krakenffnall} -s 5000 -o {output.splitdir}
		"""

rule doMaskingRel:
	"""
	Rule to mask repetitive regions in fasta file
	"""
	input:
		fastafile = "{workingdirectory}/split_fasta_rel/kraken.tax.{num}.fa"
	output:
		maskedfile = "{workingdirectory}/split_fasta_rel/kraken.tax.{num}.masked.fa"
	conda: "envs/kraken.yaml"
	threads: 1
	shell:
		"""
		dustmasker -in {input.fastafile} -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' > {output.maskedfile}
		"""

def aggregate_masking_relatives(wildcards):
	checkpoint_output=checkpoints.SplitFastaRel.get(**wildcards).output[0]
	return expand ("{workingdirectory}/split_fasta_rel/kraken.tax.{num}.masked.fa", workingdirectory=config["workingdirectory"], num=glob_wildcards(os.path.join(checkpoint_output, 'kraken.tax.{num}.fa')).num)

rule concatenate_masking_relatives:
	input:
		aggregate_masking_relatives
	output:
		"{workingdirectory}/kraken.relatives.masked.ffn"
	shell:
		"""
		cat {input} > {output}
		"""

rule CreateKrakenDB:
	"""
	Create Kraken DB for all downloaded refseq genomes
	"""
	input:
		donefile = "{workingdirectory}/taxdownload.done.txt",
		krakenffnall = "{workingdirectory}/kraken.tax.masked.ffn",
		krakenffnrel = "{workingdirectory}/kraken.relatives.masked.ffn",
		splitdirrel = directory("{workingdirectory}/split_fasta_rel/"),
		splitdir = directory("{workingdirectory}/split_fasta/"),
		krakenfasta = "{workingdirectory}/kraken.tax.ffn",
		krakenrelfasta = "{workingdirectory}/relatives/relatives.kraken.tax.ffn"
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
		fi
		rm -r {input.splitdirrel}
		rm -r {input.splitdir}
		rm {input.krakenfasta}
		rm {input.krakenrelfasta}
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
		krakenfa = "{workingdirectory}/{genus}/kraken.fa"
	output:
		paffile = "{workingdirectory}/{genus}/{genus}.paf",
		mapping = "{workingdirectory}/{genus}/{genus}.ctgs",
		contiglist = "{workingdirectory}/{genus}/{genus}.ctgs.list",
		reads = "{workingdirectory}/{genus}/{genus}.reads",
		fasta = "{workingdirectory}/{genus}/{genus}.ctgs.fa"
	threads: 10
	conda: "envs/minimap.yaml"
	shell:
		"""
		minimap2 -x map-pb -t {threads} {genome} {input.krakenfa}  > {output.paffile}
		python {scriptdir}/PafAlignment.py -p {output.paffile} -o {output.mapping} -r {output.reads}
		grep -v 'NOT COMPLETE' {output.mapping} | cut -f1 | sort | uniq > {output.contiglist} || true
		seqtk subseq {genome} {output.contiglist} > {output.fasta}
		"""

rule RunBusco:
	"""
	Detect number of BUSCO genes per contig
	"""
	input:
		circgenome = "{workingdirectory}/{genus}/{genus}.ctgs.fa",
		taxnames = expand("{datadir}/taxonomy/names.dmp",datadir=config["datadir"]),
		taxnodes = expand("{datadir}/taxonomy/nodes.dmp",datadir=config["datadir"]),
	params:
		buscodir = directory("{workingdirectory}/{genus}/busco")
	output:
		buscodbs = "{workingdirectory}/{genus}/info_dbs.txt",
		buscoini = "{workingdirectory}/{genus}/config_busco.ini",
		#proteins = "{workingdirectory}/{genus}/busco/busco/prodigal_output/predicted_genes/predicted.faa",
		completed = "{workingdirectory}/{genus}/busco/done.txt"
	conda: "envs/busco.yaml"
	threads:
		10
	shell:
		"""
		if [ -s {input.circgenome} ]; then
			busco --list-datasets > {output.buscodbs}
			python {scriptdir}/BuscoConfig.py -na {input.taxnames} -no {input.taxnodes} -f {input.circgenome} -d {params.buscodir} -dl {datadir}/busco_data/ -c {threads} -db {output.buscodbs} -o {output.buscoini}
			busco --config {output.buscoini} -f
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
		completed = "{workingdirectory}/{genus}/nucmer_contigs.done.txt",
		nucmerdelta = "{workingdirectory}/{genus}/{genus}_vs_contigs.delta",
        nucmercoords = "{workingdirectory}/{genus}/{genus}_vs_contigs.coords.txt",
		nucmercontigs = "{workingdirectory}/{genus}/{genus}_vs_contigs.overview.txt"
	conda: "envs/nucmer.yaml"
	shell:
		"""
		if [ -s {input.circgenome} ]; then
			nucmer --delta {output.nucmerdelta} {input.circgenome} {input.refseqmasked}
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
		contigsid = "{workingdirectory}/{genus}/{genus}.ids.txt",
		readids = "{workingdirectory}/{genus}/{genus}.readsids.txt",
		finalreads = "{workingdirectory}/{genus}/{genus}.final_reads.fa",
		nucmercontiglist = "{workingdirectory}/{genus}/{genus}.nucmer.contigs.txt",
		buscocontiglist = "{workingdirectory}/{genus}/{genus}.busco.contigs.txt",
		unmapped = "{workingdirectory}/{genus}/{genus}.unmapped.reads",
		unmappedfa = "{workingdirectory}/{genus}/{genus}.unmapped.fa",
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

rule NucmerRefSeqReads:
	"""
	Alignment all contigs against reference genomes
	"""
	input:
		circgenome = "{workingdirectory}/{genus}/{genus}.unmapped.fa",
		refseqmasked = "{workingdirectory}/genera/{genus}.kraken.tax.ffn",
		krakenfa = "{workingdirectory}/{genus}/kraken.fa",
	output:
		completed = "{workingdirectory}/{genus}/nucmer_reads.done.txt",
		nucmerdelta = "{workingdirectory}/{genus}/{genus}_vs_reads.delta",
        nucmercoords = "{workingdirectory}/{genus}/{genus}_vs_reads.coords.txt",
		nucmercontigs = "{workingdirectory}/{genus}/{genus}_vs_reads.overview.txt",
		nucmercontiglist = "{workingdirectory}/{genus}/{genus}.nucmer.reads.txt",
		finalreads = "{workingdirectory}/{genus}/{genus}.putative_reads.fa",
	conda: "envs/nucmer.yaml"
	shell:
		"""
		if [ -s {input.circgenome} ]; then
			nucmer --delta {output.nucmerdelta} {input.circgenome} {input.refseqmasked}
			show-coords -c -l -L 100 -r -T {output.nucmerdelta} > {output.nucmercoords}
			python {scriptdir}/ParseNucmer.py -n {output.nucmercoords} -o {output.nucmercontigs}
			grep -v 'NOT COMPLETE' {output.nucmercontigs} | cut -f1 | sort | uniq > {output.nucmercontiglist} || true
			seqtk subseq {input.krakenfa} {output.nucmercontiglist} > {output.finalreads}
		else
			touch {output.nucmerdelta}
			touch {output.nucmercoords}
			touch {output.nucmercontigs}
			touch {output.nucmercontiglist}
			touch {output.finalreads}
		fi
		touch {output.completed}
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
		"cat {input} > {output}"

def aggregate_readsets(wildcards):
	checkpoint_output=checkpoints.GetGenera.get(**wildcards).output[0]
	return expand ("{workingdirectory}/{genus}/{genus}.final_reads.fa", workingdirectory=config["workingdirectory"], genus=glob_wildcards(os.path.join(checkpoint_output, 'genus.{genus}.txt')).genus)

rule concatenate_reads:
	input:
		aggregate_readsets
	output:
		"{workingdirectory}/final_reads_removal.fa"
	shell:
		"cat {input} > {output}"

def aggregate_readsets_putative(wildcards):
	checkpoint_output=checkpoints.GetGenera.get(**wildcards).output[0]
	return expand ("{workingdirectory}/{genus}/{genus}.putative_reads.fa", workingdirectory=config["workingdirectory"], genus=glob_wildcards(os.path.join(checkpoint_output, 'genus.{genus}.txt')).genus)

rule concatenate_reads_putative:
	input:
		aggregate_readsets_putative
	output:
		"{workingdirectory}/putative_reads_removal.fa"
	shell:
		"cat {input} > {output}"


rule create_report:
	input:
		finalrem = "{workingdirectory}/final_reads_removal.fa",
		krakenout = "{workingdirectory}/kraken.output",
		putrem = "{workingdirectory}/putative_reads_removal.fa"
	output:
		rep = "{workingdirectory}/{shortname}.report.pdf"
	conda: "envs/fpdf.yaml"
	shell:
		"""
		python {scriptdir}/ReportFile.py -o {output.rep} -r {input.finalrem}
		gzip {input.krakenout}
		"""
