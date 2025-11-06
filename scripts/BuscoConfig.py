from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument(
    "-na",
    type=str,
    action="store",
    dest="namesfile",
    metavar="NAMES",
    help="NCBI names.dmp",
)
parser.add_argument(
    "-no",
    type=str,
    action="store",
    dest="nodesfile",
    metavar="NODES",
    help="NCBI nodes.dmp",
)
parser.add_argument(
    "-f",
    type=str,
    action="store",
    dest="genome",
    metavar="GENOME FASTA",
    help="fasta genome assembly file",
)
parser.add_argument(
    "-d",
    type=str,
    action="store",
    dest="dir",
    metavar="WORKDIR",
    help="define working directory for busco",
)
parser.add_argument(
    "-db", type=str, action="store", dest="db", help="define available dbs file"
)
parser.add_argument(
    "-dl",
    type=str,
    action="store",
    dest="download",
    help="define directory to store busco dbs",
)
parser.add_argument("-c", type=int, action="store", dest="cpu", help="define cpus")
parser.add_argument(
    "-o",
    type=str,
    action="store",
    dest="out",
    metavar="OUTFILE",
    help="define configfile name",
)
parser.add_argument("--version", action="version", version="%(prog)s 1.0")
args = parser.parse_args()


def readNames(names_tax_file):
    """
    input:
    - name.dmp (NCBI Taxonomy)
    output:
    - dictionary of form {node: name}
    - dictionary of form {sci name: node}
    """
    tax_names = {}
    tax_names_reverse = {}
    with open(names_tax_file, "r") as nodes_tax:
        for line in nodes_tax:
            node = [field.strip() for field in line.split("|")]
            if "scientific" in line:
                tax_names[node[1]] = node[0]
                tax_names_reverse[node[0]] = node[1]
    return tax_names_reverse, tax_names


def readNodes(nodes_tax_file):
    """
    input:
    - nodes.dmp (NCBI Taxonomy)
    output:
    - dictionary of form {parent: node}
    - dictionary of form {node: type}
    """

    tax_nodes = {}
    tax_types = {}
    with open(nodes_tax_file, "r") as nodes_tax:
        for line in nodes_tax:
            node = [field.strip() for field in line.split("|")]  # make list of line
            tax_nodes[node[0]] = node[1]  # couple node with parent
            tax_types[node[0]] = node[2]  # couple node with rank
    return tax_nodes


taxparents = readNodes(args.nodesfile)
taxnames, namestax = readNames(args.namesfile)

genus = args.out.split("/config")[0].split("/")[-1]
if "_" in genus:
    genus = genus.replace("_", " ")

busco_dbs = []
busco_short = []
m = open(args.db, "r")
for line in m:
    line = line.strip()
    if "db" in line:
        # if 'eukaryota' in line:
        #    break
        dbname = line.split(" ")[-1].split("_odb")[0]
        # print(dbname)
        dbname2 = line.split(" ")[-1].split("_")[0]
        busco_short.append(dbname2)
        busco_dbs.append(dbname)

buscoset = "Bacteria"
if genus in namestax:
    taxid = namestax[genus]
    parent = taxparents[taxid]
    while parent != taxparents[parent]:
        print(taxnames[parent])
        if taxnames[parent].lower() in busco_short:
            buscoset = taxnames[parent]
            print(buscoset)
            if taxnames[parent].lower() + "_phylum" in busco_dbs:
                buscoset = taxnames[parent].lower() + "_phylum"
            break
        parent = taxparents[parent]

print(genus + "\t" + buscoset)

condadir = os.environ["CONDA_DEFAULT_ENV"]

config_content = f"""[busco_run]
# Input file
in = {args.genome}
# Run name, used in output files and folder
out = busco
# Where to store the output directory
out_path = {args.dir}
# Path to the BUSCO dataset
lineage_dataset = {buscoset.lower()}
# Which mode to run (genome / proteins / transcriptome)
mode = genome
# How many threads to use for multithreaded steps
cpu = {args.cpu}
# Force rewrite if files already exist (True/False)
;force = False
# Local destination path for downloaded lineage datasets
download_path = {args.download}
;[tblastn]
;path = {condadir}/bin/
;command = tblastn
;[makeblastdb]
;path = {condadir}/bin/
;command = makeblastdb
;[augustus]
;path = {condadir}/bin/
;command = augustus
;[etraining]
;path = {condadir}/bin/
;command = etraining
;[gff2gbSmallDNA.pl]
;path = {condadir}/bin/
;command = gff2gbSmallDNA.pl
;[new_species.pl]
;path = {condadir}/bin/
;command = new_species.pl
;[optimize_augustus.pl]
;path = {condadir}/bin/
;command = optimize_augustus.pl
;[hmmsearch]
;path = {condadir}/bin/
;command = hmmsearch
;[sepp]
;path = {condadir}/bin/
;command = run_sepp.py
;[prodigal]
;path = {condadir}/bin/
;command = prodigal
"""

with open(args.out, "w") as l:
    l.write(config_content)
