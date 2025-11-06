from __future__ import division

import argparse
import os
import time
from pathlib import Path

from NCBIApiTools import NcbiApi

parser = argparse.ArgumentParser()
parser.add_argument(
    "--taxname",
    action="store",
    dest="tax",
    type=str,
    help="scientific name of species of interest",
)
parser.add_argument(
    "--dir", action="store", dest="dir", type=str, help="base directory"
)
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
    "--refseq", action="store", dest="refs", type=str, help="all or refseq database"
)
parser.add_argument(
    "-k",
    type=str,
    action="store",
    dest="key",
    metavar="NCBI API KEY",
    help="NCBI API Key",
    required=False,
    default=os.environ.get("NCBI_API_KEY"),
)
args = parser.parse_args()


def Average(lst):
    return sum(lst) / len(lst)


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
            if "scientific" in line or "synonym" in line:
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


api_instance = NcbiApi(args.key)

if not os.path.exists(args.dir):
    os.makedirs(args.dir)

taxparents = readNodes(args.nodesfile)
taxnames, namestax = readNames(args.namesfile)

is_refseq = args.refs.lower() == "yes" if args.refs else False

if args.tax in namestax:
    taxid = namestax[args.tax]
    parent = taxparents[taxid]
    parentname = taxnames[parent]
    foundlevel = False
    print(str(taxid) + "\t" + str(parent) + "\t" + parentname)
    while parent != taxparents[parent]:
        time.sleep(1)
        parentname = taxnames[parent]
        parentname_combi = parentname
        if " " in parentname:
            parentname_combi = parentname.replace(" ", "_")
        print(
            str(parent) + "\t" + parentname + "\t" + parentname_combi + "\t" + args.tax
        )
        # fetch data from NCBI via 'datasets' of all species from that clade
        assemblies = api_instance.get_assemblies_for_taxon(
            taxon=str(parent), filters_reference_only=is_refseq, page_size=1000
        )

        print(assemblies)
        i = 0
        if len(assemblies) > 0:
            for assembly in assemblies:
                # print(assembly.org.sci_name)
                sciname_orig = assembly.get("organism").get("organism_name")
                strainname = (
                    assembly.get("organism", {})
                    .get("infraspecific_names", {})
                    .get("strain", None)
                )
                if strainname:
                    if strainname in sciname_orig and "sp" not in sciname_orig:
                        sciname = sciname_orig.replace(strainname, "").strip()
                    else:
                        sciname = sciname_orig
                else:
                    sciname = sciname_orig
                if sciname != args.tax:
                    i = i + 1
        if i > 0:
            break
        parent = taxparents[parent]

    SpeciesDictionary = {}
    for assembly in assemblies:
        acc = assembly.get("accession")
        date = assembly.get("assembly_info", {}).get("release_date")
        sciname_orig = assembly.get("organism").get("organism_name")
        strainname = (
            assembly.get("organism", {})
            .get("infraspecific_names", {})
            .get("strain", None)
        )
        contiguity = int(assembly.get("assembly_stats").get("contig_n50"))
        size = int(assembly.get("assembly_stats").get("total_sequence_length"))

        if strainname:
            if strainname in sciname_orig and "sp" not in sciname_orig:
                sciname = sciname_orig.replace(strainname, "").strip()
            else:
                sciname = sciname_orig
        else:
            sciname = sciname_orig

        if sciname != args.tax:
            print("SECOND STEP")
            print(sciname)
            if sciname not in SpeciesDictionary:
                SpeciesDictionary[sciname] = {}
                SpeciesDictionary[sciname]["ReleaseDate"] = date
                SpeciesDictionary[sciname]["Identifier"] = acc
                SpeciesDictionary[sciname]["GenomeSize"] = size
                SpeciesDictionary[sciname]["N50"] = contiguity
            else:
                novel_submission = time.strptime(date, "%Y-%m-%d")
                old_submission = time.strptime(
                    SpeciesDictionary[sciname]["ReleaseDate"], "%Y-%m-%d"
                )
                if (
                    novel_submission > old_submission
                    or contiguity > SpeciesDictionary[sciname]["N50"]
                ):
                    SpeciesDictionary[sciname]["ReleaseDate"] = date
                    SpeciesDictionary[sciname]["Identifier"] = acc
                    SpeciesDictionary[sciname]["GenomeSize"] = size
                    SpeciesDictionary[sciname]["N50"] = contiguity

    accs = []
    for species in SpeciesDictionary:
        accs.append(SpeciesDictionary[species]["Identifier"])
    print(f"Download a package for {accs}.")
    print("Begin download of genome data package ...")
    for t in range(0, len(accs), 3):
        zipfile_name_part = Path(str(args.dir)) / f"relatives.RefSeq.part{str(t)}.zip"
        accshort = accs[t : t + 3]
        api_response = api_instance.download_genomes(
            accshort,
            outfile=zipfile_name_part,
        )
        print("Download complete part " + str(t))
        cmd = f"unzip -d {zipfile_name_part.stem} {zipfile_name_part}"
        os.system(cmd)

    cmd = "mkdir " + str(args.dir) + "/relatives.Refseq"
    os.system(cmd)
    cmd = "mkdir " + str(args.dir) + "/relatives.Refseq/ncbi_dataset"
    os.system(cmd)
    cmd = "mkdir " + str(args.dir) + "/relatives.Refseq/ncbi_dataset/data"
    os.system(cmd)
    cmd = (
        "cp -r "
        + str(args.dir)
        + "/relatives.RefSeq.part*/ncbi_dataset/data/G* "
        + str(args.dir)
        + "/relatives.Refseq/ncbi_dataset/data/"
    )
    os.system(cmd)
    cmd = (
        "cat "
        + str(args.dir)
        + "/relatives.RefSeq.part*/ncbi_dataset/data/assembly_data_report.jsonl >> "
        + str(args.dir)
        + "/relatives.Refseq/ncbi_dataset/data/assembly_data_report.jsonl"
    )
    os.system(cmd)
    cmd = (
        "cat "
        + str(args.dir)
        + "/relatives.RefSeq.part*/ncbi_dataset/data/dataset_catalog.json >> "
        + str(args.dir)
        + "/relatives.Refseq/ncbi_dataset/data/dataset_catalog.json"
    )
    os.system(cmd)
    cmd = "rm -r " + str(args.dir) + "/" + "relatives.RefSeq.part*"
    os.system(cmd)
    cmd = (
        "cat "
        + str(args.dir)
        + "/relatives.Refseq/ncbi_dataset/data/*/*fna > "
        + str(args.dir)
        + "/relatives.Refseq/acc.fasta"
    )
    os.system(cmd)
    cmd = (
        "dustmasker -in "
        + str(args.dir)
        + "/relatives.Refseq/acc.fasta -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' > "
        + str(args.dir)
        + "/relatives.Refseq/masked.fna"
    )
    os.system(cmd)
    cmd = "rm " + str(args.dir) + "/relatives.Refseq/acc.fasta "
    os.system(cmd)
else:
    print("No taxonomic name was not found")
