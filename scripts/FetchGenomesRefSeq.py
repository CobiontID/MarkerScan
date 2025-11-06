from __future__ import division
import time
import argparse
import os
from pathlib import Path

from NCBIApiTools import NcbiApi

parser = argparse.ArgumentParser()
parser.add_argument(
    "--taxname",
    action="store",
    dest="tax",
    type=str,
    help="a genus taxname to download refseq genomes for",
)
parser.add_argument(
    "--dir", action="store", dest="dir", type=str, help="base directory"
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


api_instance = NcbiApi(args.key)

taxname = str(args.tax).split("genus.")[1].split(".")[0]
taxname_orig = taxname
if "_" in taxname:
    taxname = taxname.replace("_", " ")

if not os.path.exists(args.dir):
    os.makedirs(args.dir)

# fetch data from NCBI via 'datasets' of all species from that clade
is_refseq = args.refs.lower() == "yes" if args.refs else False

assemblies = api_instance.get_assemblies_for_taxon(
    taxon=str(taxname), filters_reference_only=is_refseq, page_size=1000
)

SpeciesDictionary = {}
for assembly in assemblies:
    acc = assembly.get("accession")
    date = assembly.get("assembly_info", {}).get("release_date")
    sciname_orig = assembly.get("organism").get("organism_name")
    strainname = (
        assembly.get("organism", {}).get("infraspecific_names", {}).get("strain", None)
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
genomesizes = []
for species in SpeciesDictionary:
    accs.append(SpeciesDictionary[species]["Identifier"])
    genomesizes.append(SpeciesDictionary[species]["GenomeSize"])

if len(accs) > 0:
    nrgfs = int(round(float(Average(genomesizes) / 1000000) * 20))
    print("Genomes for " + str(len(accs)) + " species")
    print("Number of necessary GFs " + str(nrgfs))
    print(f"Download a package for {accs}.")
    print("Begin download of genome data package ...")
    print(str(len(accs)) + " genomes")
    zipfile_name = str(args.dir) + "/" + "/RefSeq." + str(taxname_orig) + ".zip"
    for t in range(0, len(accs), 100):
        accshort = accs[t : t + 100]
        zipfile_name_part = Path(args.dir) / f"/{taxname_orig}.RefSeq.part{str(t)}.zip"
        api_response = api_instance.download_genomes(
            accshort,
            outfile=zipfile_name_part,
        )
        print("Download complete part " + str(t))
        cmd = f"unzip -d {zipfile_name_part.stem} {zipfile_name_part}"
        os.system(cmd)

    cmd = "mkdir " + str(args.dir) + "/" + str(taxname_orig) + ".Refseq"
    os.system(cmd)
    cmd = "mkdir " + str(args.dir) + "/" + str(taxname_orig) + ".Refseq/ncbi_dataset"
    os.system(cmd)
    cmd = (
        "mkdir " + str(args.dir) + "/" + str(taxname_orig) + ".Refseq/ncbi_dataset/data"
    )
    os.system(cmd)
    cmd = (
        "cp -r "
        + str(args.dir)
        + "/"
        + str(taxname_orig)
        + ".RefSeq.part*/ncbi_dataset/data/G* "
        + str(args.dir)
        + "/"
        + str(taxname_orig)
        + ".Refseq/ncbi_dataset/data/"
    )
    os.system(cmd)
    cmd = (
        "cat "
        + str(args.dir)
        + "/"
        + str(taxname_orig)
        + ".RefSeq.part*/ncbi_dataset/data/assembly_data_report.jsonl >> "
        + str(args.dir)
        + "/"
        + str(taxname_orig)
        + ".Refseq/ncbi_dataset/data/assembly_data_report.jsonl"
    )
    os.system(cmd)
    cmd = (
        "cat "
        + str(args.dir)
        + "/"
        + str(taxname_orig)
        + ".RefSeq.part*/ncbi_dataset/data/dataset_catalog.json >> "
        + str(args.dir)
        + "/"
        + str(taxname_orig)
        + ".Refseq/ncbi_dataset/data/dataset_catalog.json"
    )
    os.system(cmd)
    cmd = "rm -r " + str(args.dir) + f"/{str(taxname_orig)}.RefSeq.part*"
    os.system(cmd)
else:
    print("No genomes available")
    cmd = "touch " + str(args.dir) + "/" + str(taxname_orig) + ".download.log"
    os.system(cmd)
