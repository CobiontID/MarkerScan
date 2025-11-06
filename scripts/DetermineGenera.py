from __future__ import division

import argparse
import os
from NCBIApiTools import NcbiApi

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    type=str,
    action="store",
    dest="tax",
    metavar="TAX",
    help="define tax genus file",
)
parser.add_argument(
    "-t",
    type=str,
    action="store",
    dest="type",
    metavar="TYPE",
    help="define ranking type",
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
    "-suf",
    type=str,
    action="store",
    dest="suffix",
    metavar="SUFFIX",
    help="suffix outputfile",
)
parser.add_argument(
    "-od",
    type=str,
    action="store",
    dest="outdir",
    metavar="OUTDIR",
    help="output directory",
)
parser.add_argument(
    "-g",
    type=str,
    action="store",
    dest="spoi",
    metavar="SPOI",
    help="species of interest",
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
    tax_names_sci = {}
    multiple_names = {}
    with open(names_tax_file, "r") as nodes_tax:
        for line in nodes_tax:
            node = [field.strip() for field in line.split("|")]
            if "scientific" in line or "synonym" in line:
                if "synonym" in line and "Bacteria" == node[1]:
                    # apparently there exists as class of walking sticks called Bacteria Latreilla (629395), that have as synonym name Bacteria
                    print("wrong Bacteria")
                else:
                    if node[1] in tax_names:
                        orig = tax_names[node[1]]
                        if not node[1] in multiple_names:
                            multiple_names[node[1]] = []
                            multiple_names[node[1]].append(orig)
                        multiple_names[node[1]].append(node[0])
                        tax_names[node[1]] = node[0]
                        tax_names_reverse[node[0]] = node[1]
                    else:
                        tax_names[node[1]] = node[0]
                        tax_names_reverse[node[0]] = node[1]
            if "scientific" in line:
                tax_names_sci[node[0]] = node[1]
    return tax_names_reverse, tax_names, tax_names_sci, multiple_names


"""
def readNodes(nodes_tax_file):

    input:
    - nodes.dmp (NCBI Taxonomy)
    output:
    - dictionary of form {parent: node}
    - dictionary of form {node: type}

    tax_nodes = {}
    tax_types = {}
    with open(nodes_tax_file, 'r') as nodes_tax:
        for line in nodes_tax:
            node = [field.strip() for field in line.split('|')]         #make list of line
            if node[1] not in tax_nodes:
                tax_nodes[node[1]]=[]
            tax_nodes[node[1]].append(node[0])                                #couple node with parent
            tax_types[node[0]] = node[2]                                #couple node with rank
    return tax_nodes,tax_types
"""


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
    return tax_nodes, tax_types


def getTaxChildren(tax_nodes, tax_types, taxid, ranking):
    """
    input:
    - dictionary of form {parent: node} (readNodes output)
    - dictionary of form {node: type} (readNodes output)
    - taxid
    output:
    - dictionary of form {tax_id: [descendants]}
    """

    tax_descendants = []
    children = []
    node = str(taxid)
    if node not in tax_types:  # check if node in nodes.dmp
        print(
            "[Warning] Could not find {} in nodes.dmp while parsing taxonomy hierarchy\n".format(
                node
            )
        )

    else:
        children.append(node)
        for i in children:
            if tax_types[i] == ranking:
                tax_descendants.append(i)
            else:
                if i in tax_nodes:
                    for j in tax_nodes[i]:
                        children.append(j)

    return tax_descendants


def getTaxParent(tax_nodes, tax_types, taxid, ranking):
    """
    input:
    - dictionary of form {parent: node} (readNodes output)
    - dictionary of form {node: type} (readNodes output)
    - taxid
    output:
    - dictionary of form {tax_id: [descendants]}
    """
    tax_parents = {}
    node = str(taxid)
    if node not in tax_types:  # check if node in nodes.dmp
        tax_parents[node] = None
        print(
            "[Warning] Could not find {} in nodes.dmp while parsing taxonomy hierarchy\n".format(
                node
            )
        )

    else:
        parent = tax_nodes[node]  # get parent for current node
        tax_parents[node] = [parent]  # add node to dictionary

        while (
            parent != tax_nodes[parent] and tax_types[parent] != ranking
        ):  # stop when parent = parent(parent) (i.e. 1 = 1)
            parent = tax_nodes[parent]  # get parent of parent
            tax_parents[node].append(parent)  # add parent to node in dictionary

    return tax_parents


api_instance = NcbiApi(args.key)

# determine the lineage where your tax id belongs to (lineage taken until upper level = args.type)
taxparents, taxtypes = readNodes(args.nodesfile)
taxnames, namestax, taxnames_sci, multiple_names = readNames(args.namesfile)

eukgens = []
prokgens = []

spoifamily = ""
spoiclade = ""
spoigenus = args.spoi.split()[0]
spoispecies = args.spoi
if spoispecies in namestax:
    lineage = getTaxParent(taxparents, taxtypes, namestax[spoispecies], args.type)
    lineage2 = getTaxParent(taxparents, taxtypes, namestax[spoispecies], "order")
    if lineage[namestax[spoispecies]] != None:
        spoifamily = taxnames[lineage[namestax[spoispecies]][-1]]
        spoiclade = taxnames[lineage2[namestax[spoispecies]][-1]]

print(spoigenus + "\t" + spoifamily + "\t" + spoiclade)

k = open(args.tax, "r")
for line in k:
    line = line.strip()
    sciname = line.split(";")[-2]
    if "environmental" in sciname:
        sciname = line.split(";")[-3]
    if sciname == "uncultured":
        sciname = line.split(";")[-3]
    if sciname == "endosymbionts":
        sciname = line.split(";")[-3]
    if "Hafnia-Obesumbacterium" in sciname:
        sciname = sciname.split("-")[0]
    if "Escherichia-Shigella" in sciname:
        sciname = sciname.split("-")[0]
    if (
        sciname not in namestax
        and sciname != "Unclassified"
        and sciname != "Chloroplast"
        and sciname != "Mitochondrion"
    ):
        sciname = line.split(";")[-3]
        print("NOW: " + sciname)
    print(sciname)
    if sciname in namestax:
        # print(sciname)
        taxid_line = 0
        if sciname in multiple_names:
            print(sciname)
            found_true_lineage = False
            besttaxid = ""
            bestcounter = 0
            for elem in multiple_names[sciname]:
                lineage = getTaxParent(taxparents, taxtypes, elem, args.type)
                fulllineage = getTaxParent(taxparents, taxtypes, elem, "superkingdom")
                counterhere = 0
                for x in fulllineage[elem]:
                    # print(x)
                    cnt = line.split(";").count(taxnames[x])
                    if cnt > 0:
                        print(taxnames[x])
                        counterhere = counterhere + 1
                if int(counterhere) > int(bestcounter):
                    bestcounter = counterhere
                    besttaxid = elem
            taxid_line = besttaxid
            lineage = getTaxParent(taxparents, taxtypes, taxid_line, args.type)
            fulllineage = getTaxParent(taxparents, taxtypes, taxid_line, "superkingdom")
            cladelineage = getTaxParent(taxparents, taxtypes, taxid_line, "order")
            rootlevelname = taxnames[fulllineage[taxid_line][-1]]
            cladelevelname = taxnames[cladelineage[taxid_line][-1]]
        else:
            taxid_line = namestax[sciname]
            lineage = getTaxParent(taxparents, taxtypes, namestax[sciname], args.type)
            fulllineage = getTaxParent(
                taxparents, taxtypes, namestax[sciname], "superkingdom"
            )
            cladelineage = getTaxParent(
                taxparents, taxtypes, namestax[sciname], "order"
            )
            rootlevelname = taxnames[fulllineage[namestax[sciname]][-1]]
            cladelevelname = taxnames[cladelineage[namestax[sciname]][-1]]
        if taxtypes[taxid_line] == args.type:
            print("FAMILY:" + sciname + " CLADE:" + cladelevelname)
            taxlevelname = sciname
            if "Eukaryota" == rootlevelname:
                n_genomes = api_instance.assembly_count_for_taxon(
                    taxon=str(taxlevelname)
                )
            else:
                n_genomes = api_instance.assembly_count_for_taxon(
                    taxon=str(taxlevelname),
                    filters_assembly_source="refseq",
                )
            foundlevel = False
            if n_genomes > 0:
                foundlevel = True
            if foundlevel:
                if (
                    "Eukaryota" == rootlevelname
                    and taxlevelname not in eukgens
                    and taxlevelname != spoifamily
                    and cladelevelname != spoiclade
                ):
                    fulllineage_euk = taxlevelname
                    for elem in getTaxParent(
                        taxparents, taxtypes, namestax[taxlevelname], "superkingdom"
                    )[namestax[taxlevelname]]:
                        fulllineage_euk = fulllineage_euk + "," + taxnames[elem]
                    print(fulllineage_euk)
                    eukgens.append(fulllineage_euk)
                elif (
                    taxlevelname not in prokgens
                    and taxlevelname != spoifamily
                    and cladelevelname != spoiclade
                ):
                    prokgens.append(taxlevelname)
            else:
                print("No genomes in databases")
                if cladelevelname != "root":
                    taxlevelname = cladelevelname
                    if (
                        "Eukaryota" == rootlevelname
                        and taxlevelname not in eukgens
                        and taxlevelname != spoifamily
                        and cladelevelname != spoiclade
                    ):
                        fulllineage_euk = taxlevelname
                        for elem in getTaxParent(
                            taxparents, taxtypes, namestax[taxlevelname], "superkingdom"
                        )[namestax[taxlevelname]]:
                            fulllineage_euk = fulllineage_euk + "," + taxnames[elem]
                        print(fulllineage_euk)
                        n_genomes = api_instance.assembly_count_for_taxon(
                            taxon=str(taxlevelname.split(",")[0])
                        )
                        if n_genomes > 0:
                            eukgens.append(fulllineage_euk)
                    elif (
                        taxlevelname not in prokgens
                        and taxlevelname != spoifamily
                        and cladelevelname != spoiclade
                    ):
                        n_genomes = api_instance.assembly_count_for_taxon(
                            taxon=str(taxlevelname)
                        )
                        if n_genomes > 0:
                            prokgens.append(taxlevelname)
        elif lineage[taxid_line] != None:
            print("DIFFERENT THAN FAMILY:" + sciname + " CLADE:" + cladelevelname)
            taxlevelname = taxnames_sci[lineage[taxid_line][-1]]
            if int(lineage[taxid_line][-1]) != 1:
                if "Eukaryota" == rootlevelname:
                    n_genomes = api_instance.assembly_count_for_taxon(
                        taxon=str(taxlevelname)
                    )
                else:
                    n_genomes = api_instance.assembly_count_for_taxon(
                        taxon=str(taxlevelname),
                        filters_assembly_source="refseq",
                    )
                foundlevel = False
                if n_genomes > 0:
                    foundlevel = True
                if foundlevel == True:
                    # if taxtypes[namestax[sciname]] == args.type:
                    if (
                        "Eukaryota" == rootlevelname
                        and taxlevelname not in eukgens
                        and taxlevelname != spoifamily
                        and cladelevelname != spoiclade
                    ):
                        fulllineage_euk = taxlevelname
                        for elem in getTaxParent(
                            taxparents, taxtypes, namestax[taxlevelname], "superkingdom"
                        )[namestax[taxlevelname]]:
                            fulllineage_euk = fulllineage_euk + "," + taxnames[elem]
                        print(fulllineage_euk)
                        eukgens.append(fulllineage_euk)
                    elif (
                        taxlevelname not in prokgens
                        and taxlevelname != spoifamily
                        and cladelevelname != spoiclade
                    ):
                        prokgens.append(taxlevelname)
                else:
                    print("No genomes in databases")
                    if cladelevelname != "root":
                        taxlevelname = cladelevelname
                        if (
                            "Eukaryota" == rootlevelname
                            and taxlevelname not in eukgens
                            and taxlevelname != spoifamily
                            and cladelevelname != spoiclade
                        ):
                            fulllineage_euk = taxlevelname
                            for elem in getTaxParent(
                                taxparents,
                                taxtypes,
                                namestax[taxlevelname],
                                "superkingdom",
                            )[namestax[taxlevelname]]:
                                fulllineage_euk = fulllineage_euk + "," + taxnames[elem]
                            print(fulllineage_euk)
                            n_genomes = api_instance.assembly_count_for_taxon(
                                taxon=str(taxlevelname.split(",")[0])
                            )
                            if n_genomes > 0:
                                eukgens.append(fulllineage_euk)
                        elif (
                            taxlevelname not in prokgens
                            and taxlevelname != spoifamily
                            and cladelevelname != spoiclade
                        ):
                            genome_summary = api_instance.assembly_count_for_taxon(
                                taxon=str(taxlevelname)
                            )
                            i = 0
                            if n_genomes > 0:
                                prokgens.append(taxlevelname)
    else:
        print("NOT FOUND:" + sciname + "\t" + line)

file1 = args.outdir + "/prok." + args.suffix
k = open(file1, "w")
for elem in prokgens:
    k.write(elem + "\n")
k.close()

file1 = args.outdir + "/euk." + args.suffix
k = open(file1, "w")
for elem in eukgens:
    k.write(elem + "\n")
k.close()
