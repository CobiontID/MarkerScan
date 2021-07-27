from __future__ import division
import json
import time
import argparse
import os
try:
    import ncbi.datasets
except ImportError:
    print('ncbi.datasets module not found. To install, run `pip install ncbi-datasets-pylib`.')
from ncbi.datasets.package import dataset

parser = argparse.ArgumentParser()
parser.add_argument('--taxname', action="store", dest="tax", type=str, help='scientific name of species of interest')
parser.add_argument('--dir', action="store", dest="dir", type=str, help='base directory')
parser.add_argument("-na", type=str, action='store', dest='namesfile', metavar='NAMES',help='NCBI names.dmp')
parser.add_argument("-no", type=str, action='store', dest='nodesfile', metavar='NODES',help='NCBI nodes.dmp')
parser.add_argument('--refseq', action="store", dest="refs", type=str, help='all or refseq database')
args = parser.parse_args()

def Average(lst): 
    return sum(lst) / len(lst) 

def readNames(names_tax_file):
    '''
    input:
    - name.dmp (NCBI Taxonomy)
    output:
    - dictionary of form {node: name}
    - dictionary of form {sci name: node}
    '''
    tax_names = {}
    tax_names_reverse= {}
    with open(names_tax_file, 'r') as nodes_tax:
        for line in nodes_tax:
            node = [field.strip() for field in line.split('|')]
            if 'scientific' in line or 'synonym' in line:
                tax_names[node[1]] = node[0]
                tax_names_reverse[node[0]] = node[1]
    return tax_names_reverse,tax_names

def readNodes(nodes_tax_file):

    '''
    input:
    - nodes.dmp (NCBI Taxonomy)
    output:
    - dictionary of form {parent: node}
    - dictionary of form {node: type}
    '''

    tax_nodes = {}
    tax_types = {}
    with open(nodes_tax_file, 'r') as nodes_tax:
        for line in nodes_tax:
            node = [field.strip() for field in line.split('|')]         #make list of line
            tax_nodes[node[0]] = node[1]                                #couple node with parent
            tax_types[node[0]] = node[2]                                #couple node with rank
    return tax_nodes


api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())

if not os.path.exists(args.dir):
    os.makedirs(args.dir)

taxparents=readNodes(args.nodesfile)
taxnames,namestax=readNames(args.namesfile)

if args.tax in namestax:
    taxid = namestax[args.tax]
    parent = taxparents[taxid]
    parentname = taxnames[parent]
    foundlevel = False
    print(str(taxid)+'\t'+str(parent)+'\t'+parentname)
    while parent != taxparents[parent]:
        time.sleep(1)
        parentname = taxnames[parent]
        parentname_combi = parentname
        if ' ' in parentname:
            parentname_combi = parentname.replace(' ','_')
        print(str(parent)+'\t'+parentname+'\t'+parentname_combi+'\t'+args.tax)
        # fetch data from NCBI via 'datasets' of all species from that clade
        if args.refs == 'yes':
            genome_summary = api_instance.assembly_descriptors_by_taxon(taxon=str(parent),filters_assembly_source='refseq')
        else:
            genome_summary = api_instance.assembly_descriptors_by_taxon(taxon=str(parent))
        i=0
        for assembly in map(lambda d: d.assembly, genome_summary.assemblies):
            #print(assembly.org.sci_name)
            sciname_orig=assembly.org.sci_name
            strainname=assembly.org.strain
            #print(sciname_orig)
            if strainname != None:
                if strainname in sciname_orig and 'sp' not in sciname_orig:
                    sciname=sciname_orig.replace(strainname,'').strip()
                else:
                    sciname=sciname_orig
            else:
                sciname=sciname_orig
            #print(sciname)
            if sciname != args.tax:
                i=i+1
        if i > 0:
            break
        parent = taxparents[parent]
    
    SpeciesDictionary={}
    i=0
    for assembly in map(lambda d: d.assembly, genome_summary.assemblies):
        contiguity=int(assembly.contig_n50)
        i=i+1
        acc=assembly.assembly_accession
        sciname_orig=assembly.org.sci_name
        strainname=assembly.org.strain
        if strainname != None:
            if strainname in sciname_orig and 'sp' not in sciname_orig:
                sciname=sciname_orig.replace(strainname,'').strip()
            else:
                sciname=sciname_orig
        else:
            sciname=sciname_orig
        if sciname != args.tax:
            print('SECOND STEP')
            #print(assembly.org.sci_name)
            #print(sciname_orig)
            print(sciname) 
            if sciname not in SpeciesDictionary:
                #print(sciname)
                SpeciesDictionary[sciname]={}
                SpeciesDictionary[sciname]['ReleaseDate']=assembly.submission_date
                SpeciesDictionary[sciname]['Identifier']=acc
                SpeciesDictionary[sciname]['GenomeSize']=int(assembly.seq_length)
                SpeciesDictionary[sciname]['N50']=contiguity
            else:
                novel_submission = time.strptime(assembly.submission_date, "%Y-%m-%d")
                old_submission = time.strptime(SpeciesDictionary[sciname]['ReleaseDate'], "%Y-%m-%d")
                if novel_submission > old_submission or contiguity > SpeciesDictionary[sciname]['N50']:
                    SpeciesDictionary[sciname]['ReleaseDate']=assembly.submission_date
                    SpeciesDictionary[sciname]['Identifier']=acc
                    SpeciesDictionary[sciname]['GenomeSize']=int(assembly.seq_length)
                    SpeciesDictionary[sciname]['N50']=contiguity
    
    accs=[]
    for species in SpeciesDictionary:
        print(accs)
        #if not os.path.exists(args.dir2+"/"+SpeciesDictionary[species]['Identifier']):
        accs.append(SpeciesDictionary[species]['Identifier'])
    print(f'Download a package for {accs}.')
    print('Begin download of genome data package ...')
    zipfile_name = str(args.dir)+"/"+"/RefSeq.relatives.zip"
    api_response = api_instance.download_assembly_package(accs,exclude_sequence=False, hydrated='FULLY_HYDRATED',_preload_content=False,filename=zipfile_name)
    with open(zipfile_name, 'wb') as f:
        f.write(api_response.data)
    print('Download complete')
    package = dataset.AssemblyDataset(zipfile_name)
    print(package.get_catalog())
    cmd="unzip -d "+str(args.dir)+"/relatives.Refseq "+str(args.dir)+"/RefSeq.relatives.zip"
    os.system(cmd)
    cmd="cat "+str(args.dir)+"/relatives.Refseq/ncbi_dataset/data/*/*fna > "+str(args.dir)+"/relatives.Refseq/acc.fasta"
    os.system(cmd)
    cmd="dustmasker -in "+str(args.dir)+"/relatives.Refseq/acc.fasta -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' > "+str(args.dir)+"/relatives.Refseq/masked.fna"
    os.system(cmd)
    cmd="rm "+str(args.dir)+"/RefSeq.relatives.zip "+str(args.dir)+"/relatives.Refseq/acc.fasta "
    os.system(cmd)
else:
    print('No taxonomic name was not found')