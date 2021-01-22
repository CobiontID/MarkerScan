from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-na", type=str, action='store', dest='namesfile', metavar='NAMES',help='NCBI names.dmp')
parser.add_argument("-no", type=str, action='store', dest='nodesfile', metavar='NODES',help='NCBI nodes.dmp')
parser.add_argument("-f", type=str, action='store', dest='genome', metavar='GENOME FASTA',help='fasta genome assembly file')
parser.add_argument("-d", type=str, action='store', dest='dir', metavar='WORKDIR',help='define working directory for busco')
parser.add_argument("-db", type=str, action='store', dest='db', help='define available dbs file')
parser.add_argument("-dl", type=str, action='store', dest='download',help='define directory to store busco dbs')
parser.add_argument("-c", type=int, action='store', dest='cpu',help='define cpus')
parser.add_argument("-o", type=str, action='store', dest='out', metavar='OUTFILE',help='define configfile name')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

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
            if 'scientific' in line:
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

taxparents=readNodes(args.nodesfile)
taxnames,namestax=readNames(args.namesfile)

genus=args.out.split('/config')[0].split('/')[-1]
if '_' in genus:
    genus=genus.replace('_',' ')

busco_dbs=[]
m=open(args.db,'r')
for line in m:
    line=line.strip()
    if 'db' in line:
        #if 'eukaryota' in line:
        #    break
        dbname=line.split(' ')[-1].split('_')[0]
        busco_dbs.append(dbname)

buscoset = 'Bacteria'
if genus in namestax:
    taxid = namestax[genus]
    parent = taxparents[taxid]
    while parent != taxparents[parent]:
        if taxnames[parent].lower() in busco_dbs:
            buscoset = taxnames[parent]
            break
        parent = taxparents[parent]

print(genus+'\t'+buscoset)

condadir = os.environ['CONDA_DEFAULT_ENV']

l=open(args.out,'w')
l.write('[busco_run]'+'\n')
l.write('# Input file'+'\n')
l.write('in = '+args.genome+'\n')
l.write('# Run name, used in output files and folder'+'\n')
l.write('out = busco'+'\n')
l.write('# Where to store the output directory'+'\n')
l.write('out_path = '+args.dir+'\n')
l.write('# Path to the BUSCO dataset'+'\n')
l.write('lineage_dataset =  '+buscoset.lower()+'\n')
l.write('# Which mode to run (genome / proteins / transcriptome)'+'\n')
l.write('mode = genome'+'\n')
l.write('# How many threads to use for multithreaded steps'+'\n')
l.write('cpu = '+str(args.cpu)+'\n')
l.write('# Force rewrite if files already exist (True/False)'+'\n')
l.write(';force = False'+'\n')
l.write('# Local destination path for downloaded lineage datasets'+'\n')
l.write('download_path = '+args.download+'\n')
l.write('[tblastn]'+'\n')
l.write('path = '+condadir+'/bin/'+'\n')
l.write('command = tblastn'+'\n')
l.write('[makeblastdb]'+'\n')
l.write('path = '+condadir+'/bin/'+'\n')
l.write('command = makeblastdb'+'\n')
l.write('[augustus]'+'\n')
l.write('path = '+condadir+'/bin/'+'\n')
l.write('command = augustus'+'\n')
l.write('[etraining]'+'\n')
l.write('path = '+condadir+'/bin/'+'\n')
l.write('command = etraining'+'\n')
l.write('[gff2gbSmallDNA.pl]'+'\n')
l.write('path = '+condadir+'/bin/'+'\n')
l.write('command = gff2gbSmallDNA.pl'+'\n')
l.write('[new_species.pl]'+'\n')
l.write('path = '+condadir+'/bin/'+'\n')
l.write('command = new_species.pl'+'\n')
l.write('[optimize_augustus.pl]'+'\n')
l.write('path = '+condadir+'/bin/'+'\n')
l.write('command = optimize_augustus.pl'+'\n')
l.write('[hmmsearch]'+'\n')
l.write('path = '+condadir+'/bin/'+'\n')
l.write('command = hmmsearch'+'\n')
l.write('[sepp]'+'\n')
l.write('path = '+condadir+'/bin/'+'\n')
l.write('command = run_sepp.py'+'\n')
l.write('[prodigal]'+'\n')
l.write('path = '+condadir+'/bin/'+'\n')
l.write('command = prodigal'+'\n')
l.close()
