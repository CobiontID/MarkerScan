from __future__ import division
import argparse
import configparser
import os
import sys
import glob
import gzip
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str, action='store', dest='directory', metavar='DIR',help='define refseq directory')
parser.add_argument("-na", type=str, action='store', dest='namesfile', metavar='NAMES',help='NCBI names.dmp')
parser.add_argument("-o", type=str, action='store', dest='output', help='output file')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

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
    return tax_names

taxnames=readNames(results.namesfile)
out=open(results.output,'w')
files = glob.glob(results.directory + '/*.gbff.gz')
for gfffile in files:
    print(gfffile)
    with gzip.open(gfffile, "rt") as fh:
        gb = SeqIO.parse(fh, "gb")
        for entry in gb:
            taxon=entry.annotations["organism"]
            lineage= "; ".join(entry.annotations["taxonomy"])
            accession = entry.id
            if taxon in taxnames:
                taxid=taxnames[taxon]
                out.write(accession+'\t'+lineage+'\t'+str(taxid)+'\n')
            else:
                print(taxon)
