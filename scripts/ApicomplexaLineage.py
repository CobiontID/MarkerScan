from __future__ import division
import argparse
import configparser
import os
import sys
import glob

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

g=open(results.output,'w')
taxnames=readNames(results.namesfile)
out=open(results.output,'w')
files = glob.glob(results.directory + '/*.fasta')
for fastafile in files:
    print(fastafile)
    f =open(fastafile,'r')
    for record in f:
        if record.startswith(">"):
            record=record.strip()
            accession=' '.join(record.split('>')[1].split(' ')[1:3])
            if 'strain' in record:
                accession=' '.join(record.split('>')[1].split(' ')[1:]).split('strain')[0].strip()
            if 'apicoplast' in record:
                accession=' '.join(record.split('>')[1].split(' ')[1:]).split('apicoplast')[0].strip()
            if 'isolate' in record:
                accession=' '.join(record.split('>')[1].split(' ')[1:]).split('isolate')[0].strip()
            if 'strain' in accession:
                accession = accession.split('strain')[0].strip()
            #print(accession)
            if accession in taxnames:
                found=True
                newline=record.split(' ')[0]+"|kraken:taxid|"+str(taxnames[accession])+" "+" ".join(record.split(' ')[1:])
                g.write(newline+"\n")
            else:
                print('NOT FOUND '+accession)
                found=False
        else:
            if found == True:
                record=record.strip()
                g.write(record+"\n")
    f.close()
g.close()
