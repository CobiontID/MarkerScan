from __future__ import division
import argparse
import configparser
import os
import sys
import glob
import json

'''
This scripts takes as an input a data directory where the ncbi dataset is downloaded 
and will reformat all underlying fasta files to include their taxid in each '>' header e.g. >contig1|kraken:taxid|9606
'''

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str, action='store', dest='directory', metavar='DIR',help='define refseq directory')
parser.add_argument("-o", type=str, action='store', dest='output', metavar='OUT',help='define output file')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

## loop over assembly_data_report.jsonl file to have all taxids per asmid
taxfile=results.directory +'/ncbi_dataset/data/assembly_data_report.jsonl'

with open(taxfile, 'r') as json_file:
    json_list = list(json_file)

SpeciesDictionary={}
for json_str in json_list:
    distro = json.loads(json_str)
    if 'refseqAssmAccession' in distro['assemblyInfo']:
        acc=distro['assemblyInfo']['refseqAssmAccession']
        #if acc == 'na':
        if 'genbankAssmAccession' in distro['assemblyInfo']:
            acc2=distro['assemblyInfo']['genbankAssmAccession']
        else:
            acc2=acc
    elif 'genbankAssmAccession' in distro['assemblyInfo']:
        acc=distro['assemblyInfo']['genbankAssmAccession']
        acc2=distro['assemblyInfo']['genbankAssmAccession']
    taxid=distro['taxId']
    SpeciesDictionary[acc]=taxid
    SpeciesDictionary[acc2]=taxid
    print(acc+'\t'+str(taxid))        

## loop over all fasta files and adapt their seq headers
g=open(results.output,'w')
files = glob.glob(results.directory + '/*/*/*/*.fna')
for fastafile in files:
    print(fastafile)
    taxID=0    
    f =open(fastafile,'r')
    for record in f:
        if record.startswith(">"):
            record=record.strip()
            accession=fastafile.split('/')[-2]
            if accession in SpeciesDictionary:
                newline=record.split(' ')[0]+"|kraken:taxid|"+str(SpeciesDictionary[accession])+" "+" ".join(record.split(' ')[1:])
            else:
                newline=record.split(' ')[0]+"|kraken:taxid|"+str(taxID).strip()+" "+" ".join(record.split(' ')[1:])
            g.write(newline+"\n")
        else:
            record=record.strip()
            g.write(record+"\n")
    f.close()
g.close()