from __future__ import division
import argparse
import configparser
import os
import sys
import glob
import json

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str, action='store', dest='directory', metavar='DIR',help='define refseq directory')
parser.add_argument("-o", type=str, action='store', dest='output', metavar='OUT',help='define output file')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

taxfile=results.directory +'/ncbi_dataset/data/assembly_data_report.jsonl'

with open(taxfile, 'r') as json_file:
    json_list = list(json_file)

SpeciesDictionary={}
for json_str in json_list:
    distro = json.loads(json_str)
    if 'refseqAssmAccession' in distro['assemblyInfo']:
        acc=distro['assemblyInfo']['refseqAssmAccession']
        #if acc == 'na':
        acc2=distro['assemblyInfo']['genbankAssmAccession']
    else:
        acc=distro['assemblyInfo']['genbankAssmAccession']
    taxid=distro['taxId']
    SpeciesDictionary[acc]=taxid
    SpeciesDictionary[acc2]=taxid
    print(acc+'\t'+str(taxid))        

g=open(results.output,'w')
files = glob.glob(results.directory + '/*/*/*/*.fna', recursive=True)
for fastafile in files:
    print(fastafile)
    #if os.path.isfile(os.path.join(os.path.dirname(fastafile),"data_report.yaml")):
    #    taxfile=os.path.join(os.path.dirname(fastafile),"data_report.yaml")
    #    print(fastafile+'\t'+taxfile)
    taxID=0    
    #   l=open(taxfile,'r')
    #   for record in l:
    #        record=record.rstrip()
    #        if 'tax_id' in record:
    #            taxID=record.split(':')[1].strip()
    #    l.close()
    f =open(fastafile,'r')
    for record in f:
        if record.startswith(">"):
            record=record.strip()
            accession=fastafile.split('/')[-2]
            #print(accession)
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

#f =open(results.fasta,'r')
#for record in f:
#        if record.startswith(">"):
#            record=record.strip()
#            newline=record.split(' ')[0]+"|kraken:taxid|"+str(results.tax).strip()+" "+" ".join(record.split(' ')[1:])
#            print(newline)
#        else:
#            record=record.strip()
#            print(record)
