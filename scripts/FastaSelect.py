from __future__ import division
import argparse
import configparser
import os
import sys
import glob
import json

parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, action='store', dest='fasta', metavar='FASTA',help='define fasta file')
parser.add_argument("-l", type=str, action='store', dest='list', metavar='LIST',help='define list')
parser.add_argument("-o", type=str, action='store', dest='output', metavar='OUT',help='define output file')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

accids={}
m=open(results.list,'r')
for record in m:
    record=record.strip()
    accids[record.split('\t')[0]]=record.split('\t')[2]
m.close()

g=open(results.output,'w')
found=False
f =open(results.fasta,'r')
for record in f:
    if record.startswith(">"):
        record=record.strip()
        accession=record.split('>')[1].split(' ')[0]
        if accession in accids:
            found=True
            newline=record.split(' ')[0]+"|kraken:taxid|"+str(accids[accession])+" "+" ".join(record.split(' ')[1:])
            g.write(newline+"\n")
        else:
            found=False
    else:
        if found == True:
            record=record.strip()
            g.write(record+"\n")
f.close()
g.close()
