from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, action='store', dest='input', metavar='INPUT',help='define list file')
parser.add_argument("-f", type=str, action='store', dest='fasta', metavar='INPUT',help='define fasta file')
parser.add_argument("-o", type=str, action='store', dest='out', metavar='INPUT',help='define outfile')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

ids=[]
m =open(results.input,'r')
for record in m:
    record=record.strip()
    ids.append(record)

f=open(results.out,"w")
l=open(results.fasta,"r")
writebool=False
for record in l:
    record=record.strip()
    if record.startswith('>'):
        if record.split('>')[1] in ids:
            writebool=True
            f.write(record+'\n')
        else:
            writebool=False
    else:
        if writebool == True:
            f.write(record+'\n')