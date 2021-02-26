from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, action='store', dest='input', metavar='INPUT',help='define fasta file')
parser.add_argument("-o", type=str, action='store', dest='output', metavar='OUTPUT',help='define conversion file')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

k=open(results.output,'w')
i=1
f =open(results.input,'r')
for record in f:
    record=record.strip()
    if record.startswith('>'):
        print('>'+str(i))
        k.write(record.split('>')[1]+'\t'+str(i)+'\n')
        i=i+1
    else:
        print(record)