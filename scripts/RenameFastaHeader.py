from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, action='store', dest='input', metavar='INPUT',help='define fasta file')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

i=1
f =open(results.input,'r')
for record in f:
    record=record.strip()
    if record.startswith('>'):
        print('>'+str(i))
        i=i+1
    else:
        print(record)