from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-m", type=str, action='store', dest='map', metavar='MAP',help='define mapping file kraken reads to genome')
parser.add_argument("-r", type=str, action='store', dest='reads', metavar='READS',help='define all reads mapping per ctg')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

contiglist=[]
f =open(results.map,'r')
for record in f:
    record=record.rstrip()
    if not 'NOT COMPLETE' in record:
        contiglist.append(record.split('\t')[0])
    else:
        pctg=float(record.split('\t')[3].split('%')[0])
        if pctg > 20:
            contiglist.append(record.split('\t')[1])
f.close()

m = open(results.reads,'r')
for record in m:
    record=record.strip()
    ctgname=record.split('\t')[0]
    if ctgname in contiglist:
        for read in record.split('\t')[1].split(','):
            print(read)
m.close()