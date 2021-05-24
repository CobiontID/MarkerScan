from __future__ import division
import argparse
import configparser
import os
import glob
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-r", type=str, action='store', dest='readinf',help='define info on reads mapped')
parser.add_argument("-c", type=str, action='store', dest='contigs',help='define contigs')
parser.add_argument("-o", type=str, action='store', dest='out', help='define outputfile')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

finalcontigs=[]
o=open(args.contigs,'r')
for line in o:
    line=line.strip()
    finalcontigs.append(line)
o.close()

readlist=[]
m=open(args.readinf,'r')
for line in m:
    line=line.strip()
    ctg=line.split('\t')[0]
    if ctg in finalcontigs:
        readnames=line.split('\t')[1].split(',')
        for read in readnames:
           #if read not in readlist:
           readlist.append(read)
m.close()

n=open(args.out,'w')
for rd in set(readlist):
    n.write(rd+'\n')
n.close()
