from __future__ import division
import argparse
import configparser
import os
import glob
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str, action='store', dest='dir', help='define busco dir')
parser.add_argument("-c", type=str, action='store', dest='conv', help='define conversion file')
parser.add_argument("-o", type=str, action='store', dest='out', help='define outputfile')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

contig_info={}
busco_contig={}
dirname=args.dir.split('/done.txt')[0]+'/busco'
for filename in glob.glob(dirname+'/run_*/full_table.tsv'):
    k=open(filename,'r')
    for line in k:
        if not line.startswith('#'):
            buscogene=line.split()[0]
            if buscogene not in busco_contig:
                busco_contig[buscogene]=[]
            if not 'Missing' in line:
                if int(line.split('\t')[2].count("_")) == 2:
                    contig='_'.join(line.split('\t')[2].split('_')[:-1])
                elif int(line.split('\t')[2].count("_")) == 1:
                    if line.split('\t')[2].startswith('a_'):
                        contig=line.split('\t')[2]
                    else:
                        contig=line.split('\t')[2].split('_')[0]
                else:
                    contig=line.split('\t')[2]
                if contig not in contig_info:
                    contig_info[contig]=[]
                #if buscogene not in contig_info[contig]:
                contig_info[contig].append(buscogene)
                busco_contig[buscogene].append(contig)
    k.close()

readinfo={}
k=open(args.conv,'r')
for line in k:
    line=line.strip()
    readinfo[line.split('\t')[1]]=line.split('\t')[0]
k.close()

l=open(args.out,'w')
for read in contig_info:
    l.write(readinfo[read]+'\n')