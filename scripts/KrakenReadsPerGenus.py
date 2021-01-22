from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, action='store', dest='input', metavar='INPUT',help='define kraken out file')
parser.add_argument("-rep", type=str, action='store', dest='report', metavar='REPORT',help='define kraken report file')
parser.add_argument("-g", type=str, action='store', dest='genus', metavar='GENUS',help='define prokaryotic genus of interest')
parser.add_argument("-r", type=str, action='store', dest='reads', metavar='READS',help='define reads list output file')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

taxname=str(results.genus).split("genus.")[1].split('.')[0]
if '_' in taxname:
    taxname=taxname.replace('_',' ')

f =open(results.report,'r')
UNDER=False
ABOVE=True
taxlist=[]
for record in f:
    record=record.rstrip()
    elems=record.split('\t')
    if elems[5].strip() == taxname: #and elems[3] == 'G':
        UNDER=True
        spacelength=int(len(elems[5]) - len(elems[5].lstrip(' ')))
        #print("spacelen:"+str(spacelength)+'\t'+str(len(elems[5].lstrip(' ')))+','+elems[5])
        taxlist.append(elems[4])
    elif UNDER == True and ABOVE == True:
        spacehere=int(len(elems[5]) - len(elems[5].lstrip(' ')))
        if spacehere <= spacelength:
        #if not elems[3].startswith('S') and elems[3] != "G1":
            ABOVE = False
        else:
            taxlist.append(elems[4])
            #print(elems[4])
f.close()

l = open(results.input,'r')
m = open(results.reads,'w')
for record in l:
    record=record.strip()
    readname=record.split('\t')[1]
    if record.split('\t')[2].isdigit():
        taxid=record.split('\t')[2]
    else:
        taxid=record.split('\t')[2].split('taxid ')[1].split(')')[0].strip()
    if str(taxid) in taxlist:
        m.write(readname+'\n')
    #elif record.startswith('U'):
    #    m.write(readname+'\n')
