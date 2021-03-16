from __future__ import division
import argparse
import configparser
import os
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, action='store', dest='input', metavar='INPUT',help='define HMM parsed coordinates file')
parser.add_argument("-f", type=str, action='store', dest='fasta', metavar='INPUT',help='define fasta file')
parser.add_argument("-o", type=str, action='store', dest='out', metavar='INPUT',help='define outfile')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

coords={}
m =open(results.input,'r')
for record in m:
   readname=record.split()[0]
   start=record.split()[4]
   stop=record.split()[5]
   coords[readname]=(start,stop)

with open(results.out,"w") as f:
        for seq_record in SeqIO.parse(results.fasta, "fasta"):
                f.write(">"+str(seq_record.id) + "\n")
                #print(coords[seq_record.id][0]+'\t'+coords[seq_record.id][1])
                if int(coords[seq_record.id][0]) > int(coords[seq_record.id][1]):
                    print(str(seq_record.id))
                    f.write(str(seq_record.seq[(int(coords[seq_record.id][1])-1):(int(coords[seq_record.id][0])-1)].reverse_complement())+ "\n")
                else:
                    f.write(str(seq_record.seq[int(coords[seq_record.id][0]):int(coords[seq_record.id][1])]) + "\n")
