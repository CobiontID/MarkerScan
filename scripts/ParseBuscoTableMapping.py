from __future__ import division
import argparse
import configparser
import os
import glob
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str, action='store', dest='dir', help='define busco dir')
parser.add_argument("-i", type=str, action='store', dest='assinf', help='define info on contigs')
parser.add_argument("-r", type=str, action='store', dest='readinf',help='define info on reads mapped')
parser.add_argument("-o", type=str, action='store', dest='out', help='define outputfile')
parser.add_argument("-o2", type=str, action='store', dest='out2', help='define outputfile readlist')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

total_genes=[]
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
            if buscogene not in total_genes:
                total_genes.append(buscogene)
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

contig_assembly={}
k=open(args.assinf,'r')
for line in k:
    if not line.startswith('#') and not 'NOT COMPLETE' in line:
        contig=line.split('\t')[0]
        contig_assembly[contig]={}
        contig_assembly[contig]['length']=int(line.split('\t')[1])
        #contig_assembly[contig]['coverage']=int(line.split('\t')[2])
        #contig_assembly[contig]['circ']=line.split('\t')[3]
        #contig_assembly[contig]['multipl']=int(line.split('\t')[5])
k.close()

finalcontigs=[]
l=open(args.out,'w')
l.write('#contig\tfound\ttotal\tcompleteness\tdensity per 100kb\tlen\n')
for ctg in contig_assembly:
    if ctg in contig_info:
        finalcontigs.append(ctg)
        percentage=float(len(set(contig_info[ctg]))/len(total_genes)*100)
        contig_assembly[ctg]['completeness']=percentage
        GenesPerLen=float(len(set(contig_info[ctg]))/contig_assembly[ctg]['length']*100000)
        l.write(ctg+'\t'+str(len(set(contig_info[ctg])))+'\t'+str(len(total_genes))+'\t'+"{:.2f}".format(percentage)+'%\t'+"{:.2f}".format(GenesPerLen)+'\t'+str(contig_assembly[ctg]['length'])+'\n')
    #else:
    #    l.write(ctg+'\t0\t'+str(len(total_genes))+'\t0%\t0\t'+str(contig_assembly[ctg]['length'])+'\n')
l.close()

readlist=[]
m=open(args.readinf,'r')
for line in m:
    line=line.strip()
    ctg=line.split('\t')[0]
    if ctg in finalcontigs:
        readnames=line.split('\t')[1].split(',')
        for read in readnames:
            if read not in readlist:
                readlist.append(read)
m.close()

n=open(args.out2,'w')
for rd in readlist:
    n.write(rd+'\n')
n.close()
