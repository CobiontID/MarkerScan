from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-p", type=str, action='store', dest='paf', metavar='PAF',help='define PAF file')
parser.add_argument("-o", type=str, action='store', dest='out',help='define contig file')
parser.add_argument("-r",type=str, action='store', dest='readfile',help='define contig read mapping file output')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

def sort_condense(ivs):
    if len(ivs) == 0:
        return []
    if len(ivs) == 1:
        if ivs[0][0] > ivs[0][1]:
            return [(ivs[0][1], ivs[0][0])]
        else:
            return ivs
    eps = []
    for iv in ivs:
        ivl = min(iv)
        ivr = max(iv)
        eps.append((ivl, False))
        eps.append((ivr, True))
    eps.sort()
    ret = []
    level = 0
    i = 0
    while i < len(eps)-1:
        if not eps[i][1]:
            level = level+1
            if level == 1:
                left = eps[i][0]
        else:
            if level == 1:
                if not eps[i+1][1] and eps[i+1][0] == eps[i][0]+1:
                    i = i+2
                    continue
                right = eps[i][0]
                ret.append((left, right))
            level = level-1
        i = i+1
    ret.append((left, eps[len(eps)-1][0]))
    return ret

alns={}
f =open(args.paf,'r')
for record in f:
        record=record.strip()
        readname=record.split('\t')[0]
        contig=record.split('\t')[5]
        readlen=int(record.split('\t')[1])
        contiglen=int(record.split('\t')[6])
        start=int(record.split('\t')[2])
        stop=int(record.split('\t')[3])
        coverage=float((stop-start)/readlen)
        coords=(int(record.split('\t')[7]),int(record.split('\t')[8]))
        ending=int(contiglen)-int(record.split('\t')[8])
        endread=int(readlen)-int(stop)
        typemapping=record.split('\t')[12].split(':')[2]
        print(readname+'\t'+str(coverage))
        if coverage >= 0.75 and typemapping == 'P':
            if contig not in alns:
                alns[contig]={}
                alns[contig]['length']=contiglen
                alns[contig]['alns']=[]
                alns[contig]['reads']=[]
            alns[contig]['alns'].append(coords)
            alns[contig]['reads'].append(readname)
        #elif typemapping == 'P':
        #    print(readname)
        elif typemapping == 'P':
            if ending < 20 or int(record.split('\t')[7]) < 20:
                strand=record.split('\t')[4]
                coverage=0
                print(readname)
                #print(str(float(stop/readlen))+'\t'+str(start))
                if int(record.split('\t')[7]) < 20 and strand == '+' and float(stop/readlen) > 0.95:
                    coverage=float((stop-start)/(readlen-start))
                elif ending < 20 and strand == '+' and float(start/readlen) < 0.05:
                    coverage=float((stop-start)/stop)
                elif int(record.split('\t')[7]) < 20 and strand == '-' and float(start/readlen) < 0.05:
                    coverage=float((stop-start)/stop)
                    print(readname)
                elif ending < 20 and strand == '-' and float(stop/readlen) > 0.95:
                    coverage=float((stop-start)/(readlen-start))
                    print(readname)
                print(coverage)
                if coverage >= 0.75:
                    if contig not in alns:
                        alns[contig]={}
                        alns[contig]['length']=contiglen
                        alns[contig]['alns']=[]
                        alns[contig]['reads']=[]
                    alns[contig]['alns'].append(coords)
                    alns[contig]['reads'].append(readname)
f.close()

k=open(args.out,'w')
finalcontigs=[]
for ctg in alns:
    #print(ctg+'\t'+str(alns[ctg]['length']))
    l1=sorted(alns[ctg]['alns'])
    finalcoords=sort_condense(l1)
    totallen=0
    for coord in finalcoords:
        lencoord=coord[1]-coord[0]
        totallen=totallen+lencoord
    percentagectg=(float(totallen/alns[ctg]['length'])*100)
    if percentagectg >= 80:
        finalcontigs.append(ctg)
        #totalreads.extend(alns[ctg]['reads'])
        k.write(ctg+'\t'+str(alns[ctg]['length'])+'\t'+str(percentagectg)+'%\n')
        print(finalcoords)
    else:
        k.write('NOT COMPLETE:\t'+ctg+'\t'+str(alns[ctg]['length'])+'\t'+str(percentagectg)+'%\n')
        print(finalcoords)
k.close()

#s=set(totalreads)
l=open(args.readfile,'w')
#for readname in s:
#    l.write(readname+'\n')
for contig in finalcontigs:
    l.write(contig+'\t'+','.join(alns[contig]['reads'])+'\n')
l.close()
