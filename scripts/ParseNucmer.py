from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-n", type=str, action='store', dest='nucmer', help='define nucmer coords file')
parser.add_argument("-o", type=str, action='store', dest='out',help='define contig file')
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
f =open(args.nucmer,'r')
for record in f:
        if record[0].isdigit():
            query=record.split('\t')[11]
            querylen=int(record.split('\t')[7])
            start=int(record.split('\t')[0])
            stop=int(record.split('\t')[1])
            coords=(int(record.split('\t')[0]),int(record.split('\t')[1]))
            if query not in alns:
                alns[query]={}
                alns[query]['length']=querylen
                alns[query]['alns']=[]
            alns[query]['alns'].append(coords)
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
    if percentagectg >= 50:
        finalcontigs.append(ctg)
        #totalreads.extend(alns[ctg]['reads'])
        k.write(ctg+'\t'+str(alns[ctg]['length'])+'\t'+str(percentagectg)+'%\n')
        #print(finalcoords)
    else:
        k.write('NOT COMPLETE:\t'+ctg+'\t'+str(alns[ctg]['length'])+'\t'+str(percentagectg)+'%\n')
        #print(finalcoords)
k.close()