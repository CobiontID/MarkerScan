from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, action='store', dest='input', metavar='INPUT',help='define hmmer tabular domain output file')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
results = parser.parse_args()

f =open(results.input,'r')
besthit={}
bestdomain={}
for record in f:
        if not record.startswith("#"):
            readname=record.split()[2]
            evalue=record.split()[12]
            score=record.split()[13]
            SSUtype=record.split()[0]
            domaineval=record.split()[12]
            domainstart=record.split()[6]
            domainstop=record.split()[7]
            #print(readname)
            if readname not in besthit:
                besthit[readname]={}
                besthit[readname]['eval']=evalue
                besthit[readname]['type']=SSUtype
                besthit[readname]['score']=score
                bestdomain[readname]={}
                bestdomain[readname][SSUtype]={}
                bestdomain[readname][SSUtype]['eval']=domaineval
                bestdomain[readname][SSUtype]['start']=domainstart
                bestdomain[readname][SSUtype]['stop']=domainstop
            else:
                if float(evalue) < float(besthit[readname]['eval']):
                    besthit[readname]['eval']=evalue
                    besthit[readname]['score']=score
                    besthit[readname]['type']=SSUtype
                    bestdomain[readname][SSUtype]={}
                    bestdomain[readname][SSUtype]['eval']=domaineval
                    bestdomain[readname][SSUtype]['start']=domainstart
                    bestdomain[readname][SSUtype]['stop']=domainstop
                if SSUtype in bestdomain[readname]:
                    if float(domaineval) < float(bestdomain[readname][SSUtype]['eval']):
                        bestdomain[readname][SSUtype]['eval']=domaineval
                        bestdomain[readname][SSUtype]['start']=domainstart
                        bestdomain[readname][SSUtype]['stop']=domainstop
                else:
                    bestdomain[readname][SSUtype]={}
                    bestdomain[readname][SSUtype]['eval']=domaineval
                    bestdomain[readname][SSUtype]['start']=domainstart
                    bestdomain[readname][SSUtype]['stop']=domainstop
                    
for read in besthit:
    if float(besthit[read]['eval']) < 1E-150:
        rnatype=besthit[read]['type']
        domainlength=abs(int(bestdomain[read][rnatype]['stop'])-int(bestdomain[read][rnatype]['start']))
        if float(bestdomain[read][rnatype]['eval']) < 1E-150 or domainlength > 1000:
            print(read+"\t"+rnatype+'\t'+str(besthit[read]['eval'])+'\t'+str(bestdomain[read][rnatype]['eval'])+'\t'+str(bestdomain[read][rnatype]['start'])+'\t'+str(bestdomain[read][rnatype]['stop'])+'\t'+str(domainlength))
