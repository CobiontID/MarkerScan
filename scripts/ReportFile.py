from __future__ import division
from fpdf import FPDF
import argparse
import configparser
import os
import sys
import glob

parser = argparse.ArgumentParser()
parser.add_argument("-o", type=str, action='store', dest='out',help='define report file')
parser.add_argument("-r", type=str, action='store', dest='rep', help='define removed reads list')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

pdf = FPDF()
pdf.add_page()

wd=os.path.dirname(args.rep)
shortname=wd.split('/')[-1]
totallines=0
ctgs=[]
for filename in glob.glob(wd+'/*readslist'):
    num_lines = sum(1 for line in open(filename))
    totallines=totallines+num_lines
    ctgs = tuple(open(filename, 'r'))

pdf.set_font("Arial", "B", size=12)
#o.write(shortname+'\n')
pdf.cell(200,12,txt=shortname, ln=1, align="C")
pdf.set_font("Arial", size=10)
#o.write("There are "+str(totallines)+" contigs detected containing the SSU locus:\n")
pdf.cell(200, 6, txt="There are "+str(totallines)+" contigs detected containing the SSU locus:", ln=1, align="L")
ctgstring=""
for ctg in ctgs:
    ctg=ctg.strip()
    ctgstring=ctgstring+','+ctg
    #o.write(ctg+'\n')
ctgstringfinal=ctgstring[1:]
pdf.cell(200, 6, txt=ctgstringfinal, ln=1, align="L")

#o.write("\nThese loci are annotated as:\n")
pdf.cell(200, 6,ln=1, align="L")
pdf.cell(200, 6, txt="These loci are annotated as:", ln=1, align="L")
annotated=wd+'/'+shortname+'.ProkSSU.reduced.SILVA.genus.txt'
k=open(annotated,'r')
for line in k:
    line=line.strip()
    #o.write(line+'\n')
    pdf.cell(200, 6, txt=line, ln=1, align="L")
k.close()

for filename in glob.glob(wd+'/*/kraken.reads'):
    genusname=filename.split('/')[-2]
    pdf.set_font("Arial", "B", size=12)
    pdf.cell(200,12,txt=genusname, ln=1, align="L")
    pdf.set_font("Arial", size=10)
    num_lines = sum(1 for line in open(filename))
    percentage=0
    krakenrep=wd+'/kraken.report'
    m=open(krakenrep,'r')
    for line in m:
        line=line.strip()
        if genusname in line:
            percentage=float(line.split('\t')[0])
    m.close()
    pdf.cell(200, 6, txt="There are "+str(num_lines)+" reads ("+str(percentage)+"%) classified by Kraken as "+genusname+"." , ln=1, align="L")
    pdf.cell(200, 6,ln=1, align="L")
    completeness=wd+'/'+genusname+'/busco/completeness_per_contig.txt'
    num_contigs = sum(1 for line in open(completeness)) - 1
    totallen=0
    l=open(completeness,'r')
    for line in l:
        line=line.strip()
        if not line.startswith('#'):
            totallen=totallen+int(line.split('\t')[5])
        line=line.replace('\t',';')
        pdf.cell(200, 6, txt=line, ln=1, align="L")
    l.close()
    mblen="{:.2f}".format(float(totallen/1000000))+"Mb"
    for buscooutput in glob.glob(wd+'/'+genusname+'/busco/busco/short_summary.specific.*.busco.txt'):
        k=open(buscooutput,'r')
        for line in k:
            line=line.strip()
            if line.startswith('C'):
                pdf.cell(200, 6, txt=line, ln=1, align="L")
        k.close()
    readids=wd+'/'+genusname+'/'+genusname+'.readsids.txt'
    num_lines = sum(1 for line in open(readids))
    pdf.cell(200, 6,ln=1, align="L")
    pdf.cell(200, 6, txt="There are "+str(num_lines)+" reads mapping to the full length of "+str(num_contigs)+" contigs ("+mblen+") containing BUSCO genes." , ln=1, align="L")
pdf.output(args.out)
