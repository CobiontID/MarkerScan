from __future__ import division
from fpdf import FPDF
import argparse
import configparser
import os
import sys
import glob
import json

parser = argparse.ArgumentParser()
parser.add_argument("-o", type=str, action='store', dest='out',help='define report file')
parser.add_argument("-r", type=str, action='store', dest='rep', help='define removed reads list')
parser.add_argument("-d", type=str, action='store', dest='datadir', help='define data dir')
parser.add_argument("-t", type=str, action='store', dest='table', help='define table file')
parser.add_argument("-l", type=str, action='store', dest='list', help='define ssu read list file')
parser.add_argument("-lm", type=str, action='store', dest='listmicro', help='define microsporidia ssu read list file')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

def calculate_N50(list_of_lengths):
    """Calculate N50 for a sequence of numbers.
 
    Args:
        list_of_lengths (list): List of numbers.
 
    Returns:
        float: N50 value.
 
    """
    tmp = []
    for tmp_number in set(list_of_lengths):
            tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()
 
    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]
 
    return median

def Average(lst): 
    return sum(lst) / len(lst) 

pdf = FPDF()
pdf.add_page()
reportdict={}

wd=os.path.dirname(args.rep)
shortname=args.out.split('.report.pdf')[0].split('/')[-1]
totallines=0
ctgstotal=[]
filename= args.list
num_lines = sum(1 for line in open(filename))
totallines=totallines+num_lines
ctgs = list(open(filename, 'r'))
ctgstotal.extend(ctgs)

filenam=args.listmicro
for line in open(filename):
    ctg=line.strip().split(':')[0]
    if ctg not in ctgstotal:
        ctgstotal.append(ctg)
        totallines=totallines+1

pdf.set_font("Arial", "B", size=12)
#o.write(shortname+'\n')
pdf.cell(200,12,txt=shortname, ln=1, align="C")
pdf.set_font("Arial", size=10)
#o.write("There are "+str(totallines)+" contigs detected containing the SSU locus:\n")
pdf.cell(200, 6, txt="There are "+str(totallines)+" contigs detected containing the SSU locus:", ln=1, align="L")
ctgstring=""
for ctg in ctgstotal:
    #print(ctg)
    ctg=ctg.strip()
    ctgstring=ctgstring+','+ctg
    #o.write(ctg+'\n')
ctgstringfinal=ctgstring[1:]
pdf.cell(200, 6, txt=ctgstringfinal, ln=1, align="L")

#reportdict['ContigsWithSSU']=ctgstotal

#o.write("\nThese loci are annotated as:\n")
pdf.cell(200, 6,ln=1, align="L")
pdf.cell(200, 6, txt="These loci are annotated as:", ln=1, align="L")
annotated=wd+'/'+shortname+'.ProkSSU.reduced.SILVA.genus.txt'
#specieslist=[]
k=open(annotated,'r')
for line in k:
    line=line.strip()
    #o.write(line+'\n')
    #specieslist.append(line)
    pdf.cell(200, 6, txt=line, ln=1, align="L")
k.close()
#reportdict['SpeciesPresent']=specieslist

reportdict['SpeciesPresent']={}
annotated2=wd+'/'+shortname+'.ProkSSU.reduced.SILVA.tax'
k=open(annotated2,'r')
for line in k:
    line=line.strip()
    if line.startswith('name'):
        type=line.split('\t')[1].split('lca_tax_')[1]
    else:
        ctgname=line.split('\t')[0]
        tax=line.split('\t')[1]
        if ctgname not in reportdict['SpeciesPresent']:
            reportdict['SpeciesPresent'][ctgname]={}
        reportdict['SpeciesPresent'][ctgname][type]=tax
k.close()

annotated3=wd+'/'+shortname+'.ProkSSU.reduced.fa.clstr'
clusterctgs=[]
k=open(annotated3,'r')
for line in k:
    line=line.strip()
    if line.startswith('>'):
        if len(clusterctgs) > 0:
            reportdict['SpeciesPresent'][representative]['Cluster']=clusterctgs
        clusterctgs=[]
    else:
        ctgname=line.split('\t')[1].split('>')[1].split('...')[0]
        length=line.split('\t')[1].split('nt')[0]
        clusterctgs.append(ctgname)
        if '*' in line:
            representative=ctgname
            reportdict['SpeciesPresent'][ctgname]['SSUlength']=length
reportdict['SpeciesPresent'][representative]['Cluster']=clusterctgs
k.close()
reportdict['Families']={}
for filename in glob.glob(wd+'/*/kraken.reads'):
    genusname=filename.split('/')[-2]
    reportdict['Families'][genusname]={}
    pdf.set_font("Arial", "B", size=12)
    pdf.cell(200,12,txt=genusname, ln=1, align="L")
    pdf.set_font("Arial", size=10)
    if os.path.getsize(filename) > 0:
        num_lines = sum(1 for line in open(filename))
    else:
        num_lines = 0
    percentage=0
    krakenrep=wd+'/kraken.report'
    m=open(krakenrep,'r')
    for line in m:
        #line=line.strip()
        searchpattern=' '+genusname+'\n'
        if searchpattern in line:
            percentage=float(line.split('\t')[0])
    m.close()
    pdf.cell(200, 6, txt="There are "+str(num_lines)+" reads ("+str(percentage)+"%) classified by Kraken as "+genusname+"." , ln=1, align="L")
    reportdict['Families'][genusname]['ClassifiedReads']=num_lines
    reportdict['Families'][genusname]['ClassifiedReadsPercentage']=percentage
    readfile=wd+'/'+genusname+'/'+genusname+".reads2assemble.fa"
    counter=0
    k=open(readfile,'r')
    for line in k:
        line=line.strip()
        if line.startswith('>'):
            counter=counter+1
    k.close()
    if counter > 100000:
        reportdict['Families'][genusname]['Busco_ClassifiedReads']="Too many reads - NA "
    buscooutput =  wd+'/'+genusname+'/buscoReads/summary.txt'
    k=open(buscooutput,'r')
    for line in k:
        line=line.strip()
        if line.startswith('C'):
            reportdict['Families'][genusname]['Busco_ClassifiedReads']=line
            pdf.cell(200, 6, txt=line, ln=1, align="L")
    k.close()
    pdf.cell(200, 6,ln=1, align="L")

    completeness=wd+'/'+genusname+'/busco/completeness_per_contig.txt'
    pdf.set_font("Arial", "U", size=10)
    pdf.cell(200, 6, txt="Busco completeness", ln=1, align="L")
    pdf.set_font("Arial", size=10)
    l=open(completeness,'r')
    for line in l:
        line=line.strip()
        line=line.replace('\t',';')
        pdf.cell(200, 6, txt=line, ln=1, align="L")
    l.close()
    buscooutput =  wd+'/'+genusname+'/busco/summary.txt'
    k=open(buscooutput,'r')
    for line in k:
        line=line.strip()
        if line.startswith('C'):
            reportdict['Families'][genusname]['Busco_Assembly']=line
            pdf.cell(200, 6, txt=line, ln=1, align="L")
    k.close()
    pdf.cell(200, 6,ln=1, align="L")

    pdf.set_font("Arial", "U", size=10)
    pdf.cell(200, 6, txt="Alignment RefSeq genomes", ln=1, align="L")
    pdf.set_font("Arial", size=10)
    nucmerfile=wd+'/'+genusname+'/'+genusname+'_vs_contigs.overview.txt'
    l=open(nucmerfile,'r')
    for line in l:
        line=line.strip()
        if not 'NOT COMPLETE' in line:
            newline=line.split('\t')[0]+';'+"{:.2f}".format(float(line.split('\t')[2].split('%')[0]))+'%'
            pdf.cell(200, 6, txt=newline, ln=1, align="L")
    l.close()
    pdf.cell(200, 6,ln=1, align="L")

    readids=wd+'/'+genusname+'/'+genusname+'.final_reads.fa'
    num_lines_reads=0
    k=open(readids,'r')
    for line in k:
        line=line.strip()
        if line.startswith('>'):
            num_lines_reads=num_lines_reads+1
    k.close()
    fastafile = wd + '/' + genusname + '/' +genusname+'.finalassembly.fa'
    num_contigs=0
    totallen=0
    seqlen=0
    l=open(fastafile,'r')
    for line in l:
        if line.startswith('>'):
            num_contigs=num_contigs+1
            totallen=totallen+seqlen
            seqlen=0
        else:
            seqlen=seqlen+len(line)
    totallen=totallen+seqlen
    mblen="{:.2f}".format(float(totallen/1000000))+"Mb"
    pdf.cell(200, 6, txt="There are "+str(num_lines_reads)+" reads mapping to the full length of "+str(num_contigs)+" contigs ("+mblen+") containing BUSCO genes " , ln=1, align="L")
    pdf.cell(200, 6, txt="and/or mapping to refseq genomes." , ln=1, align="L")
    reportdict['Families'][genusname]['BuscoNucmer_Assembly_Contigs']=num_contigs
    reportdict['Families'][genusname]['BuscoNucmer_Assembly_ContigLength']=mblen
    reportdict['Families'][genusname]['BuscoNucmer_Assembly_Reads']=num_lines_reads

    if int(num_lines) > 0:
        totalfraction="{:.2f}".format(float(num_lines_reads/num_lines)*100)
        if (float(num_lines_reads/num_lines)*100) < 80:
            pdf.set_font("Arial", "B", size=10)
    else:
        totalfraction=0
    pdf.cell(200, 6, txt="A total fraction of "+str(totalfraction)+"% of the classified kraken reads are removed." , ln=1, align="L")
    pdf.set_font("Arial", size=10)
    pdf.cell(200, 6,ln=1, align="L")

    pdf.set_font("Arial", "U", size=10)
    pdf.cell(200, 6, txt="Hifiasm assembly", ln=1, align="L")
    pdf.set_font("Arial", size=10)
    buscocontigs_asm= wd + '/' + genusname + '/' + genusname + '.re-assembly.fa'
    buscocontigs=[]
    l=open(buscocontigs_asm,'r')
    for line in l:
        line=line.strip()
        if line.startswith('>'):
            buscocontigs.append(line.split('>')[1])
    l.close()
    assemblyfile = wd + '/' + genusname + '/hifiasm/hifiasm.p_ctg.fasta.fai'
    print(assemblyfile)
    num_contigs_hifiasm = 0
    if os.path.exists(assemblyfile) and os.path.getsize(assemblyfile) > 0:
        num_contigs_hifiasm = sum(1 for line in open(assemblyfile))
    lengths=[]
    totallen=0
    totallenbusco=0
    if os.path.exists(assemblyfile) and os.path.getsize(assemblyfile) > 0:
        l=open(assemblyfile,'r')
        for line in l:
            idname=line.split('\t')[0]
            length=int(line.split('\t')[1])
            if idname in buscocontigs:
                totallenbusco=totallenbusco+length
            totallen=totallen+length
            lengths.append(length)
        l.close()
    N50=0
    if len(lengths) > 0:
        N50=calculate_N50(lengths)
    kblenN50="{:.2f}".format(float(N50/1000))+"kb"
    mblen="{:.2f}".format(float(totallen/1000000))+"Mb"
    b_mblen="{:.2f}".format(float(totallenbusco/1000000))+"Mb"
    pdf.cell(200, 6, txt="Hifiasm assembled "+str(num_contigs_hifiasm)+" contigs with an N50 of "+kblenN50+" and a total length of "+mblen+".", ln=1, align="L")
    buscooutput =  wd+'/'+genusname+'/buscoAssembly/summary.txt'
    k=open(buscooutput,'r')
    for line in k:
        line=line.strip()
        if line.startswith('C'):
            pdf.cell(200, 6, txt=line, ln=1, align="L")
            reportdict['Families'][genusname]['Busco_Re-Assembly']=line
    k.close()

    putreadids=wd+'/'+genusname+'/'+genusname+'.re-assembly_reads.fa'
    num_lines_put=0
    k=open(putreadids,'r')
    for line in k:
        line=line.strip()
        if line.startswith('>'):
            num_lines_put=num_lines_put+1
    k.close()
    pdf.cell(200, 6, txt="There are "+str(num_lines_put)+" reads mapping to the full length of "+str(len(buscocontigs))+" contigs ("+b_mblen+") containing BUSCO genes " , ln=1, align="L")
    pdf.cell(200, 6, txt="and/or mapping to refseq genomes." , ln=1, align="L")

    reportdict['Families'][genusname]['BuscoNucmer_Re-Assembly_Contigs']=len(buscocontigs)
    reportdict['Families'][genusname]['BuscoNucmer_Re-Assembly_ContigLength']=b_mblen
    reportdict['Families'][genusname]['BuscoNucmer_Re-Assembly_Reads']=num_lines_put
    #refseqfile=args.datadir+'/'+genusname+'/'+genusname+'.refseq.log'
    #genomesize=[]
    #k=open(refseqfile,'r')
    #for line in k:
    #    line=line.strip()
    #    if not line.startswith('/') and not line.startswith('Genomes for') and not line.startswith('Number of'):
    #        genomesize.append(int(line.split('\t')[1]))
    #k.close()
    #average = Average(genomesize)
    #mblen="{:.2f}".format(float(average/1000000))+"Mb"
    #pdf.cell(200, 6, txt="The mean genomesize for the family "+str(genusname)+" is "+str(mblen), ln=1, align="L") 
    imagename = wd + '/' + genusname + '/circos.png'
    if os.path.getsize(imagename) > 0:
        pdf.image(imagename,w=120,h=120)
pdf.output(args.out)

filename_json=args.out.split('.report.pdf')[0]+'.json'
with open(filename_json, 'w') as fp:
    json.dump(reportdict, fp)
