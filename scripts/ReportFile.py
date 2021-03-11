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
    for buscooutput in glob.glob(wd+'/'+genusname+'/buscoReads/busco/short_summary.specific.*.busco.txt'):
        k=open(buscooutput,'r')
        for line in k:
            line=line.strip()
            if line.startswith('C'):
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
    for buscooutput in glob.glob(wd+'/'+genusname+'/busco/busco/short_summary.specific.*.busco.txt'):
        k=open(buscooutput,'r')
        for line in k:
            line=line.strip()
            if line.startswith('C'):
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

    readids=wd+'/'+genusname+'/'+genusname+'.readsids.txt'
    num_lines_reads = sum(1 for line in open(readids))
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

    totalfraction="{:.2f}".format(float(num_lines_reads/num_lines)*100)
    if (float(num_lines_reads/num_lines)*100) < 80:
        pdf.set_font("Arial", "B", size=10)
    pdf.cell(200, 6, txt="A total fraction of "+str(totalfraction)+"% of the classified kraken reads are removed." , ln=1, align="L")
    pdf.set_font("Arial", size=10)
    pdf.cell(200, 6,ln=1, align="L")

    pdf.set_font("Arial", "U", size=10)
    pdf.cell(200, 6, txt="Hifiasm assembly", ln=1, align="L")
    pdf.set_font("Arial", size=10)
    buscocontigs_asm= wd + '/' + genusname + '/' + genusname + '.Assembly.contigs.txt'
    buscocontigs=[]
    l=open(buscocontigs_asm,'r')
    for line in l:
        line=line.strip()
        buscocontigs.append(line)
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
    for buscooutput in glob.glob(wd+'/'+genusname+'/buscoAssembly/busco/short_summary.specific.*.busco.txt'):
        k=open(buscooutput,'r')
        for line in k:
            line=line.strip()
            if line.startswith('C'):
                pdf.cell(200, 6, txt=line, ln=1, align="L")
        k.close()
    putreadids=wd+'/'+genusname+'/'+genusname+'.assembly.unmapped.reads'
    num_lines_put = sum(1 for line in open(putreadids))
    pdf.cell(200, 6, txt="There are "+str(num_lines_put)+" additional reads which are classified as putative contamination", ln=1, align="L") 
    pdf.cell(200, 6, txt="mapping to the full length of "+str(len(buscocontigs))+" contigs ("+b_mblen+") containing BUSCO genes and/or mapping to refseq genomes.", ln=1, align="L")    
    imagename = wd + '/' + genusname + '/circos.png'
    pdf.image(imagename,w=120,h=120)
pdf.output(args.out)
