import argparse
import glob
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, action='store', dest='asm', metavar='asm',help='define asm')
parser.add_argument("-c", type=str, action='store', dest='contigs', metavar='contigs',help='define contig list')
parser.add_argument("-b", type=str, action='store', dest='buscodir', metavar='buscodir',help='define busco dir')
parser.add_argument("-k", type=str, action='store', dest='karyo', metavar='karyo',help='define karyotype')
parser.add_argument("-d", type=str, action='store', dest='dat', metavar='dat',help='define data file')
parser.add_argument("-l", type=str, action='store', dest='link', metavar='link',help='define link file')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

contiglist=[]
f =open(args.contigs,'r')
for record in f:
    record=record.strip()
    contiglist.append(record)
f.close()

m=open(args.karyo,'w')
faidxfile = args.asm+'.fai'
f =open(faidxfile, 'r')
i=0
for record in f:
    record=record.strip()
    contigname=record.split('\t')[0]
    if contigname in contiglist:
        i=i+1
        m.write('chr - '+contigname+' '+str(i)+'  0 '+record.split('\t')[1]+'\t'+contigname+'\n')

f.close()
m.close()

duplications={}
l=open(args.dat,'w')
m=open(args.link,'w')
dirname=args.buscodir.split('/done.txt')[0]+'/busco'
for filename in glob.glob(dirname+'/run_*/full_table.tsv'):
    k=open(filename,'r')
    for line in k:
        if not line.startswith('#'):
            buscogene=line.split()[0]
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
                l.write(contig+'	'+line.split('\t')[3]+'	'+line.split('\t')[4]+'	1	color=red'+'\n')
                if 'Duplicated' in line:
                    if buscogene not in duplications:
                        duplications[buscogene]=[]
                    contigfullname=contig+';'+line.split('\t')[3]+';'+line.split('\t')[4]
                    duplications[buscogene].append(contigfullname)
l.close()

for busco in duplications:
    combinations=list(itertools.combinations(duplications[busco],2))
    i=0
    for pair in combinations:
        busconame=busco+'.'+str(i)
        m.write(busconame+'	'+combinations[i][0].split(';')[0]+'	'+combinations[i][0].split(';')[1]+'	'+combinations[i][0].split(';')[2]+'\n')   
        m.write(busconame+'	'+combinations[i][1].split(';')[0]+'	'+combinations[i][1].split(';')[1]+'	'+combinations[i][1].split(';')[2]+'\n')
        i=i+1
m.close()             