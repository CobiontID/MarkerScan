from __future__ import division
import json
import time
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--taxname', action="store", dest="tax", type=str, help='a genus taxname to download refseq genomes for')
parser.add_argument('--dir', action="store", dest="dir", type=str, help='base directory')
parser.add_argument('--refseq', action="store", dest="refs", type=str, help='all or refseq database')
parser.add_argument('-d', action="store", dest="datasets", type=str, help='datasets')
args = parser.parse_args()

def Average(lst): 
    return sum(lst) / len(lst) 

taxname=str(args.tax).split("genus.")[1].split('.')[0]
taxname_orig=taxname
if '_' in taxname:
    taxname=taxname.replace('_',' ')

if not os.path.exists(args.dir):
    os.makedirs(args.dir)

# fetch data from NCBI via 'datasets' of all species from that clade
if args.refs == 'yes':
    cmd=str(args.datasets)+" assembly-descriptors --refseq tax-name '"+str(taxname)+"' > "+str(args.dir)+"/log."+str(taxname_orig)+".json"
else:
    cmd=str(args.datasets)+" assembly-descriptors tax-name '"+str(taxname)+"' > "+str(args.dir)+"/log."+str(taxname_orig)+".json"
print(cmd)
os.system(cmd)

# parse resulting JSON file
with open(str(args.dir)+"/log."+str(taxname_orig)+".json", 'r') as f:
    distros_dict = json.load(f)

if distros_dict:
    SpeciesDictionary={}
    i=0
    j=0 
    for distro in distros_dict["datasets"]:
        i=i+1
        contiguity=int(distro['contig_n50'])
        acc=distro['assembly_accession']
        sciname_orig=distro['org']['sci_name']
        taxid=distro['org']['key']
        if 'strain' in distro['org']:
            strainname=distro['org']['strain']
            if strainname in sciname_orig and 'sp' not in sciname_orig:
                sciname=sciname_orig.replace(strainname,'').strip()
            else:
                sciname=sciname_orig
        else:
            sciname=sciname_orig
        if sciname not in SpeciesDictionary:
            SpeciesDictionary[sciname]={}
            SpeciesDictionary[sciname]['ReleaseDate']=distro['submission_date']
            SpeciesDictionary[sciname]['Identifier']=acc
            SpeciesDictionary[sciname]['GenomeSize']=int(distro['seq_length'])
            SpeciesDictionary[sciname]['N50']=contiguity
        else:
            novel_submission = time.strptime(distro['submission_date'], "%Y-%m-%d")
            old_submission = time.strptime(SpeciesDictionary[sciname]['ReleaseDate'], "%Y-%m-%d")
            if novel_submission > old_submission or contiguity > SpeciesDictionary[sciname]['N50']:
                SpeciesDictionary[sciname]['ReleaseDate']=distro['submission_date']
                SpeciesDictionary[sciname]['Identifier']=acc
                SpeciesDictionary[sciname]['GenomeSize']=int(distro['seq_length'])
                SpeciesDictionary[sciname]['N50']=contiguity

    assemblytxt=""
    genomesizes=[]
    for species in SpeciesDictionary:
        j=j+1
        print(species+"\t"+str(SpeciesDictionary[species]['GenomeSize'])+'\t'+SpeciesDictionary[species]['Identifier']+"\t"+str(SpeciesDictionary[species]['N50']))
        assemblytxt=assemblytxt+","+SpeciesDictionary[species]['Identifier']
        genomesizes.append(SpeciesDictionary[species]['GenomeSize'])

    assemblytxt=assemblytxt[1:]
    nrgfs=int(round(float(Average(genomesizes)/1000000)*20))
    print('Genomes for '+str(j)+' species')
    #print(assemblytxt)
    print("Number of necessary GFs "+str(nrgfs))
    cmd=str(args.datasets)+" download assembly "+assemblytxt+" --filename "+str(args.dir)+"/RefSeq."+str(taxname_orig)+".zip > "+str(args.dir)+"/"+str(taxname_orig)+".download.log"
    os.system(cmd)

else:
    print('No genomes available')
    cmd="touch "+str(args.dir)+"/"+str(taxname_orig)+".download.log"
    os.system(cmd)
