from __future__ import division
import json
import time
import argparse
import os
try:
    import ncbi.datasets
except ImportError:
    print('ncbi.datasets module not found. To install, run `pip install ncbi-datasets-pylib`.')
from ncbi.datasets.package import dataset

parser = argparse.ArgumentParser()
parser.add_argument('--taxname', action="store", dest="tax", type=str, help='a genus taxname to download refseq genomes for')
parser.add_argument('--dir', action="store", dest="dir", type=str, help='base directory')
parser.add_argument('--refseq', action="store", dest="refs", type=str, help='all or refseq database')
args = parser.parse_args()

def Average(lst): 
    return sum(lst) / len(lst) 

api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())

taxname=str(args.tax).split("genus.")[1].split('.')[0]
taxname_orig=taxname
if '_' in taxname:
    taxname=taxname.replace('_',' ')

if not os.path.exists(args.dir):
    os.makedirs(args.dir)

# fetch data from NCBI via 'datasets' of all species from that clade
if args.refs == 'yes':
    genome_summary = api_instance.assembly_descriptors_by_taxon(taxon=str(taxname),filters_reference_only=True,page_size=5000)
else:
    genome_summary = api_instance.assembly_descriptors_by_taxon(taxon=str(taxname),page_size=5000)

SpeciesDictionary={}
i=0
for assembly in map(lambda d: d.assembly, genome_summary.assemblies):
    i=i+1
    contiguity=int(assembly.contig_n50)
    acc=assembly.assembly_accession
    sciname_orig=assembly.org.sci_name
    strainname=assembly.org.strain
    if strainname != None:
        if strainname in sciname_orig and 'sp' not in sciname_orig:
            sciname=sciname_orig.replace(strainname,'').strip()
        else:
            sciname=sciname_orig
    else:
        sciname=sciname_orig
    print(sciname)
    if sciname not in SpeciesDictionary:
        SpeciesDictionary[sciname]={}
        SpeciesDictionary[sciname]['ReleaseDate']=assembly.submission_date
        SpeciesDictionary[sciname]['Identifier']=acc
        SpeciesDictionary[sciname]['GenomeSize']=int(assembly.seq_length)
        SpeciesDictionary[sciname]['N50']=contiguity
    else:
        novel_submission = time.strptime(assembly.submission_date, "%Y-%m-%d")
        old_submission = time.strptime(SpeciesDictionary[sciname]['ReleaseDate'], "%Y-%m-%d")
        if novel_submission > old_submission or contiguity > SpeciesDictionary[sciname]['N50']:
            SpeciesDictionary[sciname]['ReleaseDate']=assembly.submission_date
            SpeciesDictionary[sciname]['Identifier']=acc
            SpeciesDictionary[sciname]['GenomeSize']=int(assembly.seq_length)
            SpeciesDictionary[sciname]['N50']=contiguity


accs=[]
genomesizes=[]
j=0
for species in SpeciesDictionary:
    print(accs)
    j=j+1
    #if not os.path.exists(args.dir2+"/"+SpeciesDictionary[species]['Identifier']):
    accs.append(SpeciesDictionary[species]['Identifier'])
    genomesizes.append(SpeciesDictionary[species]['GenomeSize'])

if len(accs) > 0:
    nrgfs=int(round(float(Average(genomesizes)/1000000)*20))
    print('Genomes for '+str(j)+' species')
    print("Number of necessary GFs "+str(nrgfs))
    print(f'Download a package for {accs}.')
    print('Begin download of genome data package ...')
    print(str(len(accs))+' genomes')
    zipfile_name = str(args.dir)+"/"+"/RefSeq."+str(taxname_orig)+".zip"
    for t in range(0, len(accs), 100):
        accshort=accs[t:t + 100]
        try:
            zipfile_name_part = str(args.dir)+"/"+"/RefSeq."+str(taxname_orig)+".part"+str(t)+".zip"
            api_response = api_instance.download_assembly_package(accshort,exclude_sequence=False, hydrated='FULLY_HYDRATED',_preload_content=False,filename=zipfile_name_part)
            with open(zipfile_name_part, 'wb') as f:
                f.write(api_response.data)
            print('Download complete part '+str(t))
            cmd="unzip -d "+str(args.dir)+"/"+str(taxname_orig)+".RefSeq.part"+str(t)+" "+str(args.dir)+"/"+"/RefSeq."+str(taxname_orig)+".part"+str(t)+".zip"
            os.system(cmd)
        except ncbi.datasets.openapi.ApiException as e:
            print("Exception when calling GenomeApi->download_assembly_package: %s\n" % e)
    cmd="mkdir "+str(args.dir)+"/"+str(taxname_orig)+".Refseq"
    os.system(cmd)
    cmd="mkdir "+str(args.dir)+"/"+str(taxname_orig)+".Refseq/ncbi_dataset"
    os.system(cmd)
    cmd="mkdir "+str(args.dir)+"/"+str(taxname_orig)+".Refseq/ncbi_dataset/data"
    os.system(cmd)
    cmd="cp -r "+str(args.dir)+"/"+str(taxname_orig)+".RefSeq.part*/ncbi_dataset/data/G* "+str(args.dir)+"/"+str(taxname_orig)+".Refseq/ncbi_dataset/data/"
    os.system(cmd)
    cmd="cat "+str(args.dir)+"/"+str(taxname_orig)+".RefSeq.part*/ncbi_dataset/data/assembly_data_report.jsonl >> "+str(args.dir)+"/"+str(taxname_orig)+".Refseq/ncbi_dataset/data/assembly_data_report.jsonl"
    os.system(cmd)
    cmd="cat "+str(args.dir)+"/"+str(taxname_orig)+".RefSeq.part*/ncbi_dataset/data/dataset_catalog.json >> "+str(args.dir)+"/"+str(taxname_orig)+".Refseq/ncbi_dataset/data/dataset_catalog.json"
    os.system(cmd)
    cmd="rm -r "+str(args.dir)+"/"+"/RefSeq."+str(taxname_orig)+".part*.zip "+str(args.dir)+"/"+str(taxname_orig)+".RefSeq.part*"
    os.system(cmd)
else:
    print('No genomes available')
    cmd="touch "+str(args.dir)+"/"+str(taxname_orig)+".download.log"
    os.system(cmd)
