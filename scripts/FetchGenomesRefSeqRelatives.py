from __future__ import division
import json
import time
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--taxname', action="store", dest="tax", type=str, help='scientific name of species of interest')
parser.add_argument('--dir', action="store", dest="dir", type=str, help='base directory')
parser.add_argument('--dir2', action="store", dest="dir2", type=str, help='base data directory')
parser.add_argument("-na", type=str, action='store', dest='namesfile', metavar='NAMES',help='NCBI names.dmp')
parser.add_argument("-no", type=str, action='store', dest='nodesfile', metavar='NODES',help='NCBI nodes.dmp')
parser.add_argument('-o', action="store", dest="out", type=str, help='out directory')
parser.add_argument('-d', action="store", dest="datasets", type=str, help='datasets')
parser.add_argument('--refseq', action="store", dest="refs", type=str, help='all or refseq database')
args = parser.parse_args()

def Average(lst): 
    return sum(lst) / len(lst) 


def readNames(names_tax_file):
    '''
    input:
    - name.dmp (NCBI Taxonomy)
    output:
    - dictionary of form {node: name}
    - dictionary of form {sci name: node}
    '''
    tax_names = {}
    tax_names_reverse= {}
    with open(names_tax_file, 'r') as nodes_tax:
        for line in nodes_tax:
            node = [field.strip() for field in line.split('|')]
            if 'scientific' in line or 'synonym' in line:
                tax_names[node[1]] = node[0]
                tax_names_reverse[node[0]] = node[1]
    return tax_names_reverse,tax_names

def readNodes(nodes_tax_file):

    '''
    input:
    - nodes.dmp (NCBI Taxonomy)
    output:
    - dictionary of form {parent: node}
    - dictionary of form {node: type}
    '''

    tax_nodes = {}
    tax_types = {}
    with open(nodes_tax_file, 'r') as nodes_tax:
        for line in nodes_tax:
            node = [field.strip() for field in line.split('|')]         #make list of line
            tax_nodes[node[0]] = node[1]                                #couple node with parent
            tax_types[node[0]] = node[2]                                #couple node with rank
    return tax_nodes

if not os.path.exists(args.dir):
    os.makedirs(args.dir)

taxparents=readNodes(args.nodesfile)
taxnames,namestax=readNames(args.namesfile)

if args.tax in namestax:
    taxid = namestax[args.tax]
    parent = taxparents[taxid]
    parentname = taxnames[parent]
    foundlevel = False
    print(str(taxid)+'\t'+str(parent)+'\t'+parentname)
    while parent != taxparents[parent]:
        parentname = taxnames[parent]
        parentname_combi = parentname
        if ' ' in parentname:
            parentname_combi = parentname.replace(' ','_')
        print(str(parent)+'\t'+parentname+'\t'+parentname_combi)
        # fetch data from NCBI via 'datasets' of all species from that clade
        if args.refs == 'yes':
            cmd=str(args.datasets)+" assembly-descriptors --refseq tax-id '"+str(parent)+"' > "+str(args.dir)+"/log."+str(parentname_combi)+".json"
        else:
            cmd=str(args.datasets)+" assembly-descriptors tax-id '"+str(parent)+"' > "+str(args.dir)+"/log."+str(parentname_combi)+".json"
        print(cmd)
        os.system(cmd)
        # parse resulting JSON file
        with open(str(args.dir)+"/log."+str(parentname_combi)+".json", 'r') as f:
            distros_dict = json.load(f)
        if distros_dict:
            i=0
            for distro in distros_dict["datasets"]:
                i=i+1
            if i > 2:
                foundlevel = True
        f.close()
        if foundlevel == True:
            cmd="mv "+str(args.dir)+"/log."+str(parentname_combi)+".json "+str(args.dir)+"/log."+str(args.tax.replace(' ','_'))+".json"
            os.system(cmd)
            break
        parent = taxparents[parent]

    with open(str(args.dir)+"/log."+str(args.tax.replace(' ','_'))+".json", 'r') as f:
        distros_dict = json.load(f)
    if distros_dict:
        os.makedirs(args.out+'/ncbi_dataset/data/')
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
            if sciname != args.tax:
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
            accname=SpeciesDictionary[species]['Identifier']
            if not os.path.exists(args.dir2+"/"+accname):
                os.makedirs(args.dir2+"/"+str(accname))
                cmd=str(args.datasets)+" download assembly "+str(accname)+" --filename "+str(args.dir2)+"/"+str(accname)+"/RefSeq.relatives.zip > "+str(args.dir2)+"/"+str(accname)+"/"+str(args.tax.replace(' ','_'))+".download.log"
                os.system(cmd)
                cmd="unzip -d "+str(args.dir2)+"/"+str(accname)+"/relatives.Refseq "+str(args.dir2)+"/"+str(accname)+"/RefSeq.relatives.zip"
                os.system(cmd)
                cmd="cat "+str(args.dir2)+"/"+str(accname)+"/relatives.Refseq/ncbi_dataset/data/"+str(accname)+"/*fna > "+str(args.dir2)+"/"+str(accname)+"/relatives.Refseq/acc.fasta"
                os.system(cmd)
                cmd="dustmasker -in "+str(args.dir2)+"/"+str(accname)+"/relatives.Refseq/acc.fasta -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' > "+str(args.dir2)+"/"+str(accname)+"/relatives.Refseq/masked.fna"
                os.system(cmd)
                cmd="rm -r "+str(args.dir2)+"/"+str(accname)+"/RefSeq.relatives.zip "+str(args.dir2)+"/"+str(accname)+"/relatives.Refseq/acc.fasta "+str(args.dir2)+"/"+str(accname)+"/relatives.Refseq/ncbi_dataset/data/"+str(accname)
                os.system(cmd)
            os.makedirs(args.out+'/ncbi_dataset/data/'+accname)
            cmd="cp "+str(args.dir2)+"/"+str(accname)+"/relatives.Refseq/masked.fna "+args.out+'/ncbi_dataset/data/'+accname+'/masked.fna'
            os.system(cmd)
            cmd="cp "+str(args.dir2)+"/"+str(accname)+"/relatives.Refseq/ncbi_dataset/data/assembly_data_report.jsonl "+args.out+'/ncbi_dataset/data/'+accname+'/assembly_data_report.jsonl'
            os.system(cmd)
            cmd="cat "+args.out+'/ncbi_dataset/data/*/assembly_data_report.jsonl > '+args.out+'/ncbi_dataset/data/assembly_data_report.jsonl'
            os.system(cmd)
        assemblytxt=assemblytxt[1:]
        nrgfs=int(round(float(Average(genomesizes)/1000000)*20))
        print('Genomes for '+str(j)+' species')
        #print(assemblytxt)
        print("Number of necessary GFs "+str(nrgfs))

else:
    print('No taxonomic name was not found')
    cmd="touch "+str(args.dir)+"/"+str(args.tax).replace(" ","_")+".download.log"
    os.system(cmd)
