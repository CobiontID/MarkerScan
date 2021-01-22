from __future__ import division
import json
import time
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--taxname', action="store", dest="tax", type=str, help='scientific name of species of interest')
parser.add_argument('--dir', action="store", dest="dir", type=str, help='base directory')
parser.add_argument("-na", type=str, action='store', dest='namesfile', metavar='NAMES',help='NCBI names.dmp')
parser.add_argument("-no", type=str, action='store', dest='nodesfile', metavar='NODES',help='NCBI nodes.dmp')
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
    while parent != taxparents[parent]:
        parentname = taxnames[parent]
        print(parentname)
        # fetch data from NCBI via 'datasets' of all species from that clade
        if args.refs == 'yes':
            cmd="/nfs/users/nfs_e/ev3/tools/datasets assembly-descriptors --refseq tax-name '"+str(parentname)+"' > "+str(args.dir)+"/log."+str(parentname)+".json"
        else:
            cmd="/nfs/users/nfs_e/ev3/tools/datasets assembly-descriptors tax-name '"+str(parentname)+"' > "+str(args.dir)+"/log."+str(parentname)+".json"
        os.system(cmd)
        # parse resulting JSON file
        with open(str(args.dir)+"/log."+str(parentname)+".json", 'r') as f:
            distros_dict = json.load(f)
        if distros_dict:
            i=0
            for distro in distros_dict["datasets"]:
                i=i+1
            if i > 2:
                foundlevel = True
        f.close()
        if foundlevel == True:
            cmd="mv "+str(args.dir)+"/log."+str(parentname)+".json "+str(args.dir)+"/log."+str(args.tax.replace(' ','_'))+".json"
            os.system(cmd)
            break
        parent = taxparents[parent]

    with open(str(args.dir)+"/log."+str(args.tax.replace(' ','_'))+".json", 'r') as f:
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
        cmd="/nfs/users/nfs_e/ev3/tools/datasets download assembly "+assemblytxt+" --filename "+str(args.dir)+"/RefSeq.relatives.zip > "+str(args.dir)+"/"+str(args.tax.replace(' ','_'))+".download.log"
        os.system(cmd)

else:
    print('No taxonomic name was not found')
    cmd="touch "+str(args.dir)+"/"+str(args.tax).replace(" ","_")+".download.log"
    os.system(cmd)
