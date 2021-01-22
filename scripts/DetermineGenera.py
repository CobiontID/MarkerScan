from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, action='store', dest='tax', metavar='TAX',help='define tax genus file')
parser.add_argument("-t", type=str, action='store', dest='type', metavar='TYPE',help='define ranking type')
parser.add_argument("-na", type=str, action='store', dest='namesfile', metavar='NAMES',help='NCBI names.dmp')
parser.add_argument("-no", type=str, action='store', dest='nodesfile', metavar='NODES',help='NCBI nodes.dmp')
parser.add_argument("-suf", type=str, action='store', dest='suffix', metavar='SUFFIX',help='suffix outputfile')
parser.add_argument("-od", type=str, action='store', dest='outdir', metavar='OUTDIR',help='output directory')
parser.add_argument("-g", type=str, action='store', dest='spoi', metavar='SPOI',help='species of interest')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

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

'''
def readNodes(nodes_tax_file):

    input:
    - nodes.dmp (NCBI Taxonomy)
    output:
    - dictionary of form {parent: node}
    - dictionary of form {node: type}

    tax_nodes = {}
    tax_types = {}
    with open(nodes_tax_file, 'r') as nodes_tax:
        for line in nodes_tax:
            node = [field.strip() for field in line.split('|')]         #make list of line
            if node[1] not in tax_nodes:
                tax_nodes[node[1]]=[]
            tax_nodes[node[1]].append(node[0])                                #couple node with parent
            tax_types[node[0]] = node[2]                                #couple node with rank
    return tax_nodes,tax_types
'''

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
    return tax_nodes,tax_types


def getTaxChildren(tax_nodes, tax_types, taxid, ranking):

    '''
    input:
    - dictionary of form {parent: node} (readNodes output)
    - dictionary of form {node: type} (readNodes output)
    - taxid
    output:
    - dictionary of form {tax_id: [descendants]}
    '''

    tax_descendants = []
    children = []
    node = str(taxid)
    if node not in tax_types:                                       #check if node in nodes.dmp
        print('[Warning] Could not find {} in nodes.dmp while parsing taxonomy hierarchy\n'.format(node))

    else:
        children.append(node)
        for i in children:
            if tax_types[i] == ranking:
                tax_descendants.append(i)
            else:
                if i in tax_nodes:
                    for j in tax_nodes[i]:
                        children.append(j)

    return tax_descendants

def getTaxParent(tax_nodes, tax_types, taxid, ranking):

    '''
    input:
    - dictionary of form {parent: node} (readNodes output)
    - dictionary of form {node: type} (readNodes output)
    - taxid
    output:
    - dictionary of form {tax_id: [descendants]}
    '''
    tax_parents = {}
    node = str(taxid)
    if node not in tax_types:                                       #check if node in nodes.dmp
        tax_parents[node] = None
        print('[Warning] Could not find {} in nodes.dmp while parsing taxonomy hierarchy\n'.format(node))

    else:
        parent = tax_nodes[node]                                    #get parent for current node
        tax_parents[node] = [parent]                                #add node to dictionary

        while parent != tax_nodes[parent] and tax_types[parent]!= ranking:    #stop when parent = parent(parent) (i.e. 1 = 1)
            parent = tax_nodes[parent]                              #get parent of parent
            tax_parents[node].append(parent)                        #add parent to node in dictionary

    return tax_parents

# determine the lineage where your tax id belongs to (lineage taken until upper level = args.type)
taxparents,taxtypes=readNodes(args.nodesfile)
taxnames,namestax=readNames(args.namesfile)

eukgens=[]
prokgens=[]

spoifamily=""
spoiclade=""
spoigenus=args.spoi.split()[0]
if spoigenus in namestax:
    lineage=getTaxParent(taxparents,taxtypes,namestax[spoigenus],args.type)
    lineage2=getTaxParent(taxparents,taxtypes,namestax[spoigenus],'order')
    if lineage[namestax[spoigenus]] != None:
        spoifamily=taxnames[lineage[namestax[spoigenus]][-1]]
        spoiclade=taxnames[lineage2[namestax[spoigenus]][-1]]

print(spoifamily+'\t'+spoiclade)

k=open(args.tax,'r')
for line in k:
    sciname=line.split(';')[-2]
    if 'environmental' in sciname:
        sciname=line.split(';')[-3]
    print(sciname)
    if sciname in namestax:
        #print(sciname)
        lineage=getTaxParent(taxparents,taxtypes,namestax[sciname],args.type)
        fulllineage=getTaxParent(taxparents,taxtypes,namestax[sciname],'superkingdom')
        cladelineage=getTaxParent(taxparents,taxtypes,namestax[sciname],'order')
        rootlevelname=taxnames[fulllineage[namestax[sciname]][-1]]
        cladelevelname=taxnames[cladelineage[namestax[sciname]][-1]]
        if taxtypes[namestax[sciname]] == args.type:
            print('FAMILY:'+sciname+' CLADE:'+cladelevelname)
            taxlevelname=sciname
            if 'Eukaryota' == rootlevelname and taxlevelname not in eukgens and taxlevelname != spoifamily and cladelevelname != spoiclade:
                eukgens.append(taxlevelname)
            elif taxlevelname not in prokgens and taxlevelname != spoifamily and cladelevelname != spoiclade:
                prokgens.append(taxlevelname)
        elif lineage[namestax[sciname]] != None:
            print('LOWER THAN FAMILY:'+sciname)
            taxlevelname=taxnames[lineage[namestax[sciname]][-1]]
            #print(taxlevelname)
            #print(fulllineage)
            if int(lineage[namestax[sciname]][-1]) != 1:
            #if taxtypes[namestax[sciname]] == args.type:
                if 'Eukaryota' == rootlevelname and taxlevelname not in eukgens and taxlevelname != spoifamily and cladelevelname != spoiclade:
                    eukgens.append(taxlevelname)
                elif taxlevelname not in prokgens and taxlevelname != spoifamily and cladelevelname != spoiclade:
                    prokgens.append(taxlevelname)

file1=args.outdir+"/prok."+args.suffix
k=open(file1,'w')
for elem in prokgens:
    k.write(elem+'\n')
k.close()

file1=args.outdir+"/euk."+args.suffix
k=open(file1,'w')
for elem in eukgens:
    k.write(elem+'\n')
k.close()
'''
allchildren=[]
k=open(args.tax,'r')
for line in k:
    sciname=line.split()[0]
    #print(sciname)
    descendants=getTaxChildren(taxparents,taxtypes,namestax[sciname],args.type)
    allchildren.extend(descendants)
set_allchildren = set(allchildren)
uniqchildren=list(set_allchildren)

for child in uniqchildren:
    print(taxnames[child])
'''
