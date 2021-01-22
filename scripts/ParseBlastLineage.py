from __future__ import division
import argparse
import configparser
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-b", type=str, action='store', dest='blast', metavar='BLAST',help='define blast outfmt 6')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument("-na", type=str, action='store', dest='namesfile', metavar='NAMES',help='NCBI names.dmp')
parser.add_argument("-no", type=str, action='store', dest='nodesfile', metavar='NODES',help='NCBI nodes.dmp')
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
            if 'scientific' in line:
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
    return tax_nodes,tax_types


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


taxparents,taxtypes=readNodes(args.nodesfile)
taxnames,namestax=readNames(args.namesfile)

currentread=""
i=0
familyhits={}
f =open(args.blast,'r')
for record in f:
        record=record.strip()
        readname=record.split('\t')[0]
        if readname not  in familyhits:
            familyhits[readname]={}
            familyhits[readname]['families']={}
            familyhits[readname]['total']=0
        if readname != currentread:
            currentread=readname
            i=0
        elif i < 10:
            if float(record.split('\t')[2]) > 90:
                taxid=record.split('\t')[1].split('|')[-1]
                lineage=getTaxParent(taxparents,taxtypes,taxid,'family')
                hitfamily=taxnames[lineage[taxid][-1]]
                if hitfamily != 'root':
                    familyhits[readname]['total']=familyhits[readname]['total']+1
                    if hitfamily not in familyhits[readname]['families']:
                        familyhits[readname]['families'][hitfamily] = 1
                    else:
                        familyhits[readname]['families'][hitfamily] = familyhits[readname]['families'][hitfamily] +1
                    #print(readname+'\t'+record.split('\t')[1]+'\t'+taxid+'\t'+hitfamily)
        i=i+1


finalfams=[]
for read in familyhits:
    numberhits=int(int(familyhits[read]['total'])/2)
    #print(read+'\t'+str(numberhits))
    for fam in familyhits[read]['families']:
        #print(read+'\t'+str(familyhits[read]['families'][fam])+'\t'+fam)
        if int(familyhits[read]['families'][fam]) > numberhits:
            #print(read+'\t'+str(familyhits[read]['families'][fam])+'\t'+fam)
            if fam not in finalfams:
                finalfams.append(fam)
                fulllineage=getTaxParent(taxparents,taxtypes,namestax[fam],'superkingdom')
                lineagetoconv=fulllineage[namestax[fam]]
                lineagestring=""
                for elem in reversed(lineagetoconv):
                    lineagestring=lineagestring+taxnames[elem]+";"
                lineagestring=lineagestring+fam+';'
                print(lineagestring)
