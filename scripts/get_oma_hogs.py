# -*- coding: utf-8 -*-
#!/usr/bin/env python
from __future__ import print_function
import argparse
import urllib
from cStringIO import StringIO
from lxml import etree
from bs4 import BeautifulSoup
import time
import re

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', nargs='+',
                    help='list of one or more OMA protein-IDs that will be downloaded.\
                    Multiple IDs are separated by newline.\
                    (required)',required=True)
args = parser.parse_args()

def read_xml(ID):
    '''
    Get the XML-formatted HOG file from OMA.
    Parse as XML-tree for downstream merging if needed
    '''
    try:
        oma_url = urllib.urlopen("http://omabrowser.org/oma/hogs/%s/orthoxml/" % ID)
        oma_xml = oma_url.read()
        oma_xml_tree = etree.parse(StringIO(oma_xml))
        return oma_xml_tree
    except:
        print("------------------------")
        print("Could not find %s in OMA" % (ID))
        print("Are you sure that this is a valid protein ID?")
        exit()

def merge_xml(oma1,oma2):
    '''
    Merge two OMA OrthoXML files by adding
    a) missing species with the corresponding genes
    b) the HOG definitions that are not in <groups> yet
    '''
    oma1_root = oma1.getroot()
    oma2_root = oma2.getroot()
    # get the namespace for the XML
    ns = "{"+oma1_root.nsmap[None]+"}"
    # iterate over all species in file #2
    for species in oma2.iter(ns+"species"):
        # get taxonomy id
        tax_id_oma2 = species.attrib["NCBITaxId"]
        # try finding taxon in file #1
        taxon_oma1 = oma1_root.find("{http://orthoXML.org/2011/}species[@NCBITaxId='"+tax_id_oma2+"']")
        # file #1 has the taxon from file #2, append all genes
        if taxon_oma1 != None:
            taxon_oma1_genes = taxon_oma1.find(ns+"database").find(ns+"genes")
            taxon_oma2_genes = species.find(ns+"database").find(ns+"genes")
            for i in taxon_oma2_genes.getchildren():
                taxon_oma1_genes.append(i)
        # meh, have to add the whole species to it.
        else:
            oma1_root.findall(ns+"species")[-1].addnext(species)
    # now add the orthologous group definition to the corresponding block
    for group in oma2.find(ns+"groups").getchildren():
        oma1.find(ns+"groups").getchildren()[-1].addnext(group)
    return oma1

# get first file in list and read it
with open(args.input[0]) as f:
    ids = f.read().splitlines()

merged_id = ids[0]
merged_xml = read_xml(merged_id)

# now read all other files and merge
other_ids = ids[1:]
for single_id in other_ids:
    xml_to_append = read_xml(single_id)
    merged_xml = merge_xml(merged_xml,xml_to_append)

# and print our new XML
output = args.input[0]+".orthoXML"
file = open(output,"w")
file.write(etree.tostring(merged_xml))
print("(1/2) Profile input file is created: "+output+"\n")
time.sleep(3)

# get fasta sequences
def get_fasta(ID):
	oma_url = urllib.urlopen("https://omabrowser.org/cgi-bin/gateway.pl?f=DisplayEntry&p1=%s" % ID)
	page = oma_url.read()
	oma_gene_soup = BeautifulSoup(page, 'html.parser')
	fasta_section = oma_gene_soup.find_all('div', attrs={'class': 'panel-body'})[-1].text.strip()
	seq = ''.join(''.join(i for i in fasta_section if not i.isdigit()).split())
	# fasta = ">"+ID+"\n"+seq+"\n"
	return(seq)

soup = BeautifulSoup(open(output),"xml")
gene2spec = {}
geneid2name = {}
for spec in soup.findAll("species"):
	specID = spec.get("NCBITaxId")
	for gene in spec.findAll("gene"):
		geneName = gene.get("protId")
		geneID = gene.get("id")
		gene2spec[geneID] = "ncbi"+specID
		geneid2name[geneID] = geneName

fasOut = args.input[0]+".fa"
fasOutFile = open(fasOut,"w")
print("now getting FASTA sequences for:")
for orthogroup in soup.findAll("orthologGroup"):
	groupID = orthogroup.get("id")
	if groupID:
		if groupID.isdigit():
			groupID = "OG_"+str(groupID)
		for ortho in orthogroup.findAll("geneRef"):
			geneID = ortho.get("id")
			print(geneid2name[geneID])
			seq = get_fasta(geneid2name[geneID])
			fasta = ">"+groupID+"|"+gene2spec[geneID]+"|"+geneid2name[geneID]+"\n"+seq+"\n"
			fasOutFile.write(fasta)
print("(2/2) Fasta file is created: "+fasOut+"\n")
