# -*- coding: utf-8 -*-
#!/usr/bin/env python
import coreapi
import sys
import argparse
import urllib
from cStringIO import StringIO
from lxml import etree
from bs4 import BeautifulSoup
import time
import re

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',
                    help='list of one or more OMA protein-IDs that will be downloaded.\
                    Multiple IDs are separated by newline.\
                    (required)',required=True)
parser.add_argument('-t','--type',
                    help='type of orthologs: PAIR, OG or HOG.\
                    (required)',required=True)
args = parser.parse_args()

# Initialize a client & load the OMA schema document
client = coreapi.Client()
schema = client.get("https://omabrowser.org/api/docs")

# create orthoxml profile file
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

def get_seq_info(ID):
	# Interact with the API endpoint
	action = ["protein", "read"]
	params = {
	    "entry_id": ID
	}
	result = client.action(schema, action, params=params)
	return(result)

def get_domain_info(ID):
	action = ["protein", "domains"]
	params = {
	    "entry_id": ID
	}
	result = client.action(schema, action, params=params)
	return(result)

def get_orthogroup_info(ID):
	action = ["group", "read"]
	params = {
		"group_id": ID
	}
	result = client.action(schema, action, params=params)
	return(result)

def get_orthopair_info(ID):
	action = ["protein", "orthologs"]
	params = {
		"entry_id": ID
	}
	result = client.action(schema, action, params=params)
	return(result)

def get_taxonomy_info(ID):
	action = ["taxonomy", "read"]
	params = {
		"root_id": ID		# either the taxon id, species name or the 5 letter UniProt species code for a root taxonomic level
	}
	result = client.action(schema, action, params=params)
	return(result)

def print_domain(domainInfo,pairID,geneID):
	domainBlock = ""
	for domain in domainInfo["regions"]:
		source = domain["source"]
		name = domain["name"]
		region = domain["location"].split(':')
		start = region[0]
		end = region[1]
		if(len(name) > 0):
			domainLine = pairID+"\t"+geneID+"\t"+source+" "+name+"\t"+start+"\t"+end+"\n"
			domainBlock += domainLine
	return(domainBlock)

###########################
output = args.input+".phyloprofile"
file = open(output,"w")
# outputLong = args.input+".phyloprofile"
# outputLongFile = open(outputLong,"w")
fasOut = args.input+".fa"
fasOutFile = open(fasOut,"w")
domainOut = args.input+".domains"
domainOutFile = open(domainOut,"w")

### get HOGs
if args.type == "HOG":
	# get first file in list and read it
	with open(args.input) as f:
	    ids = f.read().splitlines()

	merged_id = ids[0]
	merged_xml = read_xml(merged_id)

	# now read all other files and merge
	other_ids = ids[1:]
	for single_id in other_ids:
	    xml_to_append = read_xml(single_id)
	    merged_xml = merge_xml(merged_xml,xml_to_append)

	# and print our new XML
	file.write(etree.tostring(merged_xml))

	# get sequences and domain annotations
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

	print("now getting FASTA sequence and domain annotations for:")
	for orthogroup in soup.findAll("orthologGroup"):
		groupID = orthogroup.get("id")
		if groupID:
			if groupID.isdigit():
				groupID = "OG_"+str(groupID)
			for ortho in orthogroup.findAll("geneRef"):
				geneID = ortho.get("id")
				print(geneid2name[geneID])

				# output sequence
				seq = get_seq_info(geneid2name[geneID])["sequence"]
				fasta = ">"+groupID+"|"+gene2spec[geneID]+"|"+geneid2name[geneID]+"\n"+seq+"\n"
				fasOutFile.write(fasta)

				# output domains
				domainInfo = get_domain_info(geneid2name[geneID])
				for domain in domainInfo["regions"]:
					source = domain["source"]
					name = domain["name"]
					region = domain["location"].split(':')
					start = region[0]
					end = region[1]
					if(len(name) > 0):
						domainLine = groupID+"#"+geneid2name[geneID]+"\t"+geneid2name[geneID]+"\t"+source+" "+name+"\t"+start+"\t"+end+"\n"
						domainOutFile.write(domainLine)

elif (args.type == "OG") or (args.type == "PAIR"):
	file.write("geneID\tncbiID\torthoID\n")
	with open(args.input) as f:
	    ids = f.read().splitlines()
	for id in ids:
		# get info for seedID
		seedInfo = get_seq_info(id)
		specID = get_taxonomy_info(seedInfo["omaid"][:5])["id"]
		profileSeed = "OG_%s\tncbi%s\t%s\n" % (id,specID,id)
		file.write(profileSeed)
		headerSeed = ">OG_%s|ncbi%s|%s\n" % (id,specID,id)
		fastaSeed = headerSeed+seedInfo["sequence"]+"\n"
		fasOutFile.write(fastaSeed)
		domainSeed = get_domain_info(id)
		domainSeedBlock = print_domain(domainSeed,"OG_"+str(id)+"#"+str(id),id)
		if len(domainSeedBlock) == 0:
			domainSeedBlock = "OG_"+str(id)+"#"+str(id)+"\t"+id+"\tn/a\t0\t0\n"
		domainOutFile.write(domainSeedBlock)

		# get group ID
		orthoGroupID = get_seq_info(id)["oma_group"]
		print(orthoGroupID)

		# get protein members
		orthoGroup = get_orthopair_info(id)
		if(args.type == "OG"):
			orthoGroup = get_orthogroup_info(orthoGroupID)["members"]

		for mem in orthoGroup:
			geneID = mem["omaid"]
			print(geneID)
			specID = get_taxonomy_info(geneID[:5])["id"]

			profileLine = "OG_%s\tncbi%s\t%s\n" % (id,specID,geneID)
			file.write(profileLine)

			geneInfo = get_seq_info(geneID)
			header = ">OG_%s|ncbi%s|%s\n" % (id,specID,geneID)
			fasta = header+geneInfo["sequence"]+"\n"
			fasOutFile.write(fasta)

			domainInfo = get_domain_info(geneID)
			domainBlock = print_domain(domainInfo,"OG_"+str(id)+"#"+str(geneID),geneID)
			if len(domainBlock) == 0:
				domainBlock = "OG_"+str(id)+"#"+str(geneID)+"\t"+geneID+"\tn/a\t0\t0\n"
			domainSeedBlock = print_domain(domainSeed,"OG_"+str(id)+"#"+str(geneID),id)
			if len(domainSeedBlock) == 0:
				domainSeedBlock = "OG_"+str(id)+"#"+str(geneID)+"\t"+id+"\tn/a\t0\t0\n"
			domainOutFile.write(domainBlock)
			domainOutFile.write(domainSeedBlock)

else:
	parser.print_help()
	print("ERROR: un-accepted TYPE given!")
	sys.exit(0)

print("FINISHED!\n")
print("(1/3) Profile input file: "+output+"\n")
print("(2/3) Fasta file: "+fasOut+"\n")
print("(3/3) Domain file: "+domainOut+"\n")
