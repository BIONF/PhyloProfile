# -*- coding: utf-8 -*-
import sys
import getopt
import glob
import time
from bs4 import BeautifulSoup

def appendToDict(key,value,aDict):
	if not key in aDict.keys():
		aDict[key] = []
		aDict[key].append(value)
	else:
		if not value in aDict[key]:
			aDict[key].append(value)

def main(argv):
	inFile = ''
	try:
		opts, args = getopt.getopt(argv,"i:h",["inFile","help"])
	except getopt.GetoptError:
		print('orthoxmlParser.py -i input')
		sys.exit(2)

	for opt,arg in opts:
		if opt in ('-h','--help'):
			print('orthoxmlParser.py -i <orthoxml file>')
			sys.exit()
		elif opt in ('-i','--inFile'):
			inFile = arg

	##### read file into beatifulsoup object
	xmlIn = BeautifulSoup(open(inFile),"xml")

	##### PARSING XML FILE
	### get list of species together with NCBI taxon IDs and their corresponding genes
	taxonID = {}	# taxonID{(protID:taxonID)}
	protID = {}	# protID{(geneID:protID)}

	for spec in xmlIn.findAll("species"):
		taxID = spec.get("NCBITaxId")	# taxon ID
		for gene in spec.findAll("gene"):		# <gene id="1" protid="NP_050056.1|ATP synthase F0 subunit 8"></gene>
			geneID = gene.get("id")
			orthoID = gene.get("protId").split('|')[0]		# protid="NP_050057.1|NADH dehydrogenase subunit 5"

			taxonID[orthoID] = taxID
			protID[geneID] = orthoID

	### get score IDs and theis descriptions
	scoreDesc = {}	# scoreDesc{scoreID:scoreDesc}
	header = "geneID\tncbiID\torthoID"
	for score in xmlIn.findAll("scoreDef"):		# <scoredef desc="Distance between edge seed ortholog" id="inparalog"></scoredef>
		scoreDesc[score.get("id")] = score.get("desc")
		header = header+"\t"+score.get("id")
	print(header)

	### get ortholog info for each group (groupID ncbiID orthoID scores)
	index = 10000
	for orthogroup in xmlIn.findAll("orthologGroup"):
		groupID = orthogroup.get("id")
		if groupID.isdigit():
			groupIndex = index + int(groupID)
			groupID = "OG_"+str(groupIndex)

		for ortho in orthogroup.findAll("geneRef"):
			orthoID = protID[ortho.get("id")]
			result = [groupID,"ncbi"+taxonID[orthoID],orthoID]

			for scoreID in scoreDesc.keys():
				if(ortho.findAll("score", id = scoreID)):
					score = ortho.find("score", id = scoreID)
					result.append(score.get("value"))
				else:
					result.append("NA")

			### PRINT OUTPUT
			print('\t'.join(result))

if __name__ == "__main__":
	if len(sys.argv[1:])==0:
		print('orthoxmlParser.py -i input')
		sys.exit(2)
	else:
		main(sys.argv[1:])
