# -*- coding: utf-8 -*-
import sys
import getopt
import glob
import time
import re

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
		print('fastaParser.py -i input.fasta')
		sys.exit(2)

	for opt,arg in opts:
		if opt in ('-h','--help'):
			print('fastaParser.py -i <fasta file>')
			sys.exit()
		elif opt in ('-i','--inFile'):
			inFile = arg

	print ("geneID\tncbiID\torthoID\tvar1\tvar2")
	with open(inFile) as fp:
   		for line in fp:
			if re.match(">", line) is not None:
				hit = line.rstrip().replace(">","").split('|')
				strOut = '\t'.join(hit).strip()
				# hit = hit.replace('|', "\t")
				if len(hit) < 5:
					for i in range(len(hit),5):
						strOut = strOut +"\t"+"NA"

				print(strOut)

if __name__ == "__main__":
	if len(sys.argv[1:])==0:
		print('fastaParser.py -i input.fasta')
		sys.exit(2)
	else:
		main(sys.argv[1:])
