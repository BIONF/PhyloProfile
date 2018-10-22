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
		print('hmmscanParser.py -i hmmscanOutput.txt')
		sys.exit(2)

	for opt,arg in opts:
		if opt in ('-h','--help'):
			print('hmmscanParser.py -i <hmmscan output file>')
			sys.exit()
		elif opt in ('-i','--inFile'):
			inFile = arg

	with open(inFile) as fp:
   		for line in fp:
			if not re.match("^#", line) is not None:
				hit = line.strip("\n").split()
				geneID = hit[3]
				print(geneID.split('|')[0]+"#"+geneID.split('|')[2]+"\t"+geneID.split('|')[2]+"\t"+hit[5]+"\t"+"pfam_"+hit[0]+"\t"+hit[19]+"\t"+hit[20]+"\tNA"+"\tN")
				# time.sleep(1)

if __name__ == "__main__":
	if len(sys.argv[1:])==0:
		print('hmmscanParser.py -i hmmscanOutput.txt')
		sys.exit(2)
	else:
		main(sys.argv[1:])
