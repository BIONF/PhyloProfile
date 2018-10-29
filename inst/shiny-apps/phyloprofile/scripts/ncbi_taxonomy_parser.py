# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Parse info from taxonomy NCBI dmp files (names.dmp and nodes.dmp)
# to create a mapping file taxonNamesFull.txt that contains
# ncbiID	fullName	rank	parentID
# Author: Carla Mölbert {carla.moelbert@gmx.de}
# Date: 10.07.2018

import argparse
import subprocess

# Get the arguments for the script
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
					help = 'folder contain the names.dmp and nodes.dmp files(required)',
					required = True)
args = parser.parse_args()

# Generate dictionary with key = ncbiID and value = fullName
def get_names_dict(names):
	names_dict = dict()
	for name in names:
		if ("scientific name" in name):
			tmp = name.split("\t|\t")
			names_dict[tmp[0]] = tmp[1]
			#names_dict = names_dict.update({name[0]: name[1]})
	return(names_dict)

# Get a new line for the output file
def get_new_line (node, names_dict):
	tmp = node.split("\t|\t")
	ncbiID = tmp[0]
	parentID = tmp[1]
	rank = tmp[2]
	fullName = names_dict[ncbiID]
	new_line = new_line = ncbiID + "\t" + fullName + "\t" + rank + "\t" + parentID + "\n"
	return(new_line)

########################################
output = "taxonNamesFull.txt"
file = open(output,"w")
file.write("ncbiID\tfullName\trank\tparentID\n") # header

# names to open the files
filename_names = args.input + "/names.dmp"
filenames_nodes = args.input + "/nodes.dmp"

# Get the needed information from names.dmp
with open(filename_names) as names_file:
	print("Read names.dmp \n")
	names = names_file.readlines()
names_dict = get_names_dict(names)

# Combine the information with the information from nodes.dmp
with open (filenames_nodes) as nodes_file:
	print("Read nodes.dmp \n")
	nodes = nodes_file.read().splitlines()

print("Write taxonNamesFull.txt")
for node in nodes:
	new_line = get_new_line(node, names_dict)
	file.write(new_line)

# replace "no rank" by "norank"
cmd = "sed -i -e 's/no rank/norank/g' ./taxonNamesFull.txt"
subprocess.call(cmd, shell = True)

print("Finished! Please replace phyloprofile/data/taxonNamesFull.txt by the new taxonNamesFull.txt file!")
