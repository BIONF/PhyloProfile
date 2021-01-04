#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-x", "--xml", required=True,
                    help="OrthoXML file that should be read (required)")
parser.add_argument("-m", "--mapping", required=True,
                    help="tab-separated file mapping OMA standalone\
                    input filenames to NCBI Taxonomy IDs (required)")
args = parser.parse_args()

def read_mapping(mapping_file):
    '''
    read the tab-separated file
    1. column: filenames used for OMA (w/ or w/o .fa-suffix)
    2. column: NCBI taxonomy-ID
    '''
    taxon_ids = {}
    for i in open(mapping_file):
        taxon_name, ncbi_id = i.strip().split("\t")
        if taxon_name[-3:] == ".fa":
            taxon_name = taxon_name[:-3]
        taxon_ids[taxon_name] = ncbi_id
    return taxon_ids

def convert_orthoxml(xml_file,name_to_ncbi):
    for line in open(xml_file):
        line = line
        if "xmlns=" in line:
            lineTmp = line.replace('xmlns=', 'xmlns:xsi=')
            print(lineTmp.rstrip())
        elif "NCBITaxId=" in line:
            try:
                taxon_name = line.split("name=")[1].split(" NCBI")[0].replace('"','')
                ncbi_id = name_to_ncbi[taxon_name]
                print(' <species name="'+taxon_name+'" NCBITaxId="' + ncbi_id +'">')
            except:
                print("------------------------")
                print("Could not find %s in the mapping file" % (taxon_name))
                print("Please fix the mapping file and try again")
                exit()
        else:
            print(line.rstrip())


name_to_ncbi = read_mapping(args.mapping)
convert_orthoxml(args.xml,name_to_ncbi)
