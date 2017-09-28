#!/usr/bin/env python
import argparse
import urllib
from cStringIO import StringIO
from lxml import etree


parser = argparse.ArgumentParser()
parser.add_argument('-i','--ids', nargs='+',
                    help='one or more OMA protein-IDs that will be downloaded.\
                    Multiple IDs are separated by spaces.\
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
merged_id = args.ids[0]
merged_xml = read_xml(merged_id)

# now read all other files and merge
other_ids = args.ids[1:]
for single_id in other_ids:
    xml_to_append = read_xml(single_id)
    merged_xml = merge_xml(merged_xml,xml_to_append)

# and print our new XML
print(etree.tostring(merged_xml))
