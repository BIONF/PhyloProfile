# Demo data
*PhyloProfile* comes with some small test data that you can use to fully explore the functionality.

Detail about input files are expained in [Phyloprofile Wiki Page](https://github.com/BIONF/PhyloProfile/wiki/Input-Data).

## Description
This test data is a *subset* of the [LCA Microsporidia data](https://github.com/BIONF/phyloprofile-data/blob/master/expTestData/lca_microsporidia.zip) used for [testing the performance](https://github.com/BIONF/PhyloProfile/wiki/Performance-test) of *PhyloProfile*. It contains the phylogenetic profiles of 10 genes across 10 taxa, integrated with 2 additional information layers: **(1) Domain architecture similarity** and **(2) Traceability** scores.

*Domain architecture similarity* scores are used to compare the protein architecture between seed and its ortholog ([Koestler et al. (2010), BMC Bioinformatics](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2931517/)). While *Traceability* of a protein defines the point beyond which sequence similarity based approaches are bound to fail for a ortholog prediction (Jain et al., unpublished).

## Content
### Main input files:
contain phylogenetic profiles integrated with domain similarity and traceability scores.
- `test.main.wide`: input in wide (matrix) format.
- `test.main.long`: input in long format.
- `test.main.fasta`: input in fasta format.
- `test.main.xml`: input in OrthoXML format.

Use one of those files as the **Main input** file on the *Input & settings* page after starting *PhyloProfile*.

### Optional input files:
- `domain_files/*.domains`: This folder contains the feature architecture data that you can optionally give under the **Additional annotation input** upload on the *Input & settings* page after startup.
- `fasta_files/*.fa`: This folder contains the fasta sequences for demo input. It has separate multi-fasta files for input taxa (e.g. `A.thaliana@3702.fa` or `H.sapiens@9606.fa`) and a concatenated file (`concatenatedFile.fa`) contains all sequences. In order to show the sequence in **Detailed plot** you have 2 options, either (1) upload the concatenated file using the option *Concatenated fasta file* or (2) set path to fasta folder with the option *Fasta folder* in **FASTA config** on the *Input & settings* page after startup.
- `other/test.geneList`: After doing the initial plot with the files above you can use this file on the *Customized profile* tab to sub-select for only the genes present in this file.
- `other/test.taxaList`: This contains list of taxon names. Use this to test the function of fetching NCBI taxonomy IDs, which can be found in *Function* tab.

### Other pre-processing input files: for these files the user needs to run parsing scripts to prepare the compatible inputs for PhyloProfile.
- `oma/oma_example.orthoxml`: this is an output file from OMA standalone. Although it has OrthoXML format, but it lacks the taxonomy information. In order to convert it to an input file for PhyloProfile, please run this command

	`python scripts/convert_oma_standalone_orthoxml.py -x oma_example.orthoxml -m taxon_mapping_oma_orthoxml.csv > oma_example_phyloprofile_compatible.orthoxml`

- `oma/omaIDs.list`: use this file for testing fetching OMA orthologs (OMA pairs, OMA groups or OMA HOGs) by using the command

	`python ./scripts/get_oma_browser.py -i omaIDs.list -t PAIR`

	*see [this WIKI Section](https://github.com/BIONF/PhyloProfile/wiki/Input-Data#oma-browser) for more detail about `get_oma_browser.py`.*

- `pfamAnno/hmmscanOut.txt` and `pfamAnno/pfamscanOut.txt`: PFAM annotation output files from hmmscan and pfamscan. Run the following commands to convert them into compatible domain files:

	`python scripts/hmmscanParser.py -i test.input.hmmscanOut > test.input.domains`

	(replace *hmmscanParser.py* by *pfamscanParser.py* for parsing pfamscan output file).

### Detail about input files are expained in our [PhyloProfile Wiki Page](https://github.com/BIONF/PhyloProfile/wiki/Input-Data)
