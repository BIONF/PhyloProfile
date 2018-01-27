# Demo data
*PhyloProfile* comes with some test data that you can use to fully explore the functionality.

Detail about input files are expained in [Phyloprofile Wiki Page](https://github.com/BIONF/PhyloProfile/wiki/Input-Data).

## Main input files:
- `test.main.wide`: input in wide (matrix) format.
- `test.main.long`: input in long format.
- `test.main.fasta`: input in fasta format.
- `test.main.xml`: input in OrthoXML format.

Use one of those files as the **Main input** file on the *Input & settings* page after starting *PhyloProfile*.

## Optional input files:
- `domain_files/*.domains`: This folder contains the feature architecture data (e.g. Pfam domains) that you can optionally give under the **Additional annotation input** upload on the *Input & settings* page after startup.
- `fasta_files/*.fa`: This folder contains the fasta sequences for demo input. In order to show the sequence in **Detailed plot** you have to set path to fasta folder in **FASTA config** on the *Input & settings* page after startup.
- `other/test.geneList`: After doing the initial plot with the files above you can use this file on the *Customized profile* tab to sub-select for only the genes present in this file.
- `other/test.taxaList`: This contains list of taxon names. Use this to test the function of fetching NCBI taxonomy IDs, which can be found in *Function* tab.

## Other pre-processing input files: for these files the user needs to run parsing scripts to prepare the compatible inputs for PhyloProfile.
- `oma/oma_example.orthoxml`: this is an output file from OMA standalone. Although it has OrthoXML format, but it lacks the taxonomy information. In order to convert it to an input file for PhyloProfile, please run this command

	`python scripts/convert_oma_standalone_orthoxml.py -x oma_example.orthoxml -m taxon_mapping_oma_orthoxml.csv > oma_example_phyloprofile_compatible.orthoxml`

- oma/omaIDs.list: use this file for testing fetching OMA HOGs by using the command

	`python ./scripts/get_oma_hogs.py -i omaIDs.list`

- `pfamAnno/hmmscanOut.txt` and `pfamAnno/pfamscanOut.txt`: PFAM annotation output files from hmmscan and pfamscan. Run the following commands to convert them into compatible domain files:

	`python scripts/hmmscanParser.py -i test.input.hmmscanOut > test.input.domains`

	(replace *hmmscanParser.py* by *pfamscanParser.py* for parsing pfamscan output file).
